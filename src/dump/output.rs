use std::io::Write;
use std::sync::mpsc;
use std::sync::Arc;
use std::thread;
use std::time::Duration;

use anyhow::Result;
use parking_lot::Condvar;
use parking_lot::Mutex;

use crate::{
    cli::{FilterOptions, OutputFormat},
    output::{build_writers, Compression},
    BUFFER_SIZE,
};

/// Set the default buffer size to 1MB
const DEFAULT_BUFFER_SIZE: usize = 1024 * 1024;

/// Set the maximum overflow buffer size to 128MB
const MAXIMUM_BUFFER_SIZE: usize = 128 * 1024 * 1024;

/// Sets the default sleep time (in milliseconds)
const DEFAULT_SLEEP_MS: u64 = 100;

/// A shorthand for the type of output handles we expect to write
pub type BoxedWriter = Box<dyn Write + Send>;

/// A shorthand for the type of writers we expect to use
pub type BoxedSegmentWriter = Box<dyn SegmentWriter + Send>;

/// Reusable trait for Writer structs which handle IO of segments as a group
pub trait SegmentWriter {
    /// Number of segments expected by the writer
    fn num_segments(&self) -> usize;

    /// Write all the segments to their respective IO handles
    fn write_all_buffers(&mut self, buffers: &mut [Vec<u8>], counts: &mut [usize]) -> Result<()>;

    /// Return local buffers to mimic the expected writer buffers on-thread
    fn generate_local_buffers(&self) -> Vec<Vec<u8>> {
        vec![Vec::with_capacity(BUFFER_SIZE); self.num_segments()]
    }
}

/// Handles the creation logic and pipes the IO to the right Writer struct.
#[allow(clippy::too_many_arguments)]
pub fn build_segment_writer(
    outdir: Option<&str>,
    prefix: &str,
    compression: Compression,
    format: OutputFormat,
    num_threads: usize,
    filter_opts: &FilterOptions,
    is_fifo: bool,
    is_split: bool,
) -> Result<BoxedSegmentWriter> {
    if is_split {
        if is_fifo {
            let wtr = BufferedWriter::new(
                outdir,
                prefix,
                compression,
                format,
                num_threads,
                filter_opts,
                is_fifo,
            )?;
            Ok(Box::new(wtr))
        } else {
            let wtr = DirectWriter::new(
                outdir,
                prefix,
                compression,
                format,
                num_threads,
                filter_opts,
                is_fifo,
            )?;
            Ok(Box::new(wtr))
        }
    } else {
        let wtr = DirectWriter::new(
            None,
            prefix,
            compression,
            format,
            num_threads,
            filter_opts,
            false,
        )?;
        Ok(Box::new(wtr))
    }
}

/// A thead-local writer that owns a subprocess handling the actual writing
struct ThreadWriter {
    /// Owned reusable write buffer with a conditional variable marking when it's been written to
    buffer_pair: Arc<(Mutex<Vec<u8>>, Condvar)>,
    /// The signal to end the subprocess
    shutdown_sender: mpsc::Sender<()>,
    /// Handle to the owned subprocess
    join_handle: Option<thread::JoinHandle<Result<()>>>,
}

impl ThreadWriter {
    fn new(mut handle: BoxedWriter) -> Self {
        let buffer_pair = Arc::new((Mutex::new(Vec::new()), Condvar::new()));
        let buffer_pair_clone = Arc::clone(&buffer_pair);

        let (shutdown_sender, shutdown_receiver) = mpsc::channel();

        // Start the worker thread
        let join_handle = thread::spawn(move || -> Result<()> {
            let (buffer, cvar) = &*buffer_pair_clone;

            loop {
                // Wait for data or shutdown signal
                let mut guard = buffer.lock();

                while guard.is_empty() {
                    // Check for shutdown before waiting
                    if shutdown_receiver.try_recv().is_ok() {
                        return Ok(()); // Exit the thread
                    }

                    // Wait for notification that data is available
                    cvar.wait(&mut guard);

                    // Check again for shutdown after waking
                    if shutdown_receiver.try_recv().is_ok() {
                        return Ok(()); // Exit the thread
                    }
                }

                // We have data to process
                let mut data = std::mem::take(&mut *guard);
                drop(guard); // Release lock before I/O

                // Perform actual write (potentially blocking I/O)
                handle.write_all(data.drain(..).as_slice())?;
                handle.flush()?;
            }
        });

        ThreadWriter {
            buffer_pair,
            shutdown_sender,
            join_handle: Some(join_handle),
        }
    }

    fn ingest(&self, data: &[u8]) {
        let (buffer, cvar) = &*self.buffer_pair;
        loop {
            let mut guard = buffer.lock();
            if guard.len() <= MAXIMUM_BUFFER_SIZE {
                guard.extend_from_slice(data);
                cvar.notify_one();
                break;
            } else {
                thread::sleep(Duration::from_millis(DEFAULT_SLEEP_MS));
            }
        }
    }
}

impl Drop for ThreadWriter {
    fn drop(&mut self) {
        // Signal thread to shut down
        self.shutdown_sender
            .send(())
            .expect("Error in sending signal");

        // notify the condition variable to wake up the worker thread
        let (_buffer, cvar) = &*self.buffer_pair;
        cvar.notify_all(); // make sure to wake the thread even if the buffer is empty

        // Wait for thread to finish
        if let Some(handle) = self.join_handle.take() {
            handle
                .join()
                .expect("Error in joining thread")
                .expect("Error within thread");
        }
    }
}

/// A writer struct which writes to output threads through an intermediary buffer.
///
/// The downstream handles are run as child processes so that IO is not blocked (for FIFO)
pub struct BufferedWriter {
    segment_buffers: Vec<Vec<u8>>,
    thread_writers: Vec<ThreadWriter>,
}
impl BufferedWriter {
    pub fn new(
        outdir: Option<&str>,
        prefix: &str,
        compression: Compression,
        format: OutputFormat,
        num_threads: usize,
        filter_opts: &FilterOptions,
        is_fifo: bool,
    ) -> Result<Self> {
        let segment_handles = build_writers(
            outdir,
            prefix,
            compression,
            format,
            num_threads,
            filter_opts,
            is_fifo,
        )?;
        let segment_buffers = vec![Vec::with_capacity(DEFAULT_BUFFER_SIZE); segment_handles.len()];
        let thread_writers = segment_handles.into_iter().map(ThreadWriter::new).collect();
        Ok(Self {
            segment_buffers,
            thread_writers,
        })
    }

    fn write_to_handles(&mut self) -> Result<()> {
        for (writer, buf) in self
            .thread_writers
            .iter()
            .zip(self.segment_buffers.iter_mut())
        {
            if !buf.is_empty() {
                writer.ingest(buf.drain(..).as_slice());
            }
        }
        Ok(())
    }
}
impl SegmentWriter for BufferedWriter {
    fn num_segments(&self) -> usize {
        self.thread_writers.len()
    }

    fn write_all_buffers(&mut self, buffers: &mut [Vec<u8>], counts: &mut [usize]) -> Result<()> {
        for (shared_buf, (local_buf, local_count)) in self
            .segment_buffers
            .iter_mut()
            .zip(buffers.iter_mut().zip(counts.iter_mut()))
        {
            // Skip writing empty segments
            if *local_count == 0 {
                continue;
            }
            shared_buf.extend_from_slice(local_buf);
            local_buf.clear();
            *local_count = 0;
        }

        self.write_to_handles()
    }
}

/// A Writer struct which writes directly to output handles without any buffering
pub struct DirectWriter {
    segment_handles: Vec<BoxedWriter>,
}

impl DirectWriter {
    pub fn new(
        outdir: Option<&str>,
        prefix: &str,
        compression: Compression,
        format: OutputFormat,
        num_threads: usize,
        filter_opts: &FilterOptions,
        is_fifo: bool,
    ) -> Result<Self> {
        let segment_handles = build_writers(
            outdir,
            prefix,
            compression,
            format,
            num_threads,
            filter_opts,
            is_fifo,
        )?;
        Ok(Self { segment_handles })
    }
}

impl SegmentWriter for DirectWriter {
    fn num_segments(&self) -> usize {
        self.segment_handles.len()
    }

    fn write_all_buffers(&mut self, buffers: &mut [Vec<u8>], counts: &mut [usize]) -> Result<()> {
        for (handle, (local_buf, local_count)) in self
            .segment_handles
            .iter_mut()
            .zip(buffers.iter_mut().zip(counts.iter_mut()))
        {
            // Skip writing empty segments
            if *local_count == 0 {
                continue;
            }
            handle.write_all(local_buf.drain(..).as_slice())?;
            handle.flush()?;
            *local_count = 0;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{self, Write};
    use std::sync::{Arc, Mutex};

    // Simple in-memory writer that lets us inspect written data
    struct TestWriter {
        data: Arc<Mutex<Vec<u8>>>,
    }

    impl Write for TestWriter {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            let mut guard = self.data.lock().unwrap();
            guard.extend_from_slice(buf);
            Ok(buf.len())
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    // DirectWriter::write_all_buffers tests
    #[test]
    fn direct_writer_write_all_buffers_happy_path_with_empty_segment() {
        // Shared data holders so we can inspect after the call
        let data1 = Arc::new(Mutex::new(Vec::new()));
        let data2 = Arc::new(Mutex::new(Vec::new()));

        let writer1: Box<dyn Write + Send> = Box::new(TestWriter {
            data: data1.clone(),
        });
        let writer2: Box<dyn Write + Send> = Box::new(TestWriter {
            data: data2.clone(),
        });

        let mut dw = DirectWriter {
            segment_handles: vec![writer1, writer2],
        };

        // Prepare buffers: first segment has data, second is empty
        let mut buffers = vec![b"ACGT".to_vec(), Vec::new()];
        let mut counts = vec![1, 0];

        dw.write_all_buffers(&mut buffers, &mut counts).unwrap();

        // After writing, buffers should be drained and counts reset
        assert!(buffers[0].is_empty());
        assert!(buffers[1].is_empty());
        assert_eq!(counts, vec![0, 0]);

        // Verify data written to the correct writer
        let written1 = data1.lock().unwrap().clone();
        let written2 = data2.lock().unwrap().clone();
        assert_eq!(written1, b"ACGT");
        assert!(written2.is_empty());
    }
}
