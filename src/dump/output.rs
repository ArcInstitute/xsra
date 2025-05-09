use std::io::Write;
use std::sync::mpsc;
use std::sync::Arc;
use std::thread;

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
            let &(ref buffer, ref cvar) = &*buffer_pair_clone;

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
                handle.write_all(&data.drain(..).as_slice())?;
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
        let mut guard = buffer.lock();
        guard.extend_from_slice(data);
        cvar.notify_one();
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
        let thread_writers = segment_handles
            .into_iter()
            .map(|handle| ThreadWriter::new(handle))
            .collect();
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
