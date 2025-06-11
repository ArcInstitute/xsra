use std::io::Write;

use anyhow::Result;
use ncbi_vdb_sys::SegmentType;
use serde::{Serialize, Serializer};

#[derive(Debug, PartialEq)]
pub struct SegmentTypeWrapper(SegmentType);
impl std::fmt::Display for SegmentTypeWrapper {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.0 {
            SegmentType::Technical => write!(f, "Technical"),
            SegmentType::Biological => write!(f, "Biological"),
        }
    }
}
impl Serialize for SegmentTypeWrapper {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        // Use the Display implementation we already have to convert to a string
        serializer.serialize_str(&self.to_string())
    }
}
impl From<SegmentType> for SegmentTypeWrapper {
    fn from(ty: SegmentType) -> Self {
        Self(ty)
    }
}

#[derive(Debug, Serialize, PartialEq)]
pub struct SegmentStats {
    sid: usize,
    segment_type: SegmentTypeWrapper,
    mean_length: f64,
    mean_quality: f64,
}
impl SegmentStats {
    pub fn new(sid: usize, segment_type: SegmentType, mean_length: f64, mean_quality: f64) -> Self {
        Self {
            sid,
            segment_type: segment_type.into(),
            mean_length,
            mean_quality,
        }
    }
}

#[derive(Debug, Serialize)]
pub struct DescribeStats {
    total_spots: usize,
    processed_spots: usize,
    spot_range: (usize, usize),
    stats: Vec<SegmentStats>,
}
impl DescribeStats {
    pub fn new(
        segment_types: Vec<SegmentType>,
        segment_lengths: Vec<Vec<f64>>,
        segment_qualities: Vec<Vec<f64>>,
        processed_spots: usize,
        first_spot: usize,
        last_spot: usize,
        total_spots: usize,
    ) -> Self {
        let num_segments = segment_types.len();
        let average_segment_lengths = Self::average(segment_lengths);
        let average_segment_qualities = Self::average(segment_qualities);
        let stats = (0..num_segments)
            .map(|idx| {
                SegmentStats::new(
                    idx,
                    segment_types[idx],
                    average_segment_lengths[idx],
                    average_segment_qualities[idx],
                )
            })
            .collect();
        Self {
            total_spots,
            processed_spots,
            spot_range: (first_spot, last_spot),
            stats,
        }
    }

    fn average(vec: Vec<Vec<f64>>) -> Vec<f64> {
        vec.iter()
            .map(|lens| {
                if lens.is_empty() {
                    0.0
                } else {
                    lens.iter().sum::<f64>() / lens.len() as f64
                }
            })
            .collect()
    }

    pub fn segment_lengths(&self) -> Vec<f64> {
        self.stats.iter().map(|s| s.mean_length).collect()
    }

    pub fn pprint<W: Write>(&self, wtr: &mut W) -> Result<()> {
        serde_json::to_writer_pretty(wtr, self)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // DescribeStats::new tests
    #[test]
    fn test_describe_stats_creation() {
        let segment_types = vec![SegmentType::Biological, SegmentType::Technical];
        let segment_lengths = vec![vec![100.0, 150.0], vec![20.0, 30.0]]; // means: 125.0, 25.0
        let segment_qualities = vec![vec![30.0, 34.0], vec![25.0, 25.0]]; // means: 32.0, 25.0
        let processed_spots = 10;
        let first_spot = 1;
        let last_spot = 11;
        let total_spots = 100;

        let stats = DescribeStats::new(
            segment_types,
            segment_lengths,
            segment_qualities,
            processed_spots,
            first_spot,
            last_spot,
            total_spots,
        );

        let expected_stats = vec![
            SegmentStats {
                sid: 0,
                segment_type: SegmentType::Biological.into(),
                mean_length: 125.0,
                mean_quality: 32.0,
            },
            SegmentStats {
                sid: 1,
                segment_type: SegmentType::Technical.into(),
                mean_length: 25.0,
                mean_quality: 25.0,
            },
        ];

        // Assert that the generated stats vector matches the expected vector
        assert_eq!(
            stats.stats, expected_stats,
            "The generated SegmentStats vector should match the expected calculations"
        );
    }
}
