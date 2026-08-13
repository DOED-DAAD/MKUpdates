use bio::io::fastq;
use clap::Args;
use core::panic;
use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::default;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::PathBuf;

//TODO: should I add high precision floats

#[derive(Args, Debug)]
pub struct FastqQualityArgs {
    #[arg(short, long, num_args=1..)]
    fastq: Vec<PathBuf>,

    #[arg(short, long, num_args=1..)]
    names: Option<Vec<String>>,

    #[arg(short, long)]
    output: PathBuf,
}

const PHRED_VALUES: usize = 93;

type PhredArray = [usize; PHRED_VALUES];

struct FastQQuality {
    name: String,
    total_reads: usize,
    total_bp: usize,
    qual_sum: usize,
    qual_max: usize,
    qual_min: usize,
    qual_mean: f64,
    qual_std: f64,
    read_qual_mean: f64,
    read_qual_std: f64,
    quality: PhredArray,
    sequence_lengths: HashMap<usize, usize>,
    avg_qual_read: HashMap<usize, usize>,
    mean_sequence_length: Option<f64>,
    std_sequence_length: Option<f64>,
    min_sequence_length: Option<usize>,
    max_sequence_length: Option<usize>,
}

fn combine_maps(d1: &HashMap<usize, usize>, d2: &HashMap<usize, usize>) -> HashMap<usize, usize> {
    let mut new_map = d1.clone();
    d2.iter().for_each(|(key, value)| {
        new_map
            .entry(*key)
            .and_modify(|nv| *nv += *value)
            .or_insert(*value);
    });
    new_map
}

impl Default for FastQQuality {
    fn default() -> Self {
        Self {
            quality: [0; PHRED_VALUES],
            ..Default::default()
        }
    }
}

impl FastQQuality {
    fn combine_quality(&self, other_qual: PhredArray) -> PhredArray {
        let mut new_phred = [0; PHRED_VALUES];
        for i in 0..PHRED_VALUES {
            new_phred[i] = self.quality[i] + other_qual[i];
        }
        new_phred
    }

    fn add(&self, other: FastQQuality) -> FastQQuality {
        // Combine data from 2 reads

        let n_new_quality = self.combine_quality(other.quality);

        let n_sequence_lengths = combine_maps(&self.sequence_lengths, &other.sequence_lengths);
        let n_avg_qual_read = combine_maps(&self.avg_qual_read, &other.avg_qual_read);

        FastQQuality {
            name: String::from("combined"),
            quality: n_new_quality,
            sequence_lengths: n_sequence_lengths,
            avg_qual_read: n_avg_qual_read,
            total_reads: self.total_reads + other.total_reads,
            total_bp: self.total_bp + other.total_bp,
            qual_sum: self.qual_sum + other.qual_sum,
            ..Default::default()
        }
    }
}

fn get_reader(file_path: PathBuf) -> Box<dyn Read> {
    let f = match File::open(&file_path) {
        Ok(file) => file,
        Err(error) => panic!("Could not open file: {:?}", error),
    };

    let extension = match file_path.extension() {
        Some(ext) => ext,
        None => panic!("Could not determine file extension: {:?}", file_path),
    }
    .to_str()
    .unwrap();

    match extension {
        "gz" => Box::new(GzDecoder::new(f)),
        _ => Box::new(f),
    }
}
