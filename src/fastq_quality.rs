use clap::Args;
use core::{fmt, panic};
use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::PathBuf;

//TODO: should I add high precision floats

#[derive(Args, Debug)]
pub struct FastqQualityArgs {
    #[arg(short, long, num_args=1..)]
    fastq: Vec<PathBuf>,

    // Need to veriy names are same length as paths passed if passed
    #[arg(short, long, num_args=1..)]
    names: Option<Vec<String>>,

    #[arg(short, long)]
    output: Option<PathBuf>,
}

const PHRED_VALUES: usize = 93;
const PHRED_OFFSET: u8 = 33;

type PhredArray = [usize; PHRED_VALUES];

#[derive(Debug)]
pub struct FastQQuality {
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
    mean_sequence_length: f64,
    std_sequence_length: f64,
    min_sequence_length: usize,
    max_sequence_length: usize,
}

impl fmt::Display for FastQQuality {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "\"{name}\": {{
                \"total_bp\": {total_bp},
                \"total_reads\": {total_reads},
                \"qual_min\": {qual_min},
                \"qual_max\": {qual_max},
                \"qual_sum\": {qual_sum},
                \"qual_mean\": {qual_mean:.2},
                \"qual_std\": {qual_std:.2},
                \"read_qual_mean\": {read_qual_mean:.2},
                \"read_qual_std\": {read_qual_std:.2},
                \"mean_sequence_length\": {mean_sequence_length:.2},
                \"min_sequence_length\": {min_sequence_length},
                \"max_sequence_length\": {max_sequence_length},
                \"std_sequence_length\": {std_sequence_length:.2},
                }}",
            name = self.name,
            total_bp = self.total_bp,
            total_reads = self.total_reads,
            qual_min = self.qual_min,
            qual_max = self.qual_max,
            qual_sum = self.qual_sum,
            qual_mean = self.qual_mean,
            qual_std = self.qual_std,
            read_qual_mean = self.read_qual_mean,
            read_qual_std = self.read_qual_std,
            mean_sequence_length = self.mean_sequence_length,
            min_sequence_length = self.min_sequence_length,
            max_sequence_length = self.max_sequence_length,
            std_sequence_length = self.std_sequence_length,
        )
    }
}

impl Default for FastQQuality {
    fn default() -> Self {
        Self {
            quality: [0; PHRED_VALUES],
            name: String::new(),
            total_reads: 0,
            total_bp: 0,
            qual_sum: 0,
            qual_max: 0,
            qual_min: 0,
            qual_mean: 0.0,
            qual_std: 0.0,
            read_qual_mean: 0.0,
            read_qual_std: 0.0,
            sequence_lengths: HashMap::new(),
            avg_qual_read: HashMap::new(),
            mean_sequence_length: 0.,
            std_sequence_length: 0.,
            min_sequence_length: 0,
            max_sequence_length: 0,
        }
    }
}

impl FastQQuality {
    fn combine_maps<T: Eq + std::hash::Hash + Clone + Copy>(
        d1: &HashMap<T, usize>,
        d2: &HashMap<T, usize>,
    ) -> HashMap<T, usize> {
        let mut new_map = d1.clone();
        d2.iter().for_each(|(key, value)| {
            new_map
                .entry(*key)
                .and_modify(|nv| *nv += *value)
                .or_insert(*value);
        });
        new_map
    }

    fn calculate_mean_min_max_scalar(
        &self,
        values: &HashMap<usize, usize>,
        total: usize,
    ) -> (f64, usize, usize) {
        let mut sum: usize = 0;
        let mut min_length: usize = usize::MAX;
        let mut max_length: usize = 0;

        for (k, v) in values {
            sum += (*k) * (*v);
            if *k > max_length {
                max_length = *k;
            }
            if *k < min_length {
                min_length = *k;
            }
        }

        (sum as f64 / total as f64, min_length, max_length)
    }

    fn calculate_mean_min_max_float(
        &self,
        values: &HashMap<usize, usize>,
        total: usize,
    ) -> (f64, usize, usize) {
        let mut sum: f64 = 0.;
        let mut min_length: usize = usize::MAX;
        let mut max_length: usize = 0;

        for (k, v) in values {
            sum += f64::from_bits(*k as u64) * (*v) as f64;
            if *k > max_length {
                max_length = *k;
            }
            if *k < min_length {
                min_length = *k;
            }
        }

        (sum / total as f64, min_length, max_length)
    }

    fn set_per_read_quality_metrics(&mut self) {
        let (mean, _, _) = self.calculate_mean_min_max_float(&self.avg_qual_read, self.total_reads);
        let std = self.calculate_std_reads(&self.avg_qual_read, mean, self.total_reads);
        self.read_qual_mean = mean;
        self.read_qual_std = std;
    }

    fn set_quality_min_max(&mut self) {
        let mut max_idx: usize = 0;
        let mut min_idx: usize = usize::MAX;

        for (k, v) in self.quality.iter().enumerate() {
            if min_idx == usize::MAX && *v != 0 {
                min_idx = k;
            }
            if *v != 0 {
                max_idx = k;
            }
        }
        self.qual_min = min_idx;
        self.qual_max = max_idx;
    }

    fn set_qual_mean(&mut self) {
        self.qual_mean = self.qual_sum as f64 / self.total_bp as f64;
    }

    fn set_qual_std(&mut self) {
        let mut sum_of_squares = 0.;
        for (k, v) in self.quality.iter().enumerate() {
            if *v != 0 {
                let diff = k as f64 - self.qual_mean;
                let diff_sq = diff.powi(2);
                sum_of_squares += diff_sq;
            }
        }

        let normalized_squares = sum_of_squares / self.total_bp as f64;
        self.qual_std = normalized_squares.sqrt();
    }

    fn set_read_length_metrics(&mut self) {
        let (mean, min, max) =
            self.calculate_mean_min_max_scalar(&self.sequence_lengths, self.total_reads);
        let std = self.calculate_std_reads(&self.sequence_lengths, mean, self.total_reads);
        self.mean_sequence_length = mean;
        self.min_sequence_length = min;
        self.max_sequence_length = max;
        self.std_sequence_length = std;
    }

    fn set_quality_metrics(&mut self) {
        self.set_quality_min_max();
        self.set_qual_mean();
        self.set_qual_std();
    }

    fn calculate_std_reads(&self, values: &HashMap<usize, usize>, mean: f64, total: usize) -> f64 {
        let mut sum_of_squares: f64 = 0.;
        for (k, v) in values {
            let diff_sq = (f64::from_bits(*k as u64) - mean).powi(2);
            sum_of_squares += diff_sq * (*v) as f64;
        }

        let normalized_squares = sum_of_squares / total as f64;
        normalized_squares.sqrt()
    }

    fn combine_quality(&self, other_qual: PhredArray) -> PhredArray {
        let mut new_phred = [0; PHRED_VALUES];
        for i in 0..PHRED_VALUES {
            new_phred[i] = self.quality[i] + other_qual[i];
        }
        new_phred
    }

    fn update_quality_record(&mut self, quality: &[u8]) {
        let mut qual_sum: usize = 0;
        for phred_score in quality.iter() {
            let idx: usize = (phred_score - PHRED_OFFSET) as usize;
            qual_sum += idx;
            self.quality[idx] += 1;
        }
        let read_average_quality = qual_sum as f64 / quality.len() as f64;
        self.avg_qual_read
            .entry(read_average_quality.to_bits() as usize)
            .and_modify(|nv| *nv += 1)
            .or_insert(1);
        self.qual_sum += qual_sum;
    }

    fn add_record(&mut self, seq_analysis: &bio::io::fastq::Record) {
        self.total_reads += 1;
        self.total_bp += seq_analysis.seq().len();
        self.update_quality_record(seq_analysis.qual());
        self.sequence_lengths
            .entry(seq_analysis.seq().len())
            .and_modify(|nv| *nv += 1)
            .or_insert(1);
    }

    fn add(&mut self, other: FastQQuality) {
        // Combine data from 2 reads

        let n_new_quality = self.combine_quality(other.quality);

        let n_sequence_lengths =
            FastQQuality::combine_maps(&self.sequence_lengths, &other.sequence_lengths);
        let n_avg_qual_read = FastQQuality::combine_maps(&self.avg_qual_read, &other.avg_qual_read);

        self.quality = n_new_quality;
        self.sequence_lengths = n_sequence_lengths;
        self.avg_qual_read = n_avg_qual_read;
        self.total_reads += other.total_reads;
        self.total_bp += other.total_bp;
        self.qual_sum += other.qual_sum;
    }

    fn set_statistics(&mut self) {
        self.set_read_length_metrics();
        self.set_quality_metrics();
        self.set_per_read_quality_metrics();
    }

    pub fn get_fastq_stats(input_args: FastqQualityArgs) -> Result<(), std::io::Error> {
        // Can get either the file name or names name if passed
        // Need to unpack names if exists

        let output_file: Box<dyn Write> = match input_args.output {
            Some(arg) => Box::new(File::create(arg).unwrap()),
            None => Box::new(std::io::stdout().lock()),
        };
        let mut output_writer = BufWriter::with_capacity(1024 * 1024, output_file);

        let names = match input_args.names.as_ref() {
            Some(nv) => nv,
            None => &std::vec::Vec::new(),
        };
        let mut data_out: Vec<String> = std::vec::Vec::new();
        let mut combined_data = FastQQuality {
            name: "combined".to_string(),
            ..Default::default()
        };

        for (k, v) in input_args.fastq.iter().enumerate() {
            let fallback = v.file_prefix().unwrap().to_string_lossy().into_owned();
            let sample_name: &String = match names.get(k) {
                Some(name) => name,
                None => &fallback,
            };
            let mut data = FastQQuality {
                name: sample_name.clone(),
                ..Default::default()
            };
            let reader =
                bio::io::fastq::Reader::from_bufread(BufReader::new(get_reader(v.to_path_buf())));
            for record in reader.records() {
                data.add_record(&record.unwrap());
            }
            data.set_statistics();
            data_out.push(format!("{}", data));
            combined_data.add(data);
        }
        for record in data_out {
            writeln!(output_writer, "{}", record)?;
        }
        combined_data.set_statistics();
        writeln!(output_writer, "{}", combined_data)?;
        output_writer
            .flush()
            .expect("Could not output sequence buffer.");

        Ok(())
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
