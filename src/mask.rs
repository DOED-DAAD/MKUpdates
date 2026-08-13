use clap::Args;
use rust_htslib::faidx::{Reader, build};
use rust_htslib::{bam, bam::Read};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

fn in_range(s: &str) -> Result<f32, String> {
    let min = 0.0;
    let max = 1.0;
    let val: f32 = s
        .parse::<f32>()
        .map_err(|_| format!("Value must be between {} and {}", min, max))
        .unwrap();

    if val >= min && val <= max {
        Ok(val)
    } else {
        Err(format!("Value must be between {} and {}", min, max))
    }
}

#[derive(Args, Debug)]
pub struct MaskArgs {
    #[arg(short, long)]
    bam: PathBuf, // Bam path

    #[arg(short, long)]
    fasta: PathBuf, // Fasta path, Gzip is not supported only bgzip.

    #[arg(short = 'w', long, default_value_t = 1)]
    fwd_depth: u32, // Minimum read depth from forward reads

    #[arg(short, long, default_value_t = 1)]
    rev_depth: u32, // Minimum read-depth support from reverse reads

    #[arg(short, long, default_value_t = 0.75, value_parser = in_range)]
    pileup_support: f32, // Minimum depth of position needed to be supported by pileup to not be
    // called an N. Value is passes is a proportion between 0 and 1.
    #[arg(short, long, default_value_t = 2)]
    min_depth: u32,
}

fn read_indexed_fasta_file(input_fasta: PathBuf) -> Reader {
    Reader::from_path(input_fasta).expect("Failed to open faidx.")
}

fn read_bam_file(input_bam: PathBuf) -> bam::IndexedReader {
    bam::IndexedReader::from_path(input_bam).expect("Could not open bam file.")
}

pub fn mask(args: MaskArgs) -> Result<(), std::io::Error> {
    let fasta_input = args.fasta.clone();
    let fasta_input_size =
        std::fs::metadata(fasta_input.clone()).expect("Could not get fasta file size.");
    let fai_idx_path = fasta_input.with_added_extension("fai");

    if !fai_idx_path.exists() {
        build(&args.fasta).expect("Failed to build fasta index.")
    }

    let fasta_reader = read_indexed_fasta_file(args.fasta);

    let bam_input = args.bam.clone();
    let bam_idx_path = bam_input.with_added_extension(".bai");
    if !bam_idx_path.exists() {
        bam::index::build(bam_input, None, bam::index::Type::Bai, 2)
            .expect("Could not find or build bam index.");
    }

    let mut bam_reader = read_bam_file(args.bam);
    bam_reader
        .set_threads(2)
        .expect("Could not allocate 2 threads for bam reader.");

    let header_copy = rust_htslib::bam::header::Header::from_template(bam_reader.header());
    let header_hashmap = header_copy.to_hashmap();
    let mut output_fasta = fasta_input.clone();
    output_fasta.set_extension("masked.fasta");

    // TODO: add check to verify sequence lines present
    // TODO: Verify only primary alignments used in the future
    let fasta_out = File::create(output_fasta).expect("Could not create output file.");
    let mut writer = BufWriter::with_capacity(1024 * 1024, fasta_out);
    let mut log_writer =
        BufWriter::with_capacity(fasta_input_size.len() as usize, std::io::stdout().lock());

    if let Some(sq_records) = header_hashmap.get("SQ") {
        for record in sq_records {
            let name = &record["SN"];
            bam_reader
                .fetch(name)
                .expect("Sequence in header is not present in bam file.");

            let mut seq = fasta_reader
                .fetch_seq(name, 0, fasta_reader.fetch_seq_len(name) as usize)
                .expect("Sequence in bam is not in fasta.");

            for p in bam_reader.pileup() {
                let pileup = p.expect("Missing sequence in bam from header.");

                // Closure for writing log outputs
                let mut mask_pos = |pos: usize, seq: &mut Vec<u8>| {
                    unsafe { *seq.get_unchecked_mut(pos) = b'N' };
                    writeln!(
                        log_writer,
                        "Sequence {} position {} is not supported changing to N.",
                        name, pos
                    )
                    .unwrap();
                };

                // Need to check that the variant is of sufficient depth and forward and reverse
                // depths are met
                if pileup.depth() <= args.min_depth {
                    mask_pos(pileup.pos() as usize, &mut seq);
                    continue;
                }

                let comp_base = unsafe { *seq.get_unchecked(pileup.pos() as usize) };
                let mut fwd_count: u32 = 0;
                let mut rev_count: u32 = 0;
                let mut comp_base_count: u32 = 0;
                let minimum_level = (args.pileup_support * pileup.depth() as f32) as u32;
                for alignment in pileup.alignments() {
                    if let Some(qpos) = alignment.qpos()
                        && alignment.record().seq()[qpos] == comp_base
                    {
                        comp_base_count += 1;
                    }

                    rev_count += alignment.record().is_reverse() as u32; // Adds 1 if true 0 if false
                    fwd_count += !alignment.record().is_reverse() as u32; // adds 1 if false 0 if true
                }

                let passed_strand_direction =
                    fwd_count >= args.fwd_depth && rev_count >= args.rev_depth;

                if comp_base_count <= minimum_level || !passed_strand_direction {
                    mask_pos(pileup.pos() as usize, &mut seq);
                }
            }
            writer.write_all(b">")?;
            writer.write_all(name.as_bytes())?;
            writer.write_all(b"\n")?;
            writer.write_all(&seq)?;
            writer.write_all(b"\n")?;
        }
    }
    writer
        .flush()
        .expect("Could not flush output sequence buffer.");

    Ok(())
}
