use clap::{Parser, Subcommand};

mod fastq_quality;
mod mash;
mod mask;

#[derive(Parser, Debug)]
#[command(version, about, long_about=None, arg_required_else_help=true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    Mask(mask::MaskArgs),
    FastqQuality(fastq_quality::FastqQualityArgs),
    MashTopHit(mash::MashTopHitArgs),
    MashMetagenomicArgs(mash::MashMetagenomicArgs),
}

fn main() -> Result<(), std::io::Error> {
    let cli = Cli::parse();
    let _ = match cli.command {
        Commands::Mask(args) => mask::mask(args),
        Commands::FastqQuality(args) => fastq_quality::FastQQuality::get_fastq_stats(args),
        Commands::MashTopHit(args) => mash::get_mash_top_hit(args),
        Commands::MashMetagenomicArgs(args) => mash::check_metagenomic(args),
    };

    Ok(())
}
