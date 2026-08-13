use clap::{Parser, Subcommand};

mod fastq_quality;
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
}

fn main() -> Result<(), std::io::Error> {
    let cli = Cli::parse();
    let _ = match cli.command {
        Commands::Mask(args) => mask::mask(args),
        _ => Result::Ok(()),
    };

    Ok(())
}
