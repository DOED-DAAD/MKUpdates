use clap::Args;
use paste::paste;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

const REPORT_TAXON_DELIMITER: char = ';';
const TAXON_LEVEL_SPLIT: &str = "__";
const ALTERNATE_TAXA_ALLOWED: u32 = 1;

#[derive(Args, Debug)]
pub struct MashTopHitArgs {
    /// Input file path
    #[arg(short, long)]
    pub input: String,

    #[arg(short, long, default_value_t = 0.90)]
    pub percent_identity_cutoff: f64,
}

#[derive(Args, Debug)]
pub struct MashMetagenomicArgs {
    /// Input file path
    #[arg(short, long)]
    pub input: PathBuf,

    #[arg(short, long)]
    pub equivalent_taxa: PathBuf,

    #[arg(short, long, default_value_t = 'g')]
    taxonomic_level: char,

    #[arg(short, long, default_value_t = 0.90)]
    pub percent_identity_cutoff: f64,
}

#[derive(Debug)]
#[allow(dead_code)]
struct MashHit {
    identity: f64,
    shared_hashes: (u32, u32),
    median_multiplicity: u32,
    p_value: f64,
    query_id: String,
    query_note: String,
}

macro_rules! generate_mash_sort_fn {
    ($field_name:tt $(, $field_val:tt)? ) => {
        paste! {
            fn [<sort_mash_ $field_name>](mash_records: &mut [MashHit]) {
                mash_records.sort_by(|a, b| {
                a.$field_name $( .$field_val )?
                    .partial_cmp(&b.$field_name $( .$field_val )? )
                    .unwrap_or(std::cmp::Ordering::Equal)
                });
            }
        }
    };
}

generate_mash_sort_fn!(identity);
generate_mash_sort_fn!(p_value);
generate_mash_sort_fn!(shared_hashes, 0);

fn read_equivalent_json(
    equivalent_data: PathBuf,
) -> Result<HashMap<String, Vec<String>>, Box<dyn std::error::Error>> {
    let file = File::open(equivalent_data)?;
    let reader = BufReader::new(file);
    let json_data: HashMap<String, Vec<String>> = serde_json::from_reader(reader)?;
    Ok(json_data)
}

fn get_equivalent_taxa(
    equivalent_taxa: PathBuf,
    taxonomic_level: char,
) -> HashMap<String, HashSet<String>> {
    /*
            The below operation takes in the array of equivalent taxonomy values from  `equivalent_taxa.json`
            an example of the array can be seen below:
            [
                f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli",
                f__Enterobacteriaceae;g__Shigella;s__Shigella dysenteriae"
            ]

            The map process, splits each value in the array on the taxon delimiter value e.g. ';'
                e.g.
                    [
                        [f__Enterobacteriaceae, g__Escherichia, s__Escherichia coli],
                        [f__Enterobacteriaceae, g__Shigella, s__Shigella dysenteriae]
                    ]

            The values split on the delimited will then be in separate lists, that need to be flattened
            which itertools.chain.from_iterable does for us.

            The listed values are then iterated over and only the taxonomic level
            of interest (specified taxonomic_classification_level) is stored in the dictionary of equivalent values.

            Only the taxonomic level of interest is stored to prevent the stored sets from increasing in size.

            The final result is a dictionary with the hits of:
            taxa_equivalent: {
                Escherichia: { g__Escherichia, g__Shigella},
                Example: { g__genera1, g__genera2},
            }
    */

    let json_data = read_equivalent_json(equivalent_taxa).unwrap();
    let mut equivalent_taxa: HashMap<String, HashSet<String>> = HashMap::new();

    for (key, value) in &json_data {
        for taxa_string in value.iter() {
            // assuming below actually executes, as map is lazy.
            for taxa in taxa_string
                .split(REPORT_TAXON_DELIMITER)
                .filter(|x| x.starts_with(taxonomic_level))
            {
                if let Some(set) = equivalent_taxa.get_mut(key) {
                    set.insert(taxa.to_string());
                } else {
                    let mut set = HashSet::new();
                    set.insert(taxa.to_string());
                    equivalent_taxa.insert(key.clone(), set);
                }
            }
        }
    }

    equivalent_taxa
}

fn get_taxa_level(taxa_string: &str) -> String {
    match taxa_string.find(TAXON_LEVEL_SPLIT) {
        Some(index) => taxa_string[0..index].to_string(),
        None => panic!("Could not identify taxa level in string."),
    }
}

fn parse_flatten_queries(mash_data: Vec<MashHit>) -> HashMap<String, HashSet<String>> {
    let mut flattened_queries: HashMap<String, HashSet<String>> = HashMap::new();

    for mash_hit in mash_data {
        mash_hit
            .query_note
            .split(REPORT_TAXON_DELIMITER)
            .for_each(|x| {
                flattened_queries
                    .entry(get_taxa_level(x))
                    .or_default()
                    .insert(x.to_string());
            });
    }
    flattened_queries
    // Now flattened_queries contains a mapping from query_id to a vector of MashHit records.
}

fn normalize_taxa(
    taxa_levels: &mut HashMap<String, HashSet<String>>,
    taxa_level: String,
    equivalent_taxa: &HashMap<String, HashSet<String>>,
) {
    let levels = taxa_levels
        .get_mut(&taxa_level)
        .expect("Could not find taxa level in equivalent taxa.");
    for (key, value) in equivalent_taxa {
        if !value.is_disjoint(levels) {
            dbg!("No value: {}", value);
            levels.retain(|item| !value.contains(item));
            levels.insert(key.clone());
        }
    }
}

impl MashHit {
    fn new(
        identity: f64,
        shared_hashes: (u32, u32),
        median_multiplicity: u32,
        p_value: f64,
        query_id: String,
        query_note: String,
    ) -> Self {
        MashHit {
            identity,
            shared_hashes,
            median_multiplicity,
            p_value,
            query_id,
            query_note,
        }
    }
}

fn get_shared_hashes(s: &str) -> (u32, u32) {
    let parts: Vec<&str> = s.split('/').collect();
    if parts.len() != 2 {
        panic!("Invalid shared hashes format: {}", s);
    }
    let a = parts[0]
        .parse::<u32>()
        .expect("Could not parse numerator of shared hashes.");
    let b = parts[1]
        .parse::<u32>()
        .expect("Could not parse denominator of shared hashes");
    (a, b)
}

fn parse_mash(file: PathBuf, id_cutoff: f64) -> Vec<MashHit> {
    let file = File::open(file).unwrap();
    let mut reader = BufReader::new(file);
    let mut buf = String::new();

    let mut mash_records: Vec<MashHit> = Vec::new();

    while reader.read_line(&mut buf).unwrap() > 0 {
        let fields: Vec<&str> = buf.trim().split('\t').collect();
        if fields.len() < 5 {
            panic!("Invalid line: {}", buf);
        }

        let identity = fields[0].parse::<f64>().expect("Failed to parse identity");
        if identity < id_cutoff {
            buf.clear();
            continue;
        }
        let shared_hashes = get_shared_hashes(fields[1]);
        let median_multiplicity = fields[2]
            .parse::<u32>()
            .expect("Failed to parse median multiplicity");
        let p_value = fields[3].parse::<f64>().expect("Failed to parse p-value");
        let query_id = fields[4];
        let query_note = fields[5];
        let mash_record = MashHit::new(
            identity,
            shared_hashes,
            median_multiplicity,
            p_value,
            query_id.to_string(),
            query_note.to_string(),
        );
        mash_records.push(mash_record);
        buf.clear();
    }

    mash_records
}

pub fn get_mash_top_hit(args: MashTopHitArgs) -> Result<(), std::io::Error> {
    let file_path = PathBuf::from(args.input);
    let mut mash_data = parse_mash(file_path, args.percent_identity_cutoff);
    if mash_data.is_empty() {
        println!("No Species Identified");
    }
    sort_mash_identity(&mut mash_data);
    sort_mash_p_value(&mut mash_data);
    sort_mash_shared_hashes(&mut mash_data);
    // Verified that there is always a value in the vector in the above code
    let mut best_option = mash_data
        .first()
        .expect("Could not identify top hit option.")
        .query_note
        .split(REPORT_TAXON_DELIMITER)
        .next_back()
        .expect("Could not find taxon delimiter in output string")
        .to_string();
    best_option = best_option.replace("\"", "");
    println!("{}", best_option);
    Ok(())
}

pub fn check_metagenomic(args: MashMetagenomicArgs) -> Result<(), std::io::Error> {
    let mash_data = parse_mash(args.input, args.percent_identity_cutoff);
    let equivalent_taxa: HashMap<String, HashSet<String>> =
        get_equivalent_taxa(args.equivalent_taxa, args.taxonomic_level);
    let mut flattened_queries = parse_flatten_queries(mash_data);
    normalize_taxa(
        &mut flattened_queries,
        args.taxonomic_level.to_string(),
        &equivalent_taxa,
    );
    let taxa_at_check_level = flattened_queries
        .get(&args.taxonomic_level.to_string())
        .expect("Could not find taxa level in flattened queries.")
        .len();

    if taxa_at_check_level > ALTERNATE_TAXA_ALLOWED as usize {
        println!("true");
    } else {
        println!("false");
    }

    Ok(())
}
