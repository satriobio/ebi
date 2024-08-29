extern crate seal;
extern crate bio;
extern crate rayon;

use seal::pair::{Alignment, AlignmentSet, InMemoryAlignmentMatrix, NeedlemanWunsch, SmithWaterman, Step};
use bio::io::fasta;
use rayon::prelude::*;
use std::collections::HashMap;
use rayon::ThreadPoolBuilder;


// Function to read FASTA file and return sequences as a HashMap
fn read_fasta(filename: &str) -> HashMap<String, String> {
    let reader = fasta::Reader::from_file(filename).expect("Failed to read FASTA file");
    reader.records()
        .map(|record| {
            let record = record.expect("Error reading record");
            (
                record.id().to_owned(),
                String::from_utf8_lossy(record.seq()).to_string() // Convert &[u8] to String
            )
        })
        .collect()
}

// Function to align a read against all reference sequences and store the best alignment
fn align_reads(reads: &HashMap<String, String>, references: &HashMap<String, String>, use_local: bool) -> HashMap<String, (String, isize)> {
    let strategy = SmithWaterman::new(2, -1, -1, -1);

    // Create a custom thread pool with 8 threads
    let pool = ThreadPoolBuilder::new().num_threads(16).build().unwrap();
    pool.install(|| {
        // Use Rayon to parallelize the computation
        let best_alignments: HashMap<String, (String, isize)> = reads.par_iter().map(|(read_name, read_sequence)| {
            let read_chars: Vec<char> = read_sequence.chars().collect();
            let mut best_score = isize::MIN;
            let mut best_ref_name = String::new();

            // Parallelize the processing of references for each read
            let (best_ref_name_for_read, best_score_for_read) = references.par_iter().map(|(ref_name, ref_sequence)| {
                let ref_chars: Vec<char> = ref_sequence.chars().collect();
                let set: AlignmentSet<InMemoryAlignmentMatrix> =
                    AlignmentSet::new(read_chars.len(), ref_chars.len(), strategy.clone(), |x, y| {
                        read_chars[x] == ref_chars[y]
                    })
                    .unwrap();

                let alignment = if use_local {
                    set.local_alignment()
                } else {
                    set.global_alignment()
                };

                let score = alignment.score();
                (ref_name.clone(), score)
            }).max_by_key(|(_, score)| *score) // Get the reference with the best score
            .unwrap_or((String::new(), isize::MIN)); // Fallback if no references are found

            (read_name.clone(), (best_ref_name_for_read, best_score_for_read))
        }).collect();

        best_alignments
    })
}

// Print the best alignment results
fn print_best_alignments(best_alignments: &HashMap<String, (String, isize)>) {
    for (read_name, (ref_name, score)) in best_alignments {
        println!("Read: {}, Best Reference: {}, Score: {}", read_name, ref_name, score);
    }
}

fn main() {
    let reads_filename = "/query.fa";
    let references_filename = "gencode.v46.transcripts.200bp.fa";


    let reads = read_fasta(reads_filename);
    let references = read_fasta(references_filename);

    let best_alignments = align_reads(&reads, &references, true);  // Set to `false` for global alignment

    print_best_alignments(&best_alignments);
}