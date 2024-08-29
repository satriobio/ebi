use bio::io::fasta;
use bio::alignment::pairwise::*;
use bio::scores::blosum62;
use rayon::prelude::*;
use std::sync::{Arc, Mutex};

fn main() {
    // Load the reference sequences from a FASTA file
    let fasta_file = "/home/OXFORDNANOLABS/swibowo/git/ebi/gencode.v46.transcripts.200bp.fa";
    let reader = fasta::Reader::from_file(fasta_file).expect("Failed to read FASTA file");

    let references: Vec<_> = reader.records()
        .map(|r| {
            let record = r.expect("Error reading record");
            (record.id().to_owned(), record.seq().to_owned())
        })
        .collect();

    // Load the query sequence from a FASTA file
    let query_fasta_file = "/home/OXFORDNANOLABS/swibowo/git/ebi/query.fa";
    let query_reader = fasta::Reader::from_file(query_fasta_file).expect("Failed to read query FASTA file");

    let query_sequence: Vec<u8> = query_reader.records()
        .next()
        .expect("No query sequence found")
        .expect("Error reading query sequence")
        .seq()
        .to_owned();

    // Shared reference to store the highest score and alignment
    let highest_score_alignment = Arc::new(Mutex::new((0, String::new(), String::new())));

    // Set the number of threads for the Rayon thread pool
    let num_threads = 8; // Adjust the number of threads here
    // let pool = rayon::ThreadPoolBuilder::new()
    //     .num_threads(num_threads)
    //     .build()
    //     .unwrap();

    // pool.install(|| {
    //     // Align against each reference sequence in parallel
    //     references.par_iter().for_each(|(id, seq)| {
    //         let mut aligner = Aligner::with_capacity(query_sequence.len(), seq.len(), -5, -1, &blosum62);
    //         let alignment = aligner.local(&query_sequence, seq);

    //         let score = alignment.score;

    //         // Update the highest score and corresponding alignment
    //         let mut highest = highest_score_alignment.lock().unwrap();
    //         if score > highest.0 {
    //             *highest = (score, id.clone(), String::from_utf8_lossy(seq).into_owned());
    //         }
    //     });
    // });

    for (id, seq) in references.iter() {
        let mut aligner = Aligner::with_capacity(query_sequence.len(), seq.len(), -5, -1, &blosum62);
        let alignment = aligner.local(&query_sequence, seq);
    
        let score = alignment.score;
    
        // Update the highest score and corresponding alignment
        let mut highest = highest_score_alignment.lock().unwrap();
        if score > highest.0 {
            *highest = (score, id.clone(), String::from_utf8_lossy(seq).into_owned());
        }
    }    

    // Get the highest score alignment
    let highest = highest_score_alignment.lock().unwrap();
    let (score, ref_id, ref_seq) = &*highest;

    // Output the highest score alignment
    println!("Highest score: {}", score);
    println!("Aligned reference ID: {}", ref_id);
    println!("Reference sequence: {}", ref_seq);

    // Note: The BAM output step is not implemented, you may use appropriate libraries to write BAM files.
}
