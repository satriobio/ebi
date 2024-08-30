use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use dlopen::raw::Library;
use std::collections::HashMap;
use std::ffi::CString;
use std::os::raw::{c_int, c_uint, c_uchar};
use bio::io::fasta;
use std::env;
use std::io::Write;

#[repr(C)]
struct CAlignRes {
    nScore: c_uint,
    nScore2: c_uint,
    nRefBeg: c_int,
    nRefEnd: c_int,
    nQryBeg: c_int,
    nQryEnd: c_int,
    nRefEnd2: c_int,
    sCigar: *const c_uint,
    nCigarLen: c_int,
}

#[repr(C)]
struct CProfile {
    pByte: *const c_int,
    pWord: *const c_int,
    pRead: *const c_uchar,
    pMat: *const c_uchar,
    nReadLen: c_int,
    nN: c_int,
    nBias: c_uchar,
}

struct CSsw {
    ssw: Library,
    ssw_init: unsafe extern fn(*const c_uchar, c_int, *const c_uchar, c_int, c_uchar) -> *mut CProfile,
    init_destroy: unsafe extern fn(*mut CProfile),
    ssw_align: unsafe extern fn(*mut CProfile, *const c_uchar, c_int, c_uchar, c_uchar, c_uchar, c_uint, c_int, c_int) -> *mut CAlignRes,
    align_destroy: unsafe extern fn(*mut CAlignRes),
}

impl CSsw {
    fn new(s_lib_path: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let lib = Library::open(s_lib_path)?;

        unsafe {
            let ssw_init = lib.symbol("ssw_init")?;
            let init_destroy = lib.symbol("init_destroy")?;
            let ssw_align = lib.symbol("ssw_align")?;
            let align_destroy = lib.symbol("align_destroy")?;

            Ok(CSsw {
                ssw: lib,
                ssw_init,
                init_destroy,
                ssw_align,
                align_destroy,
            })
        }
    }
}

fn read_fasta(file_path: &str) -> impl Iterator<Item = (String, String)> {
    let reader = fasta::Reader::from_file(file_path).expect("Failed to read FASTA file");

    reader.records().map(|record| {
        let record = record.expect("Failed to read record");
        (record.id().to_string(), String::from_utf8(record.seq().to_vec()).expect("Invalid UTF-8"))
    })
}

fn to_int(seq: &str, d_ele2int: &HashMap<String, usize>, default: usize) -> Vec<i8> {
    seq.chars().map(|c| {
        d_ele2int.get(&c.to_string()).copied().unwrap_or(default) as i8
    }).collect()
}

fn align_one(
    ssw: &CSsw,
    q_profile: *mut CProfile,
    r_num: &[i8],
    n_r_len: i32,
    n_open: u8,
    n_ext: u8,
    n_flag: u8,
    n_mask_len: i32
) -> (u32, u32, i32, i32, i32, i32, i32, i32, Vec<u32>) {
    // Ensure q_profile is not null
    if q_profile.is_null() {
        eprintln!("Error: q_profile is null.");
        return (0, 0, 0, 0, 0, 0, 0, 0, Vec::new());
    }

    // Call the ssw_align function
    let res_ptr = unsafe { (ssw.ssw_align)(
        q_profile,
        r_num.as_ptr() as *const c_uchar,
        n_r_len,
        n_open,
        n_ext,
        n_flag,
        0,
        0,
        n_mask_len
    )};

    // Ensure alignment result is not null
    if res_ptr.is_null() {
        eprintln!("Error: ssw_align returned null.");
        // Call destroy functions to avoid memory leaks
        unsafe { (ssw.init_destroy)(q_profile) };
        return (0, 0, 0, 0, 0, 0, 0, 0, Vec::new());
    }

    // Safely dereference the result pointer
    let res = unsafe { &*res_ptr };
    let s_cigar_len = res.nCigarLen as usize;

    // Collect cigar values if sCigar is not null
    let l_cigar: Vec<u32> = if s_cigar_len > 0 && !res.sCigar.is_null() {
        (0..s_cigar_len).map(|i| unsafe { *res.sCigar.add(i) }).collect()
    } else {
        Vec::new()
    };

    // Prepare result tuple
    let result = (
        res.nScore,
        res.nScore2,
        res.nRefBeg,
        res.nRefEnd,
        res.nQryBeg,
        res.nQryEnd,
        res.nRefEnd2,
        s_cigar_len as i32,
        l_cigar
    );

    // Clean up resources using the mutable pointer
    unsafe { (ssw.align_destroy)(res_ptr) };

    result
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 {
        eprintln!("Usage: {} <reference_fasta> <query_fasta> --thread <num_threads>", args[0]);
        std::process::exit(1);
    }

    let reference_path = &args[1];
    let query_fasta = &args[2];
    let num_threads: usize = args.iter()
        .position(|arg| arg == "--thread")
        .and_then(|i| args.get(i + 1))
        .and_then(|s| s.parse().ok())
        .unwrap_or(8); // Default to 8 threads if not specified

    ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .expect("Failed to create Rayon thread pool");

    let ssw = unsafe { CSsw::new("./libssw.so").expect("Failed to load SSW library") };

    let l_ele = ['A', 'C', 'G', 'T', 'N'];
    let mut d_ele2int = HashMap::new();
    let mut d_int2ele = HashMap::new();

    for (i, &ele) in l_ele.iter().enumerate() {
        d_ele2int.insert(ele.to_string(), i);
        d_ele2int.insert(ele.to_ascii_lowercase().to_string(), i);
        d_int2ele.insert(i, ele);
    }

    let n_ele_num = l_ele.len();
    let l_score: Vec<i8> = (0..n_ele_num).flat_map(|i| {
        (0..n_ele_num).map(move |j| {
            if l_ele[i] == l_ele[j] { 2 } else { -2 }
        })
    }).map(|x| x as i8).collect();

    let mat = l_score;
    let r_num_vec: Vec<(String, Vec<i8>)> = read_fasta(reference_path)
        .map(|(id, seq)| (id, to_int(&seq, &d_ele2int, l_ele.len())))
        .collect();

    let best_alignments: HashMap<String, (String, u32)> = read_fasta(query_fasta)
        .par_bridge()
        .map(|(s_q_id, s_q_seq)| {
            let q_num = to_int(&s_q_seq, &d_ele2int, l_ele.len());
            let q_profile = unsafe { (ssw.ssw_init)(
                q_num.as_ptr() as *const c_uchar,
                q_num.len() as c_int,
                mat.as_ptr() as *const c_uchar,
                l_ele.len() as c_int,
                1
            )};

            if q_profile.is_null() {
                eprintln!("Error: q_profile is null for query {}", s_q_id);
                return (s_q_id, (String::new(), u32::MIN));
            }

            let n_mask_len = (s_q_seq.len() / 2) as i32;
            let mut best_score = u32::MIN;
            let mut best_ref_id = String::new();

            for (s_r_id, r_num) in &r_num_vec {
                if r_num.is_empty() {
                    eprintln!("Error: r_num is empty for target {}", s_r_id);
                    continue;
                }

                unsafe {
                    let (score, _, _, _, _, _, _, _, _) = align_one(
                        &ssw, q_profile, r_num, r_num.len() as i32, 3, 1, 0, n_mask_len
                    );
                    if score == u32::MAX {
                        eprintln!("Error: Invalid score value.");
                        continue;
                    }
                    if score > best_score {
                        best_score = score;
                        best_ref_id = s_r_id.clone();
                    }
                }
            }

            unsafe { (ssw.init_destroy)(q_profile); }

            (s_q_id, (best_ref_id, best_score))
        })
        .collect();

    let stdout = std::io::stdout();
    let mut handle = stdout.lock();
    for (query_id, (ref_id, score)) in best_alignments {
        writeln!(handle, "{}\t{}\t{}", query_id, ref_id, score).expect("Failed to write to stdout");
    }
}
