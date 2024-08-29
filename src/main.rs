use dlopen::raw::Library;
use std::collections::HashMap;
use std::ffi::CString;
use std::fs::File;
use std::io::{self, BufRead};
use std::os::raw::{c_int, c_uint, c_uchar};
use std::ptr;
use bio::io::fasta;

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
        // let s_lib_name = "./";
        // let path = if !s_lib_path.is_empty() {
        //     format!("{}/{}", s_lib_path, s_lib_name)
        // } else {
        //     s_lib_name.to_string()
        // };

        let path = "/mnt/869990e7-a61f-469f-99fe-a48d24ac44ca/git/ebi/libssw.so";

        let lib = Library::open(&path)?;

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


fn read_matrix(file_path: &str) -> (Vec<String>, HashMap<String, usize>, HashMap<usize, String>, Vec<i32>) {
    let file = File::open(file_path).expect("Unable to open file");
    let reader = io::BufReader::new(file);

    let mut lines = reader.lines();
    while let Some(Ok(line)) = lines.next() {
        if !line.starts_with('#') {
            let l_ele: Vec<String> = line.trim().split_whitespace().map(String::from).collect();
            let mut d_ele2int = HashMap::new();
            let mut d_int2ele = HashMap::new();

            for (i, ele) in l_ele.iter().enumerate() {
                d_ele2int.insert(ele.clone(), i);
                d_ele2int.insert(ele.to_lowercase(), i);
                d_int2ele.insert(i, ele.clone());
            }

            let mut l_score = Vec::new();
            for line in lines {
                let line = line.expect("Unable to read line");
                l_score.extend(line.trim().split_whitespace().skip(1).map(|x| x.parse::<i32>().unwrap()));
            }

            return (l_ele, d_ele2int, d_int2ele, l_score);
        }
    }

    panic!("Matrix file is empty or does not contain valid data");
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
    // Call the FFI function and handle the result
    let res = unsafe { (ssw.ssw_align)(q_profile, r_num.as_ptr() as *const c_uchar, n_r_len, n_open, n_ext, n_flag, 0, 0, n_mask_len) };

    if res.is_null() {
        panic!("Alignment failed: returned null pointer");
    }

    // Dereference the result pointer safely
    let res = unsafe { &*res };

    // Ensure the cigar length is valid before allocation
    let s_cigar_len = res.nCigarLen as usize;
    if s_cigar_len == 0 {
        return (res.nScore, res.nScore2, res.nRefBeg, res.nRefEnd, res.nQryBeg, res.nQryEnd, res.nRefEnd2, 0, Vec::new());
    }

    // Safely collect cigar values
    let l_cigar: Vec<u32> = (0..s_cigar_len).map(|i| unsafe { *res.sCigar.add(i) }).collect();

    // Print score - ensure no issues here
    // println!("Score: {}", res.nScore);

    // Return the tuple
    (
        res.nScore,
        res.nScore2,
        res.nRefBeg,
        res.nRefEnd,
        res.nQryBeg,
        res.nQryEnd,
        res.nRefEnd2,
        s_cigar_len as i32,
        l_cigar
    )
}



fn main() {
    let query = "/mnt/869990e7-a61f-469f-99fe-a48d24ac44ca/git/ebi/query.fa";
    let target = "/mnt/869990e7-a61f-469f-99fe-a48d24ac44ca/git/ebi/gencode.v46.transcripts.200bp.fa";

    // Define constants
    let n_match = 2;
    let n_mismatch = 2;
    let n_open = 3;
    let n_ext = 1;

    // Define elements and mappings
    let l_ele = ['A', 'C', 'G', 'T', 'N'];

    // Use `String` as key type
    let mut d_ele2int = HashMap::new();
    let mut d_int2ele = HashMap::new();

    for (i, &ele) in l_ele.iter().enumerate() {
        d_ele2int.insert(ele.to_string(), i);
        d_ele2int.insert(ele.to_ascii_lowercase().to_string(), i);
        d_int2ele.insert(i, ele);
    }

    // Compute scoring matrix
    let n_ele_num = l_ele.len();
    let mut l_score = vec![0; n_ele_num * n_ele_num];
    for i in 0..n_ele_num {
        for j in 0..n_ele_num {
            l_score[i * n_ele_num + j] = if l_ele[i] == l_ele[j] { n_match } else { -n_mismatch };
        }
    }

    // Convert scoring matrix to i8 type
    let mut mat: Vec<i8> = l_score.into_iter().map(|x| x as i8).collect();

    // Load the SSW library
    let ssw = unsafe { CSsw::new("./").expect("Failed to load SSW library") };

    let mut best_score = u32::MIN;
    // let mut best_alignment = None;
    // let mut ref_id = None;

    let mut best_alignments: HashMap<String, (String, u32)> = HashMap::new();

    for (s_q_id, s_q_seq) in read_fasta(query) {
        let q_num = to_int(&s_q_seq, &d_ele2int, l_ele.len());
        let q_profile = unsafe { (ssw.ssw_init)(
            q_num.as_ptr() as *const c_uchar,
            q_num.len() as c_int,
            mat.as_ptr() as *const c_uchar,
            l_ele.len() as c_int,
            1
        )};

        if q_profile.is_null() {
            eprintln!("Error: q_profile is null.");
            continue;
        }

        let n_mask_len = (s_q_seq.len() / 2) as i32;
        let mut best_score = u32::MIN;
        let mut best_ref_id = String::new();

        for (s_r_id, s_r_seq) in read_fasta(target) {
            let r_num = to_int(&s_r_seq, &d_ele2int, l_ele.len());
            if r_num.is_empty() {
                eprintln!("Error: r_num is empty for target {}", s_r_id);
                continue;
            }

            unsafe {
                let (score, _, _, _, _, _, _, _, _) = align_one(
                    &ssw, q_profile, &r_num, r_num.len() as i32, n_open, n_ext, 0, n_mask_len
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

        unsafe { if !q_profile.is_null() { (ssw.init_destroy)(q_profile); } }

        if best_score != u32::MIN {
            best_alignments.insert(s_q_id, (best_ref_id, best_score));
        }
    }

    for (query_id, (ref_id, score)) in best_alignments {
        println!("Query ID: {}, Best Reference ID: {}, Best Score: {}", query_id, ref_id, score);
    }

    
    
}
