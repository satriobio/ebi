// ssw_lib.rs

use std::ffi::CString;
use std::os::raw::{c_int, c_uint, c_uchar};
use std::ptr;
use libloading::Library;

#[repr(C)]
pub struct CAlignRes {
    pub nScore: c_uint,
    pub nScore2: c_uint,
    pub nRefBeg: c_int,
    pub nRefEnd: c_int,
    pub nQryBeg: c_int,
    pub nQryEnd: c_int,
    pub nRefEnd2: c_int,
    pub sCigar: *mut c_uint,
    pub nCigarLen: c_int,
}

#[repr(C)]
pub struct CProfile {
    pub pByte: *mut c_int,
    pub pWord: *mut c_int,
    pub pRead: *mut c_uchar,
    pub pMat: *mut c_uchar,
    pub nReadLen: c_int,
    pub nN: c_int,
    pub nBias: c_uchar,
}

pub struct CSsw {
    lib: Library,
    ssw_init: libloading::Symbol<'static, unsafe extern fn(*const c_uchar, c_int, *const c_uchar, c_int, c_uchar) -> *mut CProfile>,
    init_destroy: libloading::Symbol<'static, unsafe extern fn(*mut CProfile)>,
    ssw_align: libloading::Symbol<'static, unsafe extern fn(*mut CProfile, *const c_uchar, c_int, c_uchar, c_uchar, c_uchar, c_uint, c_int, c_int) -> *mut CAlignRes>,
    align_destroy: libloading::Symbol<'static, unsafe extern fn(*mut CAlignRes)>,
}

impl CSsw {
    pub fn new(s_lib_path: &str) -> Self {
        let s_lib_name = "libssw.so";
        let lib_path = if !s_lib_path.is_empty() {
            std::path::Path::new(s_lib_path).join(s_lib_name)
        } else {
            let mut found = false;
            let mut lib_path = std::path::Path::new(s_lib_name).to_path_buf();
            for path in std::env::split_paths(&std::env::var("PATH").unwrap_or_default()) {
                if path.join(&lib_path).exists() {
                    lib_path = path.join(s_lib_name);
                    found = true;
                    break;
                }
            }
            if !found {
                panic!("libssw.so does not exist in PATH");
            }
            lib_path
        };

        let lib = Library::new(lib_path).expect("Failed to load library");
        unsafe {
            CSsw {
                lib,
                ssw_init: lib.get(b"ssw_init").unwrap(),
                init_destroy: lib.get(b"init_destroy").unwrap(),
                ssw_align: lib.get(b"ssw_align").unwrap(),
                align_destroy: lib.get(b"align_destroy").unwrap(),
            }
        }
    }

    pub unsafe fn init(&self, read_seq: *const c_uchar, read_len: c_int, mat: *const c_uchar, n: c_int, bias: c_uchar) -> *mut CProfile {
        (self.ssw_init)(read_seq, read_len, mat, n, bias)
    }

    pub unsafe fn destroy_init(&self, profile: *mut CProfile) {
        (self.init_destroy)(profile)
    }

    pub unsafe fn align(&self, profile: *mut CProfile, ref_seq: *const c_uchar, ref_len: c_int, open: c_uchar, ext: c_uchar, flag: c_uchar, score: c_uint, ref_beg: c_int, qry_beg: c_int) -> *mut CAlignRes {
        (self.ssw_align)(profile, ref_seq, ref_len, open, ext, flag, score, ref_beg, qry_beg)
    }

    pub unsafe fn destroy_align(&self, align_res: *mut CAlignRes) {
        (self.align_destroy)(align_res)
    }
}
