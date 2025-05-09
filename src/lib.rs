#[macro_use]
extern crate structopt;
#[macro_use]
extern crate lazy_static;
extern crate nalgebra as na;
extern crate rayon;
extern crate rug;

mod geometry;

use std::collections::VecDeque;
use std::path::PathBuf;
use std::io::prelude::*;
use std::fs::File;
use na::Vector3;
use rayon::prelude::*;
use std::f64::consts::PI;
use rug::{Float, Assign};
use rug::ops::Pow;
use geometry::APVec3;

// Constant Declarations for whole crate
const THETA_BIN_SIZE: f64 = 0.01;

//From Mathematica 11.2.0
const PI_STR: &str = "3.141592653589793238462643383279502884197";

// Constants for doing lattice sums. Non-obvious
// numerical values computed through Mathematica 11.2.0
const Z_OFFSET: f64 = 0.816496580927726;

const X_OFFSET_B: f64 = 0.5;
const Y_OFFSET_B: f64 = 0.2886751345948129;

const X_OFFSET_C: f64 = -0.5;
const Y_OFFSET_C: f64 = -0.2886751345948129;

const X_HEX_1: f64 = 1.0;
const Y_HEX_1: f64 = 0.0;
const X_HEX_2: f64 = 0.5;
const Y_HEX_2: f64 = 0.866025403784439;

const AP_PREC: u32 = 128; // precision for arbitrary precision calculation

// Constants for doing lattice sums, with 40 digit precision. Non-obvious
// numerical values computed through Mathematica 11.2.0
const AP_Z_OFFSET: &str = "0.8164965809277260327324280249019637973220";

const AP_X_OFFSET_B: &str = "0.5";
const AP_Y_OFFSET_B: &str = "0.2886751345948128822545743902509787278238";

const AP_X_OFFSET_C: &str = "-0.5";
const AP_Y_OFFSET_C: &str = "-0.2886751345948128822545743902509787278238";

const AP_X_HEX_1: &str = "1.0";
const AP_Y_HEX_1: &str = "0.0";
const AP_X_HEX_2: &str = "0.5";
const AP_Y_HEX_2: &str = "0.8660254037844386467637231707529361834714";

// lazy_static declarations of useful vectors
lazy_static! {
    static ref OFFSET_ALL: Vector3<f64> = Vector3::new(0.0, 0.0, Z_OFFSET);
    static ref OFFSET_B: Vector3<f64> = Vector3::new(X_OFFSET_B, Y_OFFSET_B, 0.0);
    static ref OFFSET_C: Vector3<f64> = Vector3::new(X_OFFSET_C, Y_OFFSET_C, 0.0);
    static ref HEX_1: Vector3<f64> = Vector3::new(X_HEX_1, Y_HEX_1, 0.0);
    static ref HEX_2: Vector3<f64> = Vector3::new(X_HEX_2, Y_HEX_2, 0.0);

    // Used in arbitrary precision versions of methods
    static ref AP_OFFSET_ALL: APVec3 = APVec3::parse("0.0", "0.0", AP_Z_OFFSET, AP_PREC);
    static ref AP_OFFSET_B: APVec3 = APVec3::parse(AP_X_OFFSET_B, AP_Y_OFFSET_B, "0.0", AP_PREC);
    static ref AP_OFFSET_C: APVec3 = APVec3::parse(AP_X_OFFSET_C, AP_Y_OFFSET_C, "0.0", AP_PREC);
    static ref AP_HEX_1: APVec3 = APVec3::parse(AP_X_HEX_1, AP_Y_HEX_1, "0.0", AP_PREC);
    static ref AP_HEX_2: APVec3 = APVec3::parse(AP_X_HEX_2, AP_Y_HEX_2, "0.0", AP_PREC);
}

// StructOpt derivation is based on Github sample code for structopt
///A program to calculate order metrics for periodic Barlow packings.
#[derive(StructOpt, Debug)]
#[structopt(name = "Periodic Barlow Order Metric Calculator")]
pub struct Opt {
    ///Choose order metric.
    ///Panics on invalid request.
    ///Choices:
    ///theta (mostly for debugging, preprogrammed bin size),
    ///T,
    ///T*,
    ///tau,
    ///hyperuniformity
    #[structopt(short = "m", long = "metric")]
    pub metric: String,

    ///Distance cutoff, in units of particle diameter
    #[structopt(short = "c", long = "cutoff")]
    pub cutoff: f64,

    ///Input file. Format:
    ///ABCABC\nABABAB
    #[structopt(short = "i", long = "input", parse(from_os_str))]
    pub input: PathBuf,

    ///Output file (optional)
    #[structopt(short = "o", long = "output", parse(from_os_str))]
    pub output: Option<PathBuf>,

    ///Convergence parameter for hyperuniformity B coefficient
    #[structopt(short = "a", long = "alpha")]
    pub alpha: Option<f64>,

    ///Convergence parameter for AP h integral and B coeff
    #[structopt(short= "s", long = "alphastr")]
    pub alpha_str: Option<String>,

    ///Verbosity level, useful for knowing run time
    #[structopt(short = "v", long = "verbose")]
    pub verbose: bool,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Layer {
    A,
    B,
    C,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Spec ( VecDeque<Layer> );

// Spec invariants (may panic if not upheld):
// Non-empty
// Starts with A
// Ends with B or C
impl Spec {
    fn len(&self) -> usize {
        self.0.len()
    }

    fn get(&self, n: usize) -> Layer {
        self.0[n]
    }

    pub fn new_with_check(raw_spec: &str) -> Spec {
        use Layer::*;

        // Parse the string read in from file
        let spec: Spec = Spec ( raw_spec.chars()
            .map(|c| match c {
                'A' => A,
                'B' => B,
                'C' => C,
                _ => panic!("Invalid spec string"), }
            )
            .collect() );
        
        // Check invariants
        assert!(*spec.0.front().expect("Spec must be non-empty") == A);
        assert!(*spec.0.back().expect("Spec must be non-empty") == B 
                || *spec.0.back().expect("Spec must be non-empty") == C);
        
        spec
    }
    
    fn shift(&self, n: usize) -> Spec {
        use Layer::*;
        
        // Rotate list of layers
        let mut new_spec = self.clone();
        for _ in 0..n {
            let front = new_spec.0.pop_front().expect("Spec can't be empty.");
            new_spec.0.push_back(front);
        }
        
        // Rename layers so beginning is still A
        match *new_spec.0.front().expect("Spec must be non-empty") {
            A => (),
            B => {
                for layer in &mut new_spec.0 {
                    *layer = match *layer {
                        A => C,
                        B => A,
                        C => B,
                    };
                }
            },
            C => {
                for layer in &mut new_spec.0 {
                    *layer = match *layer {
                        A => B,
                        B => C,
                        C => A,
                    };
                }
            },
        }
        new_spec
    }
}

pub fn theta_hex(
    layer: Layer, // A B or C
    l_idx: i64, //How many layers above or below center are we?
    max_lat_idx: i64, //How far are we going out in lattice sum?
    theta: &mut Vec<f64>) //Holder for theta series
{
    use Layer::*;
    let l_f64 = l_idx as f64;
    
    let offset = match layer {
        A => l_f64 * *OFFSET_ALL,
        B => l_f64 * *OFFSET_ALL + *OFFSET_B,
        C => l_f64 * *OFFSET_ALL + *OFFSET_C,
    };

    for i in (-max_lat_idx)..(max_lat_idx + 1) {
        for j in (-max_lat_idx)..(max_lat_idx + 1) {
            let vec = offset + (i as f64) * *HEX_1 + (j as f64) * *HEX_2;
            let r: f64 = na::norm(&vec);
            let idx = (r/THETA_BIN_SIZE) as usize;
            if idx < theta.len() {
                theta[idx] += 1.0;
            }
        }
    }

}

pub fn compute_theta(spec: Spec, opt: &Opt) {
    // turn radial cutoff into integer lattice cutoff
    // probably overkill
    let max_lat_idx = (2.0*opt.cutoff) as i64;
    
    let spec_len = spec.len();
    let spec_i64 = spec_len as i64;

    // set up theta series bins
    let n_bins = (opt.cutoff/THETA_BIN_SIZE) as usize;
    let domain: Vec<f64> = (0..n_bins)
        .map(|x| (x as f64)*THETA_BIN_SIZE)
        .collect();
    let mut theta: Vec<f64> = vec![0.0; n_bins];
    
    // iterate over particle types
    for i in 0..spec_len {
        println!("Working on particle type: {}", i);
        let curr_spec = spec.shift(i);
        // iterate over layers
        for j in (-max_lat_idx)..(max_lat_idx + 1) {
            //println!("Working on layer: {}", j);
            // standard modulo hack, see for example
            // https://stackoverflow.com/questions/31210357/is-there-a-modulus-not-remainder-function-operation
            let layer = curr_spec.get((((j%spec_i64)+spec_i64)%spec_i64) as usize);
            theta_hex(layer, j, max_lat_idx, &mut theta);
        }
    }

    // Output theta series, optionally to file
    let mut output_file;
    if let Some(ref output_path) = opt.output {
        output_file = Some(File::create(output_path)
                           .expect("Error creating output file"));
    } else {
        output_file = None;
    }
    for (i,t) in theta.iter_mut().enumerate() {
        *t /= spec_len as f64;
        println!("{}, {}", domain[i], *t);
        match output_file {
            Some(ref mut buffer) => {
                write!(buffer, "{}, {}\n", domain[i], *t)
                    .expect("Error writing to output file");
            },
            None => (),
        }
    }
}

pub fn hyper_hex(
    layer_pos: Layer, // A B or C
    layer_neg: Layer,
    l_idx: i64, //How many layers above or below center are we?
    max_lat_idx: i64, //How far are we going out in lattice sum?
    alpha: f64 //exponential convergence parameter
    ) -> f64
{
    use Layer::*;
    let l_f64 = l_idx as f64;
    
    let offset_pos = match layer_pos {
        A => l_f64 * *OFFSET_ALL,
        B => l_f64 * *OFFSET_ALL + *OFFSET_B,
        C => l_f64 * *OFFSET_ALL + *OFFSET_C,
    };
    
    let offset_neg = match layer_neg {
        A => -l_f64 * *OFFSET_ALL,
        B => -l_f64 * *OFFSET_ALL + *OFFSET_B,
        C => -l_f64 * *OFFSET_ALL + *OFFSET_C,
    };

    let mut b: f64 = 0.0;

    for i in (0..(max_lat_idx + 1)) {
        for j in (0..(max_lat_idx + 1)) {
            // +++
            let vec = offset_pos + (i as f64) * *HEX_1 + (j as f64) * *HEX_2;
            let r: f64 = na::norm(&vec);
            b += r*((-alpha*r*r).exp());
            // ++-
            if j != 0 {
                let vec = offset_pos + (i as f64) * *HEX_1 - (j as f64) * *HEX_2;
                let r: f64 = na::norm(&vec);
                b += r*((-alpha*r*r).exp());
            }
            // +-+
            if i != 0 {
                let vec = offset_pos - (i as f64) * *HEX_1 + (j as f64) * *HEX_2;
                let r: f64 = na::norm(&vec);
                b += r*((-alpha*r*r).exp());
            }
            // +--
            if (i != 0) && (j != 0) {
                let vec = offset_pos - (i as f64) * *HEX_1 - (j as f64) * *HEX_2;
                let r: f64 = na::norm(&vec);
                b += r*((-alpha*r*r).exp());
            }
            // -++
            if l_idx != 0 {
                let vec = offset_neg + (i as f64) * *HEX_1 + (j as f64) * *HEX_2;
                let r: f64 = na::norm(&vec);
                b += r*((-alpha*r*r).exp());
            }
            // -+-
            if (l_idx != 0) && (j != 0) {
                let vec = offset_neg + (i as f64) * *HEX_1 - (j as f64) * *HEX_2;
                let r: f64 = na::norm(&vec);
                b += r*((-alpha*r*r).exp());
            }
            // --+
            if (l_idx != 0) && (i != 0) {
                let vec = offset_neg - (i as f64) * *HEX_1 + (j as f64) * *HEX_2;
                let r: f64 = na::norm(&vec);
                b += r*((-alpha*r*r).exp());
            }
            // ---
            if (l_idx != 0) && (i != 0) && (j != 0) {
                let vec = offset_neg - (i as f64) * *HEX_1 - (j as f64) * *HEX_2;
                let r: f64 = na::norm(&vec);
                b += r*((-alpha*r*r).exp());
            }
        }
    }
    b
}


pub fn compute_hyperuniformity(spec: Spec, opt: &Opt) {
    let alpha = opt.alpha
        .expect("Must supply convergence parameter alpha for B calculation.");

    let max_lat_idx = (2.0*opt.cutoff) as i64;

    let spec_len = spec.len();
    let spec_i64 = spec_len as i64;

    let mut b: f64 = 0.0; //B surficial variance coefficient
    
    for i in 0..spec_len {
        println!("Working on particle type: {}", i);
        let curr_spec = spec.shift(i);
        let temp_vec: Vec<f64> = (0..(max_lat_idx + 1))
            .collect::<Vec<i64>>()
            .par_iter()
            .map(|&j| hyper_hex(
                curr_spec.get((((j%spec_i64)+spec_i64)%spec_i64) as usize),
                curr_spec.get((((-j%spec_i64)+spec_i64)%spec_i64) as usize),
                j,
                max_lat_idx,
                alpha
                )
            )
            .collect();
            //println!("{:?}", temp_vec);
            //.sum::<f64>();
        b += temp_vec.iter().sum::<f64>();
    }

    println!("b is {}", b);
    
    // normalize by # of particle types
    b /= spec_len as f64;
    
    // account for the -1 in the definition of h
    b -= 4.0*PI*(2.0f64.sqrt())/(2.0*alpha*alpha);

    println!("B Coefficient (non-normalized) = {}", b);
}

pub fn h_hex(
    layer_pos: Layer, // A B or C
    layer_neg: Layer,
    l_idx: i64, //How many layers above or below center are we?
    max_lat_idx: i64, //How far are we going out in lattice sum?
    alpha: f64 //exponential convergence parameter
    ) -> f64
{
    use Layer::*;
    let l_f64 = l_idx as f64;
    
    let offset_pos = match layer_pos {
        A => l_f64 * *OFFSET_ALL,
        B => l_f64 * *OFFSET_ALL + *OFFSET_B,
        C => l_f64 * *OFFSET_ALL + *OFFSET_C,
    };
    
    let offset_neg = match layer_neg {
        A => -l_f64 * *OFFSET_ALL,
        B => -l_f64 * *OFFSET_ALL + *OFFSET_B,
        C => -l_f64 * *OFFSET_ALL + *OFFSET_C,
    };

    let mut h: f64 = 0.0;

    for i in (0..(max_lat_idx + 1)) {
        for j in (0..(max_lat_idx + 1)) {
            // +++
            let vec = offset_pos + (i as f64) * *HEX_1 + (j as f64) * *HEX_2;
            let r: f64 = na::norm(&vec);
            h += (-alpha*r*r).exp();
            // ++-
            if j != 0 {
                let vec = offset_pos + (i as f64) * *HEX_1 - (j as f64) * *HEX_2;
                let r: f64 = na::norm(&vec);
                h += (-alpha*r*r).exp();
            }
            // +-+
            if i != 0 {
                let vec = offset_pos - (i as f64) * *HEX_1 + (j as f64) * *HEX_2;
                let r: f64 = na::norm(&vec);
                h += (-alpha*r*r).exp();
            }
            // +--
            if (i != 0) && (j != 0) {
                let vec = offset_pos - (i as f64) * *HEX_1 - (j as f64) * *HEX_2;
                let r: f64 = na::norm(&vec);
                h += (-alpha*r*r).exp();
            }
            // -++
            if l_idx != 0 {
                let vec = offset_neg + (i as f64) * *HEX_1 + (j as f64) * *HEX_2;
                let r: f64 = na::norm(&vec);
                h += (-alpha*r*r).exp();
            }
            // -+-
            if (l_idx != 0) && (j != 0) {
                let vec = offset_neg + (i as f64) * *HEX_1 - (j as f64) * *HEX_2;
                let r: f64 = na::norm(&vec);
                h += (-alpha*r*r).exp();
            }
            // --+
            if (l_idx != 0) && (i != 0) {
                let vec = offset_neg - (i as f64) * *HEX_1 + (j as f64) * *HEX_2;
                let r: f64 = na::norm(&vec);
                h += (-alpha*r*r).exp();
            }
            // ---
            if (l_idx != 0) && (i != 0) && (j != 0) {
                let vec = offset_neg - (i as f64) * *HEX_1 - (j as f64) * *HEX_2;
                let r: f64 = na::norm(&vec);
                h += (-alpha*r*r).exp();
            }
        }
    }
    h
}


pub fn compute_h(spec: Spec, opt: &Opt) {
    let alpha = opt.alpha
        .expect("Must supply convergence parameter alpha for h integral calculation.");

    let max_lat_idx = (2.0*opt.cutoff) as i64;

    let spec_len = spec.len();
    let spec_i64 = spec_len as i64;

    let mut h: f64 = 0.0; //integral of h over R3

    for i in 0..spec_len {
        println!("Working on particle type: {}", i);
        let curr_spec = spec.shift(i);
        let temp: Vec<f64> = (0..(max_lat_idx + 1))
            .collect::<Vec<i64>>()
            .par_iter()
            .map(|&j| h_hex(
                curr_spec.get((((j%spec_i64)+spec_i64)%spec_i64) as usize),
                curr_spec.get((((-j%spec_i64)+spec_i64)%spec_i64) as usize),
                j,
                max_lat_idx,
                alpha
                )
            )
            .collect();
        h += temp.iter().sum::<f64>();
    }

    println!("h is {}", h);
    
    // normalize by # of particle types
    h /= spec_len as f64;
   
    // must subtract off origin contribution, not zeroed by inclusion
    // of r anymore
    h -= 1.0;

    // account for the -1 in the definition of h
    h -= 2.0f64.sqrt()*PI.powf(1.5)/alpha.powf(1.5);

    println!("h = {}", h);

}

pub fn ap_h_hex(
    layer_pos: Layer, // A B or C
    layer_neg: Layer,
    l_idx: i32, //How many layers above or below center are we?
    max_lat_idx: i32, //How far are we going out in lattice sum?
    alpha: &Float, //exponential convergence parameter
    opt: &Opt
    ) -> Float
{
    use Layer::*;
    
    if opt.verbose {
        println!("Working on l_idx: {}", l_idx);
    }
    
    let offset_pos = match layer_pos {
        A => l_idx * AP_OFFSET_ALL.clone(),
        B => l_idx * AP_OFFSET_ALL.clone() + &*AP_OFFSET_B,
        C => l_idx * AP_OFFSET_ALL.clone() + &*AP_OFFSET_C,
    };
    
    let offset_neg = match layer_neg {
        A => -l_idx * AP_OFFSET_ALL.clone(),
        B => -l_idx * AP_OFFSET_ALL.clone() + &*AP_OFFSET_B,
        C => -l_idx * AP_OFFSET_ALL.clone() + &*AP_OFFSET_C,
    };

    let mut h: Float = Float::with_val(AP_PREC, 0);
    let mut r: Float = Float::with_val(AP_PREC, 0);
    let mut vec = APVec3::new(0.0, 0.0, 0.0, AP_PREC);
    let mut buffer = APVec3::new(0.0, 0.0, 0.0, AP_PREC);

    for i in (0..(max_lat_idx + 1)) {
        for j in (0..(max_lat_idx + 1)) {
            // +++
            vec.assign(&offset_pos);
            buffer.scalar_mul_assign(i, &*AP_HEX_1);
            vec += &buffer; 
            buffer.scalar_mul_assign(j, &*AP_HEX_2);
            vec += &buffer;
            vec.norm_mut(&mut r);
            r.square_mut();
            r = (-(alpha*r)).exp();
            h += &r;
            
            // ++-
            if j != 0 {
                vec.assign(&offset_pos);
                buffer.scalar_mul_assign(i, &*AP_HEX_1);
                vec += &buffer; 
                buffer.scalar_mul_assign(-j, &*AP_HEX_2);
                vec += &buffer;
                vec.norm_mut(&mut r);
                r.square_mut();
                r = (-(alpha*r)).exp();
                h += &r;
            }
            
            // +-+
            if i != 0 {
                vec.assign(&offset_pos);
                buffer.scalar_mul_assign(-i, &*AP_HEX_1);
                vec += &buffer; 
                buffer.scalar_mul_assign(j, &*AP_HEX_2);
                vec += &buffer;
                vec.norm_mut(&mut r);
                r.square_mut();
                r = (-(alpha*r)).exp();
                h += &r;
            }
            
            // +--
            if (i != 0) && (j != 0) {
                vec.assign(&offset_pos);
                buffer.scalar_mul_assign(-i, &*AP_HEX_1);
                vec += &buffer; 
                buffer.scalar_mul_assign(-j, &*AP_HEX_2);
                vec += &buffer;
                vec.norm_mut(&mut r);
                r.square_mut();
                r = (-(alpha*r)).exp();
                h += &r;
            }
            
            // -++
            if l_idx != 0 {
                vec.assign(&offset_neg);
                buffer.scalar_mul_assign(i, &*AP_HEX_1);
                vec += &buffer; 
                buffer.scalar_mul_assign(j, &*AP_HEX_2);
                vec += &buffer;
                vec.norm_mut(&mut r);
                r.square_mut();
                r = (-(alpha*r)).exp();
                h += &r;
            }
            
            // -+-
            if (l_idx != 0) && (j != 0) {
                vec.assign(&offset_neg);
                buffer.scalar_mul_assign(i, &*AP_HEX_1);
                vec += &buffer; 
                buffer.scalar_mul_assign(-j, &*AP_HEX_2);
                vec += &buffer;
                vec.norm_mut(&mut r);
                r.square_mut();
                r = (-(alpha*r)).exp();
                h += &r;
            }
            
            // --+
            if (l_idx != 0) && (i != 0) {
                vec.assign(&offset_neg);
                buffer.scalar_mul_assign(-i, &*AP_HEX_1);
                vec += &buffer; 
                buffer.scalar_mul_assign(j, &*AP_HEX_2);
                vec += &buffer;
                vec.norm_mut(&mut r);
                r.square_mut();
                r = (-(alpha*r)).exp();
                h += &r;
            }
            
            // ---
            if (l_idx != 0) && (i != 0) && (j != 0) {
                vec.assign(&offset_neg);
                buffer.scalar_mul_assign(-i, &*AP_HEX_1);
                vec += &buffer; 
                buffer.scalar_mul_assign(-j, &*AP_HEX_2);
                vec += &buffer;
                vec.norm_mut(&mut r);
                r.square_mut();
                r = (-(alpha*r)).exp();
                h += &r;
            }
        }
    }
    h
}


pub fn ap_compute_h(spec: Spec, opt: &Opt) {
    // this passes borrow checker thanks to
    // https://users.rust-lang.org/t/cannot-move-out-of-borrowed-context-access-struct-field/1700
    let alpha_str = opt.alpha_str.as_ref().expect("Must supply convergence parameter");
    let alpha = Float::with_val(AP_PREC, 
                                Float::parse(alpha_str)
                                .expect("Invalid float spec for alpha_str"));

    let max_lat_idx = (2.0*opt.cutoff) as i32;

    let spec_len = spec.len();
    let spec_i32 = spec_len as i32;

    let mut h: Float = Float::with_val(AP_PREC, 0); //integral of h over R3

    for i in 0..spec_len {
        println!("Working on particle type: {}", i);
        let curr_spec = spec.shift(i);
        let temp: Vec<Float> = (0..(max_lat_idx + 1))
            .collect::<Vec<i32>>()
            .par_iter()
            .map(|&j| ap_h_hex(
                curr_spec.get((((j%spec_i32)+spec_i32)%spec_i32) as usize),
                curr_spec.get((((-j%spec_i32)+spec_i32)%spec_i32) as usize),
                j,
                max_lat_idx,
                &alpha,
                opt
                )
            )
            .collect();
        h += Float::with_val(AP_PREC, Float::sum(temp.iter()));
    }

    println!("h is {}", h);
    
    // normalize by # of particle types
    h /= Float::with_val(AP_PREC, spec_len);
   
    // must subtract off origin contribution, not zeroed by inclusion
    // of r anymore
    h -= Float::with_val(AP_PREC, 1);

    // account for the -1 in the definition of h
    h -= Float::with_val(AP_PREC, 2).sqrt()
        *Float::with_val(AP_PREC, Float::parse(PI_STR).unwrap())
        .pow(Float::with_val(AP_PREC, Float::parse("1.5").unwrap()))
        /alpha.pow(Float::with_val(AP_PREC, Float::parse("1.5").unwrap()));

    println!("h = {}", h);

}

pub fn ap_hyper_hex(
    layer_pos: Layer, // A B or C
    layer_neg: Layer,
    l_idx: i32, //How many layers above or below center are we?
    max_lat_idx: i32, //How far are we going out in lattice sum?
    alpha: &Float, //exponential convergence parameter
    opt: &Opt
    ) -> Float
{
    use Layer::*;
    
    if opt.verbose {
        println!("Working on l_idx: {}", l_idx);
    }

    let offset_pos = match layer_pos {
        A => l_idx * AP_OFFSET_ALL.clone(),
        B => l_idx * AP_OFFSET_ALL.clone() + &*AP_OFFSET_B,
        C => l_idx * AP_OFFSET_ALL.clone() + &*AP_OFFSET_C,
    };
    
    let offset_neg = match layer_neg {
        A => -l_idx * AP_OFFSET_ALL.clone(),
        B => -l_idx * AP_OFFSET_ALL.clone() + &*AP_OFFSET_B,
        C => -l_idx * AP_OFFSET_ALL.clone() + &*AP_OFFSET_C,
    };

    let mut h: Float = Float::with_val(AP_PREC, 0); //buffer for collecting integral contributions
    let mut r: Float = Float::with_val(AP_PREC, 0); //buffer for radius calc
    let mut ans: Float = Float::with_val(AP_PREC, 0); //buffer for composing contribution of layer to int
    let mut vec = APVec3::new(0.0, 0.0, 0.0, AP_PREC); //buffer for building lattice vector
    let mut buffer = APVec3::new(0.0, 0.0, 0.0, AP_PREC); //buffer for scalar multiple calcs

    for i in (0..(max_lat_idx + 1)) {
        for j in (0..(max_lat_idx + 1)) {
            // +++
            vec.assign(&offset_pos);
            buffer.scalar_mul_assign(i, &*AP_HEX_1);
            vec += &buffer; 
            buffer.scalar_mul_assign(j, &*AP_HEX_2);
            vec += &buffer;
            vec.norm_mut(&mut r);
            ans.assign(r.square_ref());
            ans = (&r)*(-(alpha*ans)).exp();
            h += &ans;
            
            // ++-
            if j != 0 {
                vec.assign(&offset_pos);
                buffer.scalar_mul_assign(i, &*AP_HEX_1);
                vec += &buffer; 
                buffer.scalar_mul_assign(-j, &*AP_HEX_2);
                vec += &buffer;
                vec.norm_mut(&mut r);
                ans.assign(r.square_ref());
                ans = (&r)*(-(alpha*ans)).exp();
                h += &ans;
            }
            
            // +-+
            if i != 0 {
                vec.assign(&offset_pos);
                buffer.scalar_mul_assign(-i, &*AP_HEX_1);
                vec += &buffer; 
                buffer.scalar_mul_assign(j, &*AP_HEX_2);
                vec += &buffer;
                vec.norm_mut(&mut r);
                ans.assign(r.square_ref());
                ans = (&r)*(-(alpha*ans)).exp();
                h += &ans;
            }
            
            // +--
            if (i != 0) && (j != 0) {
                vec.assign(&offset_pos);
                buffer.scalar_mul_assign(-i, &*AP_HEX_1);
                vec += &buffer; 
                buffer.scalar_mul_assign(-j, &*AP_HEX_2);
                vec += &buffer;
                vec.norm_mut(&mut r);
                ans.assign(r.square_ref());
                ans = (&r)*(-(alpha*ans)).exp();
                h += &ans;
            }
            
            // -++
            if l_idx != 0 {
                vec.assign(&offset_neg);
                buffer.scalar_mul_assign(i, &*AP_HEX_1);
                vec += &buffer; 
                buffer.scalar_mul_assign(j, &*AP_HEX_2);
                vec += &buffer;
                vec.norm_mut(&mut r);
                ans.assign(r.square_ref());
                ans = (&r)*(-(alpha*ans)).exp();
                h += &ans;
            }
            
            // -+-
            if (l_idx != 0) && (j != 0) {
                vec.assign(&offset_neg);
                buffer.scalar_mul_assign(i, &*AP_HEX_1);
                vec += &buffer; 
                buffer.scalar_mul_assign(-j, &*AP_HEX_2);
                vec += &buffer;
                vec.norm_mut(&mut r);
                ans.assign(r.square_ref());
                ans = (&r)*(-(alpha*ans)).exp();
                h += &ans;
            }
            
            // --+
            if (l_idx != 0) && (i != 0) {
                vec.assign(&offset_neg);
                buffer.scalar_mul_assign(-i, &*AP_HEX_1);
                vec += &buffer; 
                buffer.scalar_mul_assign(j, &*AP_HEX_2);
                vec += &buffer;
                vec.norm_mut(&mut r);
                ans.assign(r.square_ref());
                ans = (&r)*(-(alpha*ans)).exp();
                h += &ans;
            }
            
            // ---
            if (l_idx != 0) && (i != 0) && (j != 0) {
                vec.assign(&offset_neg);
                buffer.scalar_mul_assign(-i, &*AP_HEX_1);
                vec += &buffer; 
                buffer.scalar_mul_assign(-j, &*AP_HEX_2);
                vec += &buffer;
                vec.norm_mut(&mut r);
                ans.assign(r.square_ref());
                ans = (&r)*(-(alpha*ans)).exp();
                h += &ans;
            }
        }
    }
    h
}


pub fn ap_compute_hyper(spec: Spec, opt: &Opt) {
    // this passes borrow checker thanks to
    // https://users.rust-lang.org/t/cannot-move-out-of-borrowed-context-access-struct-field/1700
    let alpha_str = opt.alpha_str.as_ref().expect("Must supply convergence parameter");
    let alpha = Float::with_val(AP_PREC, 
                                Float::parse(alpha_str)
                                .expect("Invalid float spec for alpha_str"));

    let max_lat_idx = (2.0*opt.cutoff) as i32;

    let spec_len = spec.len();
    let spec_i32 = spec_len as i32;

    let mut h: Float = Float::with_val(AP_PREC, 0); //integral of h over R3

    for i in 0..spec_len {
        println!("Working on particle type: {}", i);
        let curr_spec = spec.shift(i);
        let temp: Vec<Float> = (0..(max_lat_idx + 1))
            .collect::<Vec<i32>>()
            .par_iter()
            .map(|&j| ap_hyper_hex(
                curr_spec.get((((j%spec_i32)+spec_i32)%spec_i32) as usize),
                curr_spec.get((((-j%spec_i32)+spec_i32)%spec_i32) as usize),
                j,
                max_lat_idx,
                &alpha,
                opt
                )
            )
            .collect();
        h += Float::with_val(AP_PREC, Float::sum(temp.iter()));
    }

    println!("h is {}", h);
    
    // normalize by # of particle types
    h /= Float::with_val(AP_PREC, spec_len);

    // account for the -1 in the definition of h
    h -= Float::with_val(AP_PREC, 4)
        *Float::with_val(AP_PREC, Float::parse(PI_STR).unwrap())
        *(Float::with_val(AP_PREC, 2).sqrt())
        /(Float::with_val(AP_PREC, 2)*alpha.clone()*alpha.clone());


    println!("h = {}", h);

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn valid_parse() {
        use Layer::*;
        let mut deque = VecDeque::new();
        deque.push_back(A);
        deque.push_back(B);
        deque.push_back(C);
        deque.push_back(B);
        let output = Spec ( deque );
        assert!(Spec::new_with_check("ABCB") == output);
    }
    
    #[test]
    #[should_panic]
    fn invalid_parse1() {
        Spec::new_with_check("ABCA");
    }

    #[test]
    #[should_panic]
    fn invalid_parse2() {
        Spec::new_with_check("BCA");
    }

    #[test]
    fn test_shift_2() {
        assert!(Spec::new_with_check("ABCAB").shift(2) == Spec::new_with_check("ABCBC"));
    }
}


