use rug::{Float, Assign};
use std::ops::{Add, Mul, Neg, Sub, AddAssign};

// Using rug type signatures, rust book and rust by example throughout
// as guide

// A vector with arbitrary precision values
#[derive(Clone, Debug)]
pub struct APVec3 {
    x: Float,
    y: Float,
    z: Float,
}

impl APVec3 {
    pub fn new(x: f64,
           y: f64,
           z: f64,
           prec: u32
    ) -> APVec3 {
        APVec3 {
            x: Float::with_val(prec, x),
            y: Float::with_val(prec, y),
            z: Float::with_val(prec, z),
        }
    }

    pub fn assign(&mut self, vec: &APVec3) {
        self.x.assign(&vec.x);
        self.y.assign(&vec.y);
        self.z.assign(&vec.z);
    }

    pub fn scalar_mul_assign(&mut self, scalar: i32, vec: &APVec3) {
        self.x.assign(scalar * &vec.x);
        self.y.assign(scalar * &vec.y);
        self.z.assign(scalar * &vec.z);
    }

    pub fn parse(x: &str, y: &str, z: &str, prec: u32) -> APVec3 {
        APVec3 {
            x: Float::with_val(prec, Float::parse(x).unwrap()),
            y: Float::with_val(prec, Float::parse(y).unwrap()),
            z: Float::with_val(prec, Float::parse(z).unwrap()),
        }
    }

    pub fn norm(&self) -> Float {
        let prec = self.x.prec();
        let mut result = Float::with_val(prec, &self.x * &self.x + &self.y * &self.y);
        result += Float::with_val(prec, &self.z * &self.z);
        result.sqrt_mut();
        result
    }

    pub fn norm_mut(&self, result: &mut Float) {
        result.assign(&self.x * &self.x + &self.y * &self.y);
        result.add_assign(&self.z * &self.z);
        result.sqrt_mut();
    }

}

// All overloaded operators require ownership of at
// least one allocation. All implemented
// operations commute (w.r.t types, subtraction
// obviously does not commute with result)
//
// Currently implemented:
//
// Add for: 
// both owned
// 1 owned, one immut ref
//
// Mul for:
// i32 and owned
//
// Neg for:
// owned
//
// Sub for:
// both owned
// 1 owned, one immut ref

impl Add for APVec3 {
    type Output = APVec3;

    fn add(mut self, rhs: APVec3) -> APVec3 {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
        self
    }
}

// Used operator overload in rug as an example
// for type signatures
// also see
// https://stackoverflow.com/questions/28005134/how-do-i-implement-the-add-trait-for-a-reference-to-a-struct
impl<'a> Add<&'a APVec3> for APVec3 {
    type Output = APVec3;

    fn add(mut self, rhs: &APVec3) -> APVec3 {
        self.x += &rhs.x;
        self.y += &rhs.y;
        self.z += &rhs.z;
        self
    }
}

impl<'a> Add<APVec3> for &'a APVec3 {
    type Output = APVec3;

    fn add(self, mut rhs: APVec3) -> APVec3 {
        rhs.x += &self.x;
        rhs.y += &self.y;
        rhs.z += &self.z;
        rhs
    }
}

impl Mul<i32> for APVec3 {
    type Output = APVec3;

    fn mul(mut self, rhs: i32) -> APVec3 {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
        self
    }
}

impl Mul<APVec3> for i32 {
    type Output = APVec3;

    fn mul(self, mut rhs: APVec3) -> APVec3 {
        rhs.x *= self;
        rhs.y *= self;
        rhs.z *= self;
        rhs
    }
}

impl Neg for APVec3 {
    type Output = APVec3;

    fn neg(mut self) -> APVec3 {
        self.x = -self.x;
        self.y = -self.y;
        self.z = -self.z;
        self
    }
}

impl Sub<APVec3> for APVec3 {
    type Output = APVec3;

    fn sub(mut self, rhs: APVec3) -> APVec3 {
        self.x -= rhs.x;
        self.y -= rhs.y;
        self.z -= rhs.z;
        self
    }
}

impl<'a> Sub<&'a APVec3> for APVec3 {
    type Output = APVec3;

    fn sub(mut self, rhs: &APVec3) -> APVec3 {
        self.x -= &rhs.x;
        self.y -= &rhs.y;
        self.z -= &rhs.z;
        self
    }
}

impl<'a> Sub<APVec3> for &'a APVec3 {
    type Output = APVec3;

    fn sub(self, mut rhs: APVec3) -> APVec3 {
        rhs.x = &self.x - rhs.x;
        rhs.y = &self.y - rhs.y;
        rhs.z = &self.z - rhs.z;
        rhs
    }
}

impl<'a> AddAssign<&'a APVec3> for APVec3 {
    fn add_assign(&mut self, rhs: &APVec3) {
        self.x += &rhs.x;
        self.y += &rhs.y;
        self.z += &rhs.z;
    }
}
