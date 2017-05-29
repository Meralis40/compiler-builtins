use core::mem;

pub mod conv;
pub mod add;
pub mod pow;
pub mod sub;
pub mod mul;
pub mod div;

/// Trait for some basic operations on floats
pub trait Float: Sized + Copy {
    /// A uint of the same with as the float
    type Int;

    /// Returns the bitwidth of the float type
    fn bits() -> u32;

    /// Returns the bitwidth of the significand
    fn significand_bits() -> u32;

    /// Returns the bitwidth of the exponent
    fn exponent_bits() -> u32 {
        Self::bits() - Self::significand_bits() - 1
    }
    /// Returns the maximum value of the exponent
    fn exponent_max() -> u32 {
        (1 << Self::exponent_bits()) - 1
    }

    /// Returns the exponent bias value
    fn exponent_bias() -> u32 {
        Self::exponent_max() >> 1
    }

    /// Returns a mask for the sign bit
    fn sign_mask() -> Self::Int;

    /// Returns a mask for the significand
    fn significand_mask() -> Self::Int;

    // Returns the implicit bit of the float format
    fn implicit_bit() -> Self::Int;

    /// Returns a mask for the exponent
    fn exponent_mask() -> Self::Int;

    /// Returns `self` transmuted to `Self::Int`
    fn repr(self) -> Self::Int;

    #[cfg(test)]
    /// Checks if two floats have the same bit representation. *Except* for NaNs! NaN can be
    /// represented in multiple different ways. This method returns `true` if two NaNs are
    /// compared.
    fn eq_repr(self, rhs: Self) -> bool;

    /// Returns a `Self::Int` transmuted back to `Self`
    fn from_repr(a: Self::Int) -> Self;

    /// Constructs a `Self` from its parts. Inputs are treated as bits and shifted into position.
    fn from_parts(sign: bool, exponent: Self::Int, significand: Self::Int) -> Self;

    /// Returns (normalized exponent, normalized significand)
    fn normalize(significand: Self::Int) -> (i32, Self::Int);

    /// Returns `a*b` with `(high_part, low_part)`
    fn wide_multiply(a: Self::Int, b: Self::Int) -> (Self::Int, Self::Int);

    /// Returns `(hi,lo) << count`
    fn wide_left_shift(a: Self::Int, b: Self::Int, count: i32) -> (Self::Int, Self::Int);

    /// Returns `(hi,lo) >> count` (with sticky)
    fn wide_right_shift_with_sticky(a: Self::Int, b: Self::Int, count: u32) -> (Self::Int, Self::Int);

    /// Compute reciprocal (not-corrected)
    fn compute_reciprocal(b: Self::Int) -> Self::Int;
}

// FIXME: Some of this can be removed if RFC Issue #1424 is resolved
//        https://github.com/rust-lang/rfcs/issues/1424
impl Float for f32 {
    type Int = u32;
    fn bits() -> u32 {
        32
    }
    fn significand_bits() -> u32 {
        23
    }
    fn implicit_bit() -> Self::Int {
        1 << Self::significand_bits()
    }
    fn sign_mask() -> Self::Int {
        1 << (Self::bits() - 1)
    }
    fn significand_mask() -> Self::Int {
        (1 << Self::significand_bits()) - 1
    }
    fn exponent_mask() -> Self::Int {
        !(Self::sign_mask() | Self::significand_mask())
    }
    fn repr(self) -> Self::Int {
        unsafe { mem::transmute(self) }
    }
    #[cfg(test)]
    fn eq_repr(self, rhs: Self) -> bool {
        if self.is_nan() && rhs.is_nan() {
            true
        } else {
            self.repr() == rhs.repr()
        }
    }
    fn from_repr(a: Self::Int) -> Self {
        unsafe { mem::transmute(a) }
    }
    fn from_parts(sign: bool, exponent: Self::Int, significand: Self::Int) -> Self {
        Self::from_repr(((sign as Self::Int) << (Self::bits() - 1)) |
            ((exponent << Self::significand_bits()) & Self::exponent_mask()) |
            (significand & Self::significand_mask()))
    }
    fn normalize(significand: Self::Int) -> (i32, Self::Int) {
        let shift = significand.leading_zeros()
            .wrapping_sub((1u32 << Self::significand_bits()).leading_zeros());
        (1i32.wrapping_sub(shift as i32), significand << shift as Self::Int)
    }

    fn wide_multiply(a: Self::Int, b: Self::Int) -> (Self::Int, Self::Int) {
        let product = (a as u64) * (b as u64);
        let hi : u32 = (product >> 32) as u32;
        let lo : u32 = (product & 0xFFFFFFFF) as u32;

        (hi, lo)
    }

    fn wide_left_shift(a: Self::Int, b: Self::Int, count: i32) -> (Self::Int, Self::Int) {
        let hi = a << count;
        let lohi = b >> (32 - count);
        let lolo = b << count;

        (hi | lohi, lolo)
    }

    fn wide_right_shift_with_sticky(a: Self::Int, b: Self::Int, count: u32) -> (Self::Int, Self::Int) {
        if count < 32 {
            let sticky = b << (32 - count);
            let hi = a >> count;
            let lohi = a << (32 - count);
            let lolo = b >> count;

            (hi, lohi|lolo|sticky)
        }
        else if count < 64 {
            let sticky = a << (64 - count) | b;
            let lohi = a >> (count - 32);

            (0, lohi | sticky)
        }
        else {
            (0, a | b)
        }
    }

    fn compute_reciprocal(b: Self::Int) -> Self::Int {
        let q31b = b << 8;
        let reciprocal = 0x7504F333u32 - q31b;
        let correction = Self::wide_multiply(reciprocal, q31b).0.wrapping_neg();
        let reciprocal = {
            let (hi, lo) = Self::wide_multiply(reciprocal, correction);
            (lo >> 31) | (hi << 1)
        };
        let correction = Self::wide_multiply(reciprocal, q31b).0.wrapping_neg();
        let reciprocal = {
            let (hi, lo) = Self::wide_multiply(reciprocal, correction);
            (lo >> 31) | (hi << 1)
        };
        let correction = Self::wide_multiply(reciprocal, q31b).0.wrapping_neg();
        let reciprocal = {
            let (hi, lo) = Self::wide_multiply(reciprocal, correction);
            (lo >> 31) | (hi << 1)
        };
        reciprocal
    }
}
impl Float for f64 {
    type Int = u64;
    fn bits() -> u32 {
        64
    }
    fn significand_bits() -> u32 {
        52
    }
    // Returns the implicit bit of the float format
    fn implicit_bit() -> Self::Int {
        1 << Self::significand_bits()
    }
    fn sign_mask() -> Self::Int {
        1 << (Self::bits() - 1)
    }
    fn significand_mask() -> Self::Int {
        (1 << Self::significand_bits()) - 1
    }
    fn exponent_mask() -> Self::Int {
        !(Self::sign_mask() | Self::significand_mask())
    }
    fn repr(self) -> Self::Int {
        unsafe { mem::transmute(self) }
    }
    #[cfg(test)]
    fn eq_repr(self, rhs: Self) -> bool {
        if self.is_nan() && rhs.is_nan() {
            true
        } else {
            self.repr() == rhs.repr()
        }
    }
    fn from_repr(a: Self::Int) -> Self {
        unsafe { mem::transmute(a) }
    }
    fn from_parts(sign: bool, exponent: Self::Int, significand: Self::Int) -> Self {
        Self::from_repr(((sign as Self::Int) << (Self::bits() - 1)) |
            ((exponent << Self::significand_bits()) & Self::exponent_mask()) |
            (significand & Self::significand_mask()))
    }
    fn normalize(significand: Self::Int) -> (i32, Self::Int) {
        let shift = significand.leading_zeros()
            .wrapping_sub((1u64 << Self::significand_bits()).leading_zeros());
        (1i32.wrapping_sub(shift as i32), significand << shift as Self::Int)
    }

    fn wide_multiply(a: Self::Int, b: Self::Int) -> (Self::Int, Self::Int) {
        let (ahi, alo) = ((a >> 32), (a & 0xFFFFFFFF));
        let (bhi, blo) = ((b >> 32), (b & 0xFFFFFFFF));

        let plolo = alo * blo;
        let plohi = alo * bhi;
        let philo = ahi * blo;
        let phihi = ahi * bhi;

        let r0 = plolo & 0xFFFFFFFF;
        let r1 = (plolo >> 32) + (plohi & 0xFFFFFFFF) + (philo & 0xFFFFFFFF);

        let lo = r0 + (r1 << 32);
        let hi = phihi + (plohi >> 32) + (philo >> 32) + (r1 >> 32);

        (hi, lo)
    }

    fn wide_left_shift(a: Self::Int, b: Self::Int, count: i32) -> (Self::Int, Self::Int) {
        let hi = a << count;
        let lohi = b >> (64 - count);
        let lolo = b << count;

        (hi | lohi, lolo)
    }

    fn wide_right_shift_with_sticky(a: Self::Int, b: Self::Int, count: u32) -> (Self::Int, Self::Int) {
        if count < 64 {
            let sticky = b << (64 - count);
            let hi = a >> count;
            let lohi = a << (64 - count);
            let lolo = b >> count;

            (hi, lohi|lolo|sticky)
        }
        else if count < 128 {
            let sticky = a << (128 - count) | b;
            let lohi = a >> (count - 64);

            (0, lohi | sticky)
        }
        else {
            (0, a | b)
        }
    }

    fn compute_reciprocal(b: Self::Int) -> Self::Int {
        let q31b    : u32 = (b >> 21) as u32;
        let recip32       = 0x7504F333u32 - q31b;
        
        let correction32 = <f32 as Float>::wide_multiply(recip32, q31b).0.wrapping_neg();
        let recip32 = {
            let (hi, lo) = <f32 as Float>::wide_multiply(recip32, correction32);
            (lo >> 31) | (hi << 1)
        };
        let correction32 = <f32 as Float>::wide_multiply(recip32, q31b).0.wrapping_neg();
        let recip32 = {
            let (hi, lo) = <f32 as Float>::wide_multiply(recip32, correction32);
            (lo >> 31) | (hi << 1)
        };
        let correction32 = <f32 as Float>::wide_multiply(recip32, q31b).0.wrapping_neg();
        let recip32 = {
            let (hi, lo) = <f32 as Float>::wide_multiply(recip32, correction32);
            ((lo >> 31) | (hi << 1) - 1)
        };

        let q64blo = b << 11;
        let correction = {
            let r32_64 = recip32 as u64;
            let r1 = r32_64 * q31b as u64;
            let r2 = (r32_64 * q64blo) >> 32;
            
            (r1 + r2).wrapping_neg()
        };
        let (chi, clo) = ((correction >> 32) as u32, (correction & 0xFFFFFFFFu64) as u32);

        let (r1hi, r1lo) = <f32 as Float>::wide_multiply(recip32, chi);
        let (r2hi, _) = <f32 as Float>::wide_multiply(recip32, clo);

        let r1hi_64 = r1hi as u64;
        let r1lo_64 = r1lo as u64;
        let r2hi_64 = r2hi as u64;

        ((r1hi_64 << 32 | r1lo_64) + r2hi_64)
    }
}
