use core::num::Wrapping;

use float::Float;

macro_rules! div {
    ($abi:tt, $intrinsic:ident: $ty:ty) => {
        /// Returns `a * b`
        #[allow(unused_parens)]
        #[cfg_attr(not(test), no_mangle)]
        pub extern $abi fn $intrinsic(a: $ty, b: $ty) -> $ty {
            let one = Wrapping(1 as <$ty as Float>::Int);
            let zero = Wrapping(0 as <$ty as Float>::Int);

            let bits =             Wrapping(<$ty>::bits() as <$ty as Float>::Int);
            let significand_bits = Wrapping(<$ty>::significand_bits() as <$ty as Float>::Int);
            let exponent_bits =    bits - significand_bits - one;
            let max_exponent =     (one << exponent_bits.0 as usize) - one;

            let implicit_bit =     one << significand_bits.0 as usize;
            let significand_mask = implicit_bit - one;
            let sign_bit =         one << (significand_bits + exponent_bits).0 as usize;
            let abs_mask =         sign_bit - one;
            let exponent_mask =    abs_mask ^ significand_mask;
            let inf_rep =          exponent_mask;
            let quiet_bit =        implicit_bit >> 1;
            let qnan_rep =         exponent_mask | quiet_bit;
            let exponent_bias =    max_exponent >> 1;

            let shift_div = if <$ty as Float>::bits() == 32 { 1 } else { 2 };

            let a_rep = Wrapping(a.repr());
            let b_rep = Wrapping(b.repr());
            let a_exponent = Wrapping((a_rep >> significand_bits.0 as usize & max_exponent).0 as <$ty as Float>::Int);
            let b_exponent = Wrapping((b_rep >> significand_bits.0 as usize & max_exponent).0 as <$ty as Float>::Int);
            let mut a_significand = a_rep & significand_mask;
            let mut b_significand = b_rep & significand_mask;
            let a_abs = a_rep & abs_mask;
            let b_abs = b_rep & abs_mask;
            let quotient_sign = (a_rep ^ b_rep) & sign_bit;
            let mut scale : i32 = 0;

            if a_exponent - one >= max_exponent - one ||
                b_exponent - one >= max_exponent - one {
                
                if a_abs > inf_rep {
                    return (<$ty as Float>::from_repr((a_abs | quiet_bit).0));
                }
                if b_abs > inf_rep {
                    return (<$ty as Float>::from_repr((b_abs | quiet_bit).0));
                }

                if a_abs == inf_rep {
                    if b_abs != zero {
                        return (<$ty as Float>::from_repr((a_abs | quotient_sign).0));
                    }
                    return (<$ty as Float>::from_repr(qnan_rep.0));
                }

                if b_abs == inf_rep {
                    if a_abs != zero {
                        return (<$ty as Float>::from_repr((b_abs | quotient_sign).0));
                    }
                    return (<$ty as Float>::from_repr(qnan_rep.0));
                }

                if a_abs == zero {
                    return (<$ty as Float>::from_repr(quotient_sign.0));
                }

                if b_abs == zero {
                    return (<$ty as Float>::from_repr(quotient_sign.0));
                }

                if a_abs < implicit_bit {
                    let (sc, ns) = (<$ty as Float>::normalize(a_significand.0));

                    scale += sc;
                    a_significand.0 = ns;
                }

                if b_abs < implicit_bit {
                    let (sc, ns) = (<$ty as Float>::normalize(b_significand.0));

                    scale -= sc;
                    b_significand.0 = ns;
                }
            } 

            a_significand |= implicit_bit;
            b_significand |= implicit_bit;

            let mut quotient_exponent = a_exponent - b_exponent + Wrapping(scale as <$ty as Float>::Int);
            let reciprocal = (<$ty as Float>::compute_reciprocal(b_significand.0)) - 2;

            a_significand <<= shift_div;

            let mut quotient = Wrapping(<$ty as Float>::wide_multiply(a_significand.0, reciprocal).0);
            let residual = if quotient < (implicit_bit << 1) {
                quotient_exponent -= one;
                (a_significand << (significand_bits + one).0 as usize) - quotient * b_significand
            } else {
                quotient >>= 1;
                (a_significand << significand_bits.0 as usize) - quotient * b_significand
            };

            let written_exponent = quotient_exponent + exponent_bias;

            let final_rep = if written_exponent >= max_exponent {
                inf_rep | quotient_sign
            } else if written_exponent < one {
                quotient_sign
            } else {
                let round = if (residual << 1) > b_significand { 1 } else { 0 };
                let mut abs_result = quotient & significand_mask;
                abs_result |= written_exponent << significand_bits.0 as usize;
                abs_result += Wrapping(round);

                abs_result
            };

            (<$ty as Float>::from_repr(final_rep.0))
        }
    }
}

div!("C", __divsf3: f32);
div!("C", __divdf3: f64);

