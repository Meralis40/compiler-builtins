use core::num::Wrapping;

use float::Float;

macro_rules! mul {
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

            let a_rep = Wrapping(a.repr());
            let b_rep = Wrapping(b.repr());
            let a_exponent = Wrapping((a_rep >> significand_bits.0 as usize & max_exponent).0 as <$ty as Float>::Int);
            let b_exponent = Wrapping((b_rep >> significand_bits.0 as usize & max_exponent).0 as <$ty as Float>::Int);
            let mut a_significand = a_rep & significand_mask;
            let mut b_significand = b_rep & significand_mask;
            let a_abs = a_rep & abs_mask;
            let b_abs = b_rep & abs_mask;
            let product_sign = (a_rep ^ b_rep) & sign_bit;
            let mut scale : i32 = 0;

            // Detect if a or b is zero, denormal, infinity or NaN
            if a_exponent - one >= max_exponent - one ||
                b_exponent - one >= max_exponent - one {
                
                // NaN * anything = qNaN
                if a_abs > inf_rep {
                    return (<$ty as Float>::from_repr((a_abs | quiet_bit).0));
                }
                // anything * NaN = qNaN
                if b_abs > inf_rep {
                    return (<$ty as Float>::from_repr((b_abs | quiet_bit).0));
                }

                if a_abs == inf_rep {
                    if b_abs != zero {
                        return (<$ty as Float>::from_repr((a_abs | product_sign).0));
                    }
                    return (<$ty as Float>::from_repr(qnan_rep.0));
                }

                if b_abs == inf_rep {
                    if a_abs != zero {
                        return (<$ty as Float>::from_repr((b_abs | product_sign).0));
                    }
                    return (<$ty as Float>::from_repr(qnan_rep.0));
                }

                if a_abs == zero {
                    return (<$ty as Float>::from_repr(product_sign.0));
                }

                if b_abs == zero {
                    return (<$ty as Float>::from_repr(product_sign.0));
                }

                if a_abs < implicit_bit {
                    let (sc, ns) = (<$ty as Float>::normalize(a_significand.0));

                    scale += sc;
                    a_significand.0 = ns;
                }

                if b_abs < implicit_bit {
                    let (sc, ns) = (<$ty as Float>::normalize(b_significand.0));

                    scale += sc;
                    b_significand.0 = ns;
                }
            } 

            a_significand |= implicit_bit;
            b_significand |= implicit_bit;

            let prod = <$ty as Float>::wide_multiply(a_significand.0,
                (b_significand << exponent_bits.0 as usize).0);
            let mut product_hi = Wrapping(prod.0);
            let mut product_lo = Wrapping(prod.1);

            let product_exponent = Wrapping(a_exponent.0 + b_exponent.0 + (scale as <$ty as Float>::Int) - exponent_bias.0);

            if product_hi & implicit_bit != zero {
                product_hi += one;
            } else {
                let prod = <$ty as Float>::wide_left_shift(product_hi.0, product_lo.0, 1);
                product_hi.0 = prod.0;
                product_lo.0 = prod.1;
            }

            if product_exponent >= max_exponent {
                return (<$ty as Float>::from_repr((inf_rep | product_sign).0));
            }

            if product_exponent <= zero {
                let shift = one - product_exponent;

                if shift.0 as u32 >= <$ty as Float>::bits() {
                    return (<$ty as Float>::from_repr(product_sign.0));
                }

                let prod = <$ty as Float>::wide_right_shift_with_sticky(product_hi.0, product_lo.0, shift.0 as u32);

                product_hi.0 = prod.0;
                product_lo.0 = prod.1;
            } else {
                product_hi &= significand_mask;
                product_hi |= product_exponent << significand_bits.0 as usize;
            }

            product_hi |= product_sign;

            if product_lo > sign_bit {
                product_hi += one;
            }

            if product_lo == sign_bit {
                product_hi += (product_hi & one);
            }

            <$ty as Float>::from_repr(product_hi.0)
        }
    }
}

mul!("C", __mulsf3: f32);
mul!("C", __muldf3: f64);

