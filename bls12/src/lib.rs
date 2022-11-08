#![cfg_attr(not(feature = "std"), no_std)]

use ark_ec::{
    models::{short_weierstrass::SWCurveConfig, CurveConfig},
    pairing::{MillerLoopOutput, Pairing, PairingOutput},
    AffineRepr,
};
use ark_ff::{
    fields::{
        fp12_2over3over2::{Fp12, Fp12Config},
        fp2::Fp2Config,
        fp6_3over2::Fp6Config,
        Fp2,
    },
    BitIteratorBE, CyclotomicMultSubgroup, Field, PrimeField,
};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress};
use ark_std::{cfg_chunks_mut, marker::PhantomData, vec::Vec};
use derivative::Derivative;
use num_traits::{One, Zero};
use sp_io::crypto::bls12_381_multi_miller_loop;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// A particular BLS12 group can have G2 being either a multiplicative or a
/// divisive twist.
pub enum TwistType {
    M,
    D,
}

pub trait Bls12Parameters: 'static {
    /// Parameterizes the BLS12 family.
    const X: &'static [u64];
    /// Is `Self::X` negative?
    const X_IS_NEGATIVE: bool;
    /// What kind of twist is this?
    const TWIST_TYPE: TwistType;

    type Fp: PrimeField + Into<<Self::Fp as PrimeField>::BigInt>;
    type Fp2Config: Fp2Config<Fp = Self::Fp>;
    type Fp6Config: Fp6Config<Fp2Config = Self::Fp2Config>;
    type Fp12Config: Fp12Config<Fp6Config = Self::Fp6Config>;
    type G1Parameters: SWCurveConfig<BaseField = Self::Fp>;
    type G2Parameters: SWCurveConfig<
        BaseField = Fp2<Self::Fp2Config>,
        ScalarField = <Self::G1Parameters as CurveConfig>::ScalarField,
    >;
}

pub mod g1;
pub mod g2;

pub use self::{
    g1::{G1Affine, G1Prepared, G1Projective},
    g2::{G2Affine, G2Prepared, G2Projective},
};

#[derive(Derivative)]
#[derivative(Copy, Clone, PartialEq, Eq, Debug, Hash)]
pub struct Bls12<P: Bls12Parameters>(PhantomData<fn() -> P>);

impl<P: Bls12Parameters> Bls12<P> {
    // Evaluate the line function at point p.
    fn ell(f: &mut Fp12<P::Fp12Config>, coeffs: &g2::EllCoeff<P>, p: &G1Affine<P>) {
        let mut c0 = coeffs.0;
        let mut c1 = coeffs.1;
        let mut c2 = coeffs.2;
        let (px, py) = p.xy().unwrap();

        match P::TWIST_TYPE {
            TwistType::M => {
                c2.mul_assign_by_fp(py);
                c1.mul_assign_by_fp(px);
                f.mul_by_014(&c0, &c1, &c2);
            },
            TwistType::D => {
                c0.mul_assign_by_fp(py);
                c1.mul_assign_by_fp(px);
                f.mul_by_034(&c0, &c1, &c2);
            },
        }
    }

    // Exponentiates `f` by `Self::X`, and stores the result in `result`.
    fn exp_by_x(f: &Fp12<P::Fp12Config>, result: &mut Fp12<P::Fp12Config>) {
        *result = f.cyclotomic_exp(P::X);
        if P::X_IS_NEGATIVE {
            result.cyclotomic_inverse_in_place();
        }
    }
}

impl<P: Bls12Parameters> Pairing for Bls12<P> {
    type BaseField = <P::G1Parameters as CurveConfig>::BaseField;
    type ScalarField = <P::G1Parameters as CurveConfig>::ScalarField;
    type G1 = G1Projective<P>;
    type G1Affine = G1Affine<P>;
    type G1Prepared = G1Prepared<P>;
    type G2 = G2Projective<P>;
    type G2Affine = G2Affine<P>;
    type G2Prepared = G2Prepared<P>;
    type TargetField = Fp12<P::Fp12Config>;

    fn multi_miller_loop(
        a: impl IntoIterator<Item = impl Into<Self::G1Prepared>>,
        b: impl IntoIterator<Item = impl Into<Self::G2Prepared>>,
    ) -> MillerLoopOutput<Self> {
        use itertools::Itertools;

        let mut a_vec = Vec::new();
        let mut b_vec = Vec::new();

        a.into_iter()
            .map(|left| left.into())
            .zip_eq(b.into_iter().map(|right| right.into()))
            .for_each(|(left, right)| {
                let mut serialized_a = Vec::new();
                let mut serialized_b = Vec::new();
                left.serialize_with_mode(&mut serialized_a, Compress::Yes)
                    .unwrap();
                right
                    .serialize_with_mode(&mut serialized_b, Compress::Yes)
                    .unwrap();
                a_vec.push(serialized_a);
                b_vec.push(serialized_b);
            });

        let res = bls12_381_multi_miller_loop(a_vec, b_vec);
        let f: Self::TargetField =
            Fp12::deserialize_with_mode(&res[..], Compress::Yes, ark_serialize::Validate::No)
                .unwrap();

        MillerLoopOutput(f)
    }

    fn final_exponentiation(f: MillerLoopOutput<Self>) -> Option<PairingOutput<Self>> {
        let mut out: [u8; 576] = [0; 576];
        f.0.serialize_with_mode(&mut out[..], Compress::Yes)
            .unwrap();

        let res = sp_io::crypto::bls12_381_final_exponentiation(&out[..]);

        let r: Self::TargetField =
            Fp12::deserialize_with_mode(&res[..], Compress::Yes, ark_serialize::Validate::No)
                .unwrap();

        Some(PairingOutput(r))
    }
}
