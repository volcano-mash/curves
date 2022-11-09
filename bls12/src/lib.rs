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
use ark_std::{cfg_chunks_mut, marker::PhantomData, vec, vec::Vec};
use derivative::Derivative;
use itertools::Itertools;
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
        let a_vec = a
            .into_iter()
            .map(|elem| {
                let elem: Self::G1Prepared = elem.into();
                let mut serialized = vec![0u8; elem.serialized_size(Compress::Yes)];
                elem.serialize_with_mode(&mut serialized, Compress::Yes)
                    .unwrap();
                serialized
            })
            .collect();
        let b_vec = b
            .into_iter()
            .map(|elem| {
                let elem: Self::G2Prepared = elem.into();
                let mut serialized = vec![0u8; elem.serialized_size(Compress::Yes)];
                elem.serialize_with_mode(&mut serialized, Compress::Yes)
                    .unwrap();
                serialized
            })
            .collect();

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
