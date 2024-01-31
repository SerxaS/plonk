use crate::{plonk_tools::polynomial::Polynomial, poseidon_hash::sponge::PoseidonSponge};
use halo2::{
    arithmetic::Field,
    halo2curves::{
        bn256::{Fr, G1, G2},
        ff::PrimeField,
    },
};
use rand::thread_rng;

/// Calculates exponent of a number as F element.
pub fn pow(base: Fr, exp: usize) -> Fr {
    let mut mul = Fr::ONE;

    for _ in 0..exp {
        mul *= base
    }
    mul
}

/// Evaluates polynomial in the field.
pub fn eval(poly: &Polynomial<Fr>, val: Fr) -> Fr {
    let mut eval = Fr::ZERO;

    for (i, j) in poly.iter().enumerate() {
        eval += *j * pow(val, i.try_into().unwrap())
    }
    eval
}

/// Computes commitment.
pub fn com_poly(poly: &Polynomial<Fr>, x_g1: &Vec<G1>) -> G1 {
    let mut com_poly = G1::generator() * Fr::ZERO;

    for i in 0..poly.len() {
        com_poly += x_g1[i] * poly[i];
    }
    com_poly
}

/// Computes a domain of root of unity elements.
pub fn rou_dom(n: usize) -> Vec<Fr> {
    // Computes required root of unity.
    let mut len = n;
    let mut req_rou = Fr::ROOT_OF_UNITY;
    let mut counter = 0;

    while len / 2 >= 1 {
        len = len / 2;
        counter += 1;
    }

    for _ in 0..(28 - counter) {
        req_rou = req_rou.square();
    }

    // Computes domain.
    let mut rou_dom = Vec::new();

    for i in 0..n {
        rou_dom.push(pow(req_rou, i));
    }
    rou_dom
}

/// Computes a vector of shifted rou elements for compute sigmas and permutation polynomial z(X).
pub fn shifted_rou(rou_dom: &Vec<Fr>, k: Fr) -> Vec<Fr> {
    let mut shifted_rou = vec![Fr::ONE; rou_dom.len()];

    for i in 0..rou_dom.len() {
        shifted_rou[i] = k * rou_dom[i];
    }
    shifted_rou
}

/// Concatenates preprocessed inputs and hashes them.
pub(crate) fn hash_transcript(inputs: Vec<Polynomial<Fr>>) -> Fr {
    let mut transcript = Polynomial::new(Vec::new());

    for i in inputs {
        let mut sponge = PoseidonSponge::new();
        sponge.update(&i.coeff);
        let squeeze = PoseidonSponge::squeeze(&mut sponge);
        transcript.push(squeeze);
    }
    let mut hashed_transcript = Polynomial::new(Vec::new());
    let mut sponge = PoseidonSponge::new();
    sponge.update(&transcript.coeff);
    let squeeze = PoseidonSponge::squeeze(&mut sponge);
    hashed_transcript.coeff.push(squeeze);
    hashed_transcript.coeff[0]
}

/// The secret x is generated from random number but in practice this is usually implemented
/// via a secure multiparty computation (MPC). Two sets will be distributed publicly, one for
/// ([x^0]_1,...,[x^(i + 5)]_1), and one for [x]_2. The secret s is then discarded forever.
pub fn trusted_setup(num_gate: usize) -> (Vec<G1>, G2) {
    let rng = thread_rng();
    let secret_x = Fr::random(rng.clone());
    let mut x_g1 = Vec::new();

    for i in 0..num_gate + 5 {
        let trusted_x_g1 = G1::generator() * pow(secret_x, i.try_into().unwrap());
        x_g1.push(trusted_x_g1);
    }
    let x_g2 = G2::generator() * secret_x;
    (x_g1, x_g2)
}
