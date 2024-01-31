use crate::plonk_tools::{arithmetic, polynomial::Polynomial};
use halo2::halo2curves::{
    bn256::{Fr, G1, G2},
    ff::PrimeField,
};

#[derive(Debug, Clone)]
pub struct PreprocessedInputs {
    pub(crate) n: usize,
    pub(crate) trusted_setup: (Vec<G1>, G2),
    pub(crate) q_m: Vec<Fr>,
    pub(crate) q_l: Vec<Fr>,
    pub(crate) q_r: Vec<Fr>,
    pub(crate) q_o: Vec<Fr>,
    pub(crate) q_c: Vec<Fr>,
    pub(crate) sigma_star_1: Vec<Fr>,
    pub(crate) sigma_star_2: Vec<Fr>,
    pub(crate) sigma_star_3: Vec<Fr>,
    pub(crate) q_m_poly: Polynomial<Fr>,
    pub(crate) q_l_poly: Polynomial<Fr>,
    pub(crate) q_r_poly: Polynomial<Fr>,
    pub(crate) q_o_poly: Polynomial<Fr>,
    pub(crate) q_c_poly: Polynomial<Fr>,
    pub(crate) sigma_1_poly: Polynomial<Fr>,
    pub(crate) sigma_2_poly: Polynomial<Fr>,
    pub(crate) sigma_3_poly: Polynomial<Fr>,
    pub(crate) rou_dom: Vec<Fr>,
    pub(crate) k1: Fr,
    pub(crate) k2: Fr,
    pub(crate) k1_rou: Vec<Fr>,
    pub(crate) k2_rou: Vec<Fr>,
}

pub fn preprocessed_inputs() -> PreprocessedInputs {
    // Our equation is x^3 + x + 5 = 73 and x = 4. We don't want to reveal x value.
    // Selector vectors.
    let q_m = vec![Fr::from(1), Fr::from(1), Fr::from(0), Fr::from(0)];
    let q_l = vec![Fr::from(0), Fr::from(0), Fr::from(1), Fr::from(1)];
    let q_r = vec![Fr::from(0), Fr::from(0), Fr::from(1), Fr::from(0)];
    let q_o = vec![
        Fr::from(1).neg(),
        Fr::from(1).neg(),
        Fr::from(1).neg(),
        Fr::from(1).neg(),
    ];

    let q_c = vec![Fr::from(0), Fr::from(0), Fr::from(0), Fr::from(5)];

    // Number of gates:
    let n = q_m.len();

    // Trusted value's secret "x" multiply by G1 and G2.
    let trusted_setup = arithmetic::trusted_setup(n);

    // Computes a domain of root of unity.
    let rou_dom = arithmetic::rou_dom(n);

    // Shifts roots of unity values by quadratic non-residue k1 and k2.
    let k1 = Fr::from_u128(3);
    let k2 = Fr::from_u128(5);
    let k1_rou = arithmetic::shifted_rou(&rou_dom, k1);
    let k2_rou = arithmetic::shifted_rou(&rou_dom, k2);

    // Computes σ∗ permutation vectors.
    let sigma_star_1 = vec![k1_rou[0], k2_rou[0], k2_rou[1], k2_rou[2]];
    let sigma_star_2 = vec![k1_rou[1], k1_rou[2], rou_dom[0], k1_rou[3]];
    let sigma_star_3 = vec![rou_dom[1], rou_dom[2], rou_dom[3], k2_rou[3]];

    // Computes polynomials.
    let q_l_poly = Polynomial::lagrange(&rou_dom, &q_l);
    let q_r_poly = Polynomial::lagrange(&rou_dom, &q_r);
    let q_o_poly = Polynomial::lagrange(&rou_dom, &q_o);
    let q_m_poly = Polynomial::lagrange(&rou_dom, &q_m);
    let q_c_poly = Polynomial::lagrange(&rou_dom, &q_c);

    // Encodes σ∗ by the three permutation polynomials.
    let sigma_1_poly = Polynomial::lagrange(&rou_dom, &sigma_star_1);
    let sigma_2_poly = Polynomial::lagrange(&rou_dom, &sigma_star_2);
    let sigma_3_poly = Polynomial::lagrange(&rou_dom, &sigma_star_3);

    // Computes common preprocessed input.
    let common_inputs = PreprocessedInputs {
        n,
        trusted_setup,
        q_m,
        q_l,
        q_r,
        q_o,
        q_c,
        sigma_star_1,
        sigma_star_2,
        sigma_star_3,
        q_m_poly,
        q_l_poly,
        q_r_poly,
        q_o_poly,
        q_c_poly,
        sigma_1_poly,
        sigma_2_poly,
        sigma_3_poly,
        rou_dom,
        k1,
        k2,
        k1_rou,
        k2_rou,
    };
    common_inputs
}
