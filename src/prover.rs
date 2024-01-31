use crate::{
    plonk_tools::{
        arithmetic::{self, com_poly, eval, pow},
        polynomial::Polynomial,
    },
    preprocessed_inputs::PreprocessedInputs,
};
use halo2::{
    arithmetic::Field,
    halo2curves::bn256::{Fq, Fr, G1},
};
use rand::thread_rng;

#[derive(Debug, Clone)]
pub struct Proof {
    pub(crate) a_poly_com: G1,
    pub(crate) b_poly_com: G1,
    pub(crate) c_poly_com: G1,
    pub(crate) perm_poly_z_com: G1,
    pub(crate) t_x_lo_com: G1,
    pub(crate) t_x_mid_com: G1,
    pub(crate) t_x_hi_com: G1,
    pub(crate) w_z_x_com: G1,
    pub(crate) w_z_w_x_com: G1,
    pub(crate) a_poly_op: Polynomial<Fr>,
    pub(crate) b_poly_op: Polynomial<Fr>,
    pub(crate) c_poly_op: Polynomial<Fr>,
    pub(crate) sigma_1_op: Polynomial<Fr>,
    pub(crate) sigma_2_op: Polynomial<Fr>,
    pub(crate) perm_poly_z_w_op: Polynomial<Fr>,
}

pub fn prover(common_inputs: &PreprocessedInputs) -> Proof {
    // Witness(secret) vectors. We think of a, b, c as the left, right and output.
    let a = vec![Fr::from(4), Fr::from(16), Fr::from(64), Fr::from(68)];
    let b = vec![Fr::from(4), Fr::from(4), Fr::from(4), Fr::from(0)];
    let c = vec![Fr::from(16), Fr::from(64), Fr::from(68), Fr::from(73)];

    // Round 1:
    // Generates random blinding scalars.
    let rng = thread_rng();
    let blind_sc: Vec<Fr> = (0..11).map(|_| Fr::random(rng.clone())).collect();

    // Computes wire polynomials a(X), b(X), c(X).
    let zero_poly = Polynomial::zero_poly(common_inputs.n);
    let a_poly = Polynomial::new(vec![blind_sc[1], blind_sc[0]])
        .mul_poly(&zero_poly)
        .add_poly(&Polynomial::lagrange(&common_inputs.rou_dom, &a));
    let b_poly = Polynomial::new(vec![blind_sc[3], blind_sc[2]])
        .mul_poly(&zero_poly)
        .add_poly(&Polynomial::lagrange(&common_inputs.rou_dom, &b));
    let c_poly = Polynomial::new(vec![blind_sc[5], blind_sc[4]])
        .mul_poly(&zero_poly)
        .add_poly(&Polynomial::lagrange(&common_inputs.rou_dom, &c));

    // Computes commitment of a, b and c.
    let a_poly_com = com_poly(&a_poly, &common_inputs.trusted_setup.0);
    let b_poly_com = com_poly(&b_poly, &common_inputs.trusted_setup.0);
    let c_poly_com = com_poly(&c_poly, &common_inputs.trusted_setup.0);

    // Changes Fq element to Fr element for transcript hash.
    let a_poly_com_fr =
        Polynomial::new(vec![Fr::from_bytes(&Fq::to_bytes(&a_poly_com.x)).unwrap()]);
    let b_poly_com_fr =
        Polynomial::new(vec![Fr::from_bytes(&Fq::to_bytes(&b_poly_com.x)).unwrap()]);
    let c_poly_com_fr =
        Polynomial::new(vec![Fr::from_bytes(&Fq::to_bytes(&c_poly_com.x)).unwrap()]);

    // Round 2:
    // Computes permutation challenges (β, γ) ∈ F
    let beta_inputs = vec![
        common_inputs.q_l_poly.clone(),
        common_inputs.q_r_poly.clone(),
        common_inputs.q_o_poly.clone(),
        common_inputs.q_m_poly.clone(),
        common_inputs.q_c_poly.clone(),
        common_inputs.sigma_1_poly.clone(),
        common_inputs.sigma_2_poly.clone(),
        common_inputs.sigma_3_poly.clone(),
        a_poly_com_fr.clone(),
        b_poly_com_fr.clone(),
        c_poly_com_fr.clone(),
        Polynomial::new(vec![Fr::zero()]),
    ];
    let beta = arithmetic::hash_transcript(beta_inputs);

    let gamma_inputs = vec![
        common_inputs.q_l_poly.clone(),
        common_inputs.q_r_poly.clone(),
        common_inputs.q_o_poly.clone(),
        common_inputs.q_m_poly.clone(),
        common_inputs.q_c_poly.clone(),
        common_inputs.sigma_1_poly.clone(),
        common_inputs.sigma_2_poly.clone(),
        common_inputs.sigma_3_poly.clone(),
        a_poly_com_fr.clone(),
        b_poly_com_fr.clone(),
        c_poly_com_fr.clone(),
        Polynomial::new(vec![Fr::ONE]),
    ];
    let gamma = arithmetic::hash_transcript(gamma_inputs);

    // Computes accumulator.
    let acc_poly = Polynomial::perm_poly_acc(
        &a,
        &b,
        &c,
        &common_inputs.sigma_star_1,
        &common_inputs.sigma_star_2,
        &common_inputs.sigma_star_3,
        &common_inputs.rou_dom,
        common_inputs.k1,
        common_inputs.k2,
        beta,
        gamma,
    );

    // Computes permutation polynomial z(X).
    let perm_poly_z = Polynomial::new(vec![blind_sc[8], blind_sc[7], blind_sc[6]])
        .mul_poly(&zero_poly)
        .add_poly(&acc_poly);

    // Computes commitment of permutation polynomial z(X).
    let perm_poly_z_com = com_poly(&perm_poly_z, &common_inputs.trusted_setup.0);

    // Changes Fq element to Fr element for transcript hash.
    let perm_poly_z_com_fr_poly =
        Polynomial::new(vec![
            Fr::from_bytes(&Fq::to_bytes(&perm_poly_z_com.x)).unwrap()
        ]);

    // Round 3:
    // Computes quotient challenge α ∈ F.
    let alfa_inputs = vec![
        common_inputs.q_l_poly.clone(),
        common_inputs.q_r_poly.clone(),
        common_inputs.q_o_poly.clone(),
        common_inputs.q_m_poly.clone(),
        common_inputs.q_c_poly.clone(),
        common_inputs.sigma_1_poly.clone(),
        common_inputs.sigma_2_poly.clone(),
        common_inputs.sigma_3_poly.clone(),
        a_poly_com_fr.clone(),
        b_poly_com_fr.clone(),
        c_poly_com_fr.clone(),
        perm_poly_z_com_fr_poly.clone(),
    ];
    let alfa = arithmetic::hash_transcript(alfa_inputs);

    // Computes quotient polynomial t(X).
    let t_x_1 = a_poly
        .mul_poly(&b_poly)
        .mul_poly(&common_inputs.q_m_poly)
        .add_poly(&a_poly.mul_poly(&common_inputs.q_l_poly))
        .add_poly(&b_poly.mul_poly(&common_inputs.q_r_poly))
        .add_poly(&c_poly.mul_poly(&common_inputs.q_o_poly))
        .add_poly(&common_inputs.q_c_poly);

    let t_x_2 = (a_poly
        .add_poly(&Polynomial::new(vec![gamma, beta]))
        .mul_poly(&b_poly.add_poly(&Polynomial::new(vec![gamma, beta * common_inputs.k1])))
        .mul_poly(&c_poly.add_poly(&Polynomial::new(vec![gamma, beta * common_inputs.k2])))
        .mul_poly(&perm_poly_z))
    .mul_fr_poly(&alfa);

    // Computes shifted permutation polynomial z(Xω).
    let shifted_z_x = Polynomial::shift_perm_poly_z(&perm_poly_z, &common_inputs.rou_dom);
    let t_x_3 = a_poly
        .add_poly(&common_inputs.sigma_1_poly.mul_fr_poly(&beta))
        .add_fr_poly(&gamma)
        .mul_poly(
            &b_poly
                .add_poly(&common_inputs.sigma_2_poly.mul_fr_poly(&beta))
                .add_fr_poly(&gamma),
        )
        .mul_poly(
            &c_poly
                .add_poly(&common_inputs.sigma_3_poly.mul_fr_poly(&beta))
                .add_fr_poly(&gamma),
        )
        .mul_poly(&shifted_z_x)
        .mul_fr_poly(&alfa);

    let t_x_4 = &perm_poly_z
        .sub_fr_poly(&Fr::ONE)
        .mul_poly(&Polynomial::l_1_poly(&common_inputs.rou_dom))
        .mul_fr_poly(&alfa.square());

    let t_x = t_x_1
        .add_poly(&t_x_2)
        .sub_poly(&t_x_3)
        .add_poly(&t_x_4)
        .div_poly(&zero_poly)
        .0;

    // Computes quotient polynomial t(X)'s "n". Degree of polynomial is different then index.
    let t_x_n = t_x.coeff.len() / 3;
    let split_t_x = Polynomial::split_t_x(t_x);
    let blind_polys = Polynomial::blind_poly(t_x_n, blind_sc[9], blind_sc[10]);
    let t_x_lo = split_t_x.0.add_poly(&blind_polys.0);
    let t_x_mid = split_t_x
        .1
        .add_poly(&blind_polys.1.sub_fr_poly(&blind_sc[9]));
    let t_x_hi = split_t_x.2.sub_fr_poly(&blind_sc[10]);

    // Computes commitment of t_x_lo, t_x_mid and t_x_hi.
    let t_x_lo_com = com_poly(&t_x_lo, &common_inputs.trusted_setup.0);
    let t_x_mid_com = com_poly(&t_x_mid, &common_inputs.trusted_setup.0);
    let t_x_hi_com = com_poly(&t_x_hi, &common_inputs.trusted_setup.0);

    // Changes Fq element to Fr element for transcript hash.
    let t_x_lo_com_fr_poly =
        Polynomial::new(vec![Fr::from_bytes(&Fq::to_bytes(&t_x_lo_com.x)).unwrap()]);
    let t_x_mid_com_fr_poly =
        Polynomial::new(vec![Fr::from_bytes(&Fq::to_bytes(&t_x_mid_com.x)).unwrap()]);
    let t_x_hi_com_fr_poly =
        Polynomial::new(vec![Fr::from_bytes(&Fq::to_bytes(&t_x_hi_com.x)).unwrap()]);

    // Round 4:
    // Computes evaluation challenge z ∈ F.
    let z_inputs = vec![
        common_inputs.q_l_poly.clone(),
        common_inputs.q_r_poly.clone(),
        common_inputs.q_o_poly.clone(),
        common_inputs.q_m_poly.clone(),
        common_inputs.q_c_poly.clone(),
        common_inputs.sigma_1_poly.clone(),
        common_inputs.sigma_2_poly.clone(),
        common_inputs.sigma_3_poly.clone(),
        a_poly_com_fr.clone(),
        b_poly_com_fr.clone(),
        c_poly_com_fr.clone(),
        perm_poly_z_com_fr_poly.clone(),
        t_x_lo_com_fr_poly.clone(),
        t_x_mid_com_fr_poly.clone(),
        t_x_hi_com_fr_poly.clone(),
    ];
    let eval_chal_z = arithmetic::hash_transcript(z_inputs);

    // Computes opening evaluations.
    let a_poly_op = Polynomial::new(vec![eval(&a_poly, eval_chal_z)]);
    let b_poly_op = Polynomial::new(vec![eval(&b_poly, eval_chal_z)]);
    let c_poly_op = Polynomial::new(vec![eval(&c_poly, eval_chal_z)]);
    let sigma_1_op = Polynomial::new(vec![eval(&common_inputs.sigma_1_poly, eval_chal_z)]);
    let sigma_2_op = Polynomial::new(vec![eval(&common_inputs.sigma_2_poly, eval_chal_z)]);
    let perm_poly_z_w = common_inputs.rou_dom[1] * eval_chal_z;
    let perm_poly_z_w_op = Polynomial::new(vec![eval(&perm_poly_z, perm_poly_z_w)]);

    // Round 5:
    // Computes opening challenge v ∈ F.
    let v_inputs = vec![
        common_inputs.q_l_poly.clone(),
        common_inputs.q_r_poly.clone(),
        common_inputs.q_o_poly.clone(),
        common_inputs.q_m_poly.clone(),
        common_inputs.q_c_poly.clone(),
        common_inputs.sigma_1_poly.clone(),
        common_inputs.sigma_2_poly.clone(),
        common_inputs.sigma_3_poly.clone(),
        a_poly_com_fr.clone(),
        b_poly_com_fr.clone(),
        c_poly_com_fr.clone(),
        perm_poly_z_com_fr_poly.clone(),
        t_x_lo_com_fr_poly.clone(),
        t_x_mid_com_fr_poly.clone(),
        t_x_hi_com_fr_poly.clone(),
        a_poly_op.clone(),
        b_poly_op.clone(),
        c_poly_op.clone(),
        sigma_1_op.clone(),
        sigma_2_op.clone(),
        perm_poly_z_w_op.clone(),
    ];
    let open_chal_v = arithmetic::hash_transcript(v_inputs);

    // Computes linearisation polynomial r(X).
    let r_x_1 = a_poly_op
        .mul_poly(&b_poly_op.mul_poly(&common_inputs.q_m_poly))
        .add_poly(&a_poly_op.mul_poly(&common_inputs.q_l_poly))
        .add_poly(&b_poly_op.mul_poly(&common_inputs.q_r_poly))
        .add_poly(&c_poly_op.mul_poly(&common_inputs.q_o_poly))
        .add_poly(&common_inputs.q_c_poly);

    let r_x_2 = &a_poly_op
        .add_fr_poly(&(beta * eval_chal_z))
        .add_fr_poly(&gamma)
        .mul_fr_poly(&alfa)
        .mul_poly(
            &b_poly_op
                .add_fr_poly(&(beta * common_inputs.k1 * eval_chal_z))
                .add_fr_poly(&gamma),
        )
        .mul_poly(
            &c_poly_op
                .add_fr_poly(&(beta * common_inputs.k2 * eval_chal_z))
                .add_fr_poly(&gamma),
        )
        .mul_poly(&perm_poly_z);

    let r_x_3 = &a_poly_op
        .add_poly(&sigma_1_op.mul_fr_poly(&beta))
        .add_fr_poly(&gamma)
        .mul_fr_poly(&alfa)
        .mul_poly(
            &b_poly_op
                .add_poly(&sigma_2_op.mul_fr_poly(&beta))
                .add_fr_poly(&gamma),
        )
        .mul_poly(
            &c_poly_op
                .add_poly(&common_inputs.sigma_3_poly.mul_fr_poly(&beta))
                .add_fr_poly(&gamma),
        )
        .mul_poly(&perm_poly_z_w_op);

    let l_1_poly_eval_z = Polynomial::new(vec![
        (eval(&Polynomial::l_1_poly(&common_inputs.rou_dom), eval_chal_z)),
    ]);

    let r_x_4 = perm_poly_z
        .sub_poly(&Polynomial::new(vec![Fr::ONE]))
        .mul_fr_poly(&alfa.square())
        .mul_poly(&l_1_poly_eval_z);

    let zero_poly_eval_z = eval(&zero_poly, eval_chal_z);

    let r_x_5 = &t_x_lo
        .add_poly(&t_x_mid.mul_fr_poly(&pow(eval_chal_z, t_x_n)))
        .add_poly(&t_x_hi.mul_fr_poly(&pow(eval_chal_z, 2 * t_x_n)))
        .mul_fr_poly(&zero_poly_eval_z);

    let r_x = r_x_1
        .add_poly(&r_x_2)
        .sub_poly(&r_x_3)
        .add_poly(&r_x_4)
        .sub_poly(&r_x_5);

    // Computes opening proof polynomial Wz(X).
    let w_z_x = r_x
        .add_poly(&a_poly.sub_poly(&a_poly_op).mul_fr_poly(&open_chal_v))
        .add_poly(&(&b_poly.sub_poly(&b_poly_op)).mul_fr_poly(&pow(open_chal_v, 2)))
        .add_poly(&(&c_poly.sub_poly(&c_poly_op)).mul_fr_poly(&pow(open_chal_v, 3)))
        .add_poly(
            &(&common_inputs.sigma_1_poly.sub_poly(&sigma_1_op)).mul_fr_poly(&pow(open_chal_v, 4)),
        )
        .add_poly(
            &(&common_inputs.sigma_2_poly.sub_poly(&sigma_2_op)).mul_fr_poly(&pow(open_chal_v, 5)),
        )
        .div_poly(&Polynomial::new(vec![eval_chal_z.neg(), Fr::ONE]))
        .0;

    // Computes opening proof polynomial Wzω(X).
    let w_z_w_x = perm_poly_z
        .sub_poly(&perm_poly_z_w_op)
        .div_poly(&Polynomial::new(vec![perm_poly_z_w.neg(), Fr::ONE]))
        .0;

    let w_z_x_com = com_poly(&w_z_x, &common_inputs.trusted_setup.0);
    let w_z_w_x_com = com_poly(&w_z_w_x, &common_inputs.trusted_setup.0);

    // Changes Fq element to Fr element for transcript hash.
    let w_z_x_com_fr_poly =
        Polynomial::new(vec![Fr::from_bytes(&Fq::to_bytes(&w_z_x_com.x)).unwrap()]);
    let w_z_w_x_com_fr_poly =
        Polynomial::new(vec![Fr::from_bytes(&Fq::to_bytes(&w_z_w_x_com.x)).unwrap()]);

    // Computes multipoint evaluation challenge u ∈ F
    let u_inputs = vec![
        common_inputs.q_l_poly.clone(),
        common_inputs.q_r_poly.clone(),
        common_inputs.q_o_poly.clone(),
        common_inputs.q_m_poly.clone(),
        common_inputs.q_c_poly.clone(),
        common_inputs.sigma_1_poly.clone(),
        common_inputs.sigma_2_poly.clone(),
        common_inputs.sigma_3_poly.clone(),
        a_poly_com_fr,
        b_poly_com_fr,
        c_poly_com_fr,
        perm_poly_z_com_fr_poly,
        t_x_lo_com_fr_poly.clone(),
        t_x_mid_com_fr_poly.clone(),
        t_x_hi_com_fr_poly.clone(),
        a_poly_op.clone(),
        b_poly_op.clone(),
        c_poly_op.clone(),
        sigma_1_op.clone(),
        sigma_2_op.clone(),
        perm_poly_z_w_op.clone(),
        w_z_x_com_fr_poly.clone(),
        w_z_w_x_com_fr_poly.clone(),
    ];
    let multi_eval_chal_u = arithmetic::hash_transcript(u_inputs);

    // Return πSNARK:
    let pi_snark = Proof {
        a_poly_com,
        b_poly_com,
        c_poly_com,
        perm_poly_z_com,
        t_x_lo_com,
        t_x_mid_com,
        t_x_hi_com,
        w_z_x_com,
        w_z_w_x_com,
        a_poly_op,
        b_poly_op,
        c_poly_op,
        sigma_1_op,
        sigma_2_op,
        perm_poly_z_w_op,
    };
    pi_snark
}
