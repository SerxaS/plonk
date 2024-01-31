use crate::{
    plonk_tools::{
        arithmetic::{self, com_poly, eval, pow},
        polynomial::Polynomial,
    },
    preprocessed_inputs::PreprocessedInputs,
    prover::Proof,
};
use halo2::{
    arithmetic::Field,
    halo2curves::group::Curve,
    halo2curves::{
        bn256::{pairing, Fq, Fr, G1, G2},
        ff::PrimeField,
    },
};

pub fn verifier(common_inputs: PreprocessedInputs, pi_snark: Proof) -> bool {
    // Verifier preprocessed inputs.
    let q_m_com = com_poly(&common_inputs.q_m_poly, &common_inputs.trusted_setup.0);
    let q_l_com = com_poly(&common_inputs.q_l_poly, &common_inputs.trusted_setup.0);
    let q_r_com = com_poly(&common_inputs.q_r_poly, &common_inputs.trusted_setup.0);
    let q_o_com = com_poly(&common_inputs.q_o_poly, &common_inputs.trusted_setup.0);
    let q_c_com = com_poly(&common_inputs.q_c_poly, &common_inputs.trusted_setup.0);
    let sigma_1_com = com_poly(&common_inputs.sigma_1_poly, &common_inputs.trusted_setup.0);
    let sigma_2_com = com_poly(&common_inputs.sigma_2_poly, &common_inputs.trusted_setup.0);
    let sigma_3_com = com_poly(&common_inputs.sigma_3_poly, &common_inputs.trusted_setup.0);

    // Round 1:
    // Validates that the elliptic curve points are actually on the curve E: y^2 = x^3 + 3.
    assert!(
        pi_snark.a_poly_com.to_affine().y.square()
            == pi_snark
                .a_poly_com
                .to_affine()
                .x
                .cube()
                .add(&Fq::from_u128(3)),
        "Points are not on the curve!"
    );
    assert!(
        pi_snark.b_poly_com.to_affine().y.square()
            == pi_snark
                .b_poly_com
                .to_affine()
                .x
                .cube()
                .add(&Fq::from_u128(3)),
        "Points are not on the curve!"
    );
    assert!(
        pi_snark.c_poly_com.to_affine().y.square()
            == pi_snark
                .c_poly_com
                .to_affine()
                .x
                .cube()
                .add(&Fq::from_u128(3)),
        "Points are not on the curve!"
    );
    assert!(
        pi_snark.perm_poly_z_com.to_affine().y.square()
            == pi_snark
                .perm_poly_z_com
                .to_affine()
                .x
                .cube()
                .add(&Fq::from_u128(3)),
        "Points are not on the curve!"
    );
    assert!(
        pi_snark.t_x_lo_com.to_affine().y.square()
            == pi_snark
                .t_x_lo_com
                .to_affine()
                .x
                .cube()
                .add(&Fq::from_u128(3)),
        "Points are not on the curve!"
    );
    assert!(
        pi_snark.t_x_mid_com.to_affine().y.square()
            == pi_snark
                .t_x_mid_com
                .to_affine()
                .x
                .cube()
                .add(&Fq::from_u128(3)),
        "Points are not on the curve!"
    );
    assert!(
        pi_snark.t_x_hi_com.to_affine().y.square()
            == pi_snark
                .t_x_hi_com
                .to_affine()
                .x
                .cube()
                .add(&Fq::from_u128(3)),
        "Points are not on the curve!"
    );
    assert!(
        pi_snark.w_z_x_com.to_affine().y.square()
            == pi_snark
                .w_z_x_com
                .to_affine()
                .x
                .cube()
                .add(&Fq::from_u128(3)),
        "Points are not on the curve!"
    );
    assert!(
        pi_snark.w_z_w_x_com.to_affine().y.square()
            == pi_snark
                .w_z_w_x_com
                .to_affine()
                .x
                .cube()
                .add(&Fq::from_u128(3)),
        "Points are not on the curve!"
    );

    // Round 2:
    // Validates that the field elements in the proof are in Fr.
    let fr_bound = Fr::neg(&Fr::one());
    assert!(
        pi_snark.a_poly_op.coeff[0] <= fr_bound,
        "Points are not on the scalar field (Fr)!"
    );
    assert!(
        pi_snark.b_poly_op.coeff[0] <= fr_bound,
        "Points are not on the scalar field (Fr)!"
    );
    assert!(
        pi_snark.c_poly_op.coeff[0] <= fr_bound,
        "Points are not on the scalar field (Fr)!"
    );
    assert!(
        pi_snark.sigma_1_op.coeff[0] <= fr_bound,
        "Points are not on the scalar field (Fr)!"
    );
    assert!(
        pi_snark.sigma_2_op.coeff[0] <= fr_bound,
        "Points are not on the scalar field (Fr)!"
    );
    assert!(
        pi_snark.perm_poly_z_w_op.coeff[0] <= fr_bound,
        "Points are not on the scalar field (Fr)!"
    );

    // Round 3:
    // Validates (wi)i∈[ℓ] ∈ Fℓ.
    assert!(
        common_inputs.q_m_poly.coeff[0] <= fr_bound,
        "Points are not on the scalar field (Fr)!"
    );
    assert!(
        common_inputs.q_l_poly.coeff[0] <= fr_bound,
        "Points are not on the scalar field (Fr)!"
    );
    assert!(
        common_inputs.q_r_poly.coeff[0] <= fr_bound,
        "Points are not on the scalar field (Fr)!"
    );
    assert!(
        common_inputs.q_o_poly.coeff[0] <= fr_bound,
        "Points are not on the scalar field (Fr)!"
    );
    assert!(
        common_inputs.q_c_poly.coeff[0] <= fr_bound,
        "Points are not on the scalar field (Fr)!"
    );
    assert!(
        common_inputs.sigma_1_poly.coeff[0] <= fr_bound,
        "Points are not on the scalar field (Fr)!"
    );
    assert!(
        common_inputs.sigma_2_poly.coeff[0] <= fr_bound,
        "Points are not on the scalar field (Fr)!"
    );
    assert!(
        common_inputs.sigma_3_poly.coeff[0] <= fr_bound,
        "Points are not on the scalar field (Fr)!"
    );

    // Round 4:
    // Computes challenges β, γ, α, z, v, u ∈ F, from the common inputs, public input,
    // and elements of πSNARK.

    // Changes Fq element to Fr element for transcript hash.
    let a_com_fr_poly =
        Polynomial::new(vec![
            Fr::from_bytes(&Fq::to_bytes(&pi_snark.a_poly_com.x)).unwrap()
        ]);
    let b_com_fr_poly =
        Polynomial::new(vec![
            Fr::from_bytes(&Fq::to_bytes(&pi_snark.b_poly_com.x)).unwrap()
        ]);
    let c_com_fr_poly =
        Polynomial::new(vec![
            Fr::from_bytes(&Fq::to_bytes(&pi_snark.c_poly_com.x)).unwrap()
        ]);

    // Computes permutation challenges (β, γ) ∈ F.
    let beta_inputs = vec![
        common_inputs.q_l_poly.clone(),
        common_inputs.q_r_poly.clone(),
        common_inputs.q_o_poly.clone(),
        common_inputs.q_m_poly.clone(),
        common_inputs.q_c_poly.clone(),
        common_inputs.sigma_1_poly.clone(),
        common_inputs.sigma_2_poly.clone(),
        common_inputs.sigma_3_poly.clone(),
        a_com_fr_poly.clone(),
        b_com_fr_poly.clone(),
        c_com_fr_poly.clone(),
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
        a_com_fr_poly.clone(),
        b_com_fr_poly.clone(),
        c_com_fr_poly.clone(),
        Polynomial::new(vec![Fr::one()]),
    ];
    let gamma = arithmetic::hash_transcript(gamma_inputs);

    // Changes Fq element to Fr element for transcript hash.
    let perm_poly_z_com_fr_poly = Polynomial::new(vec![Fr::from_bytes(&Fq::to_bytes(
        &pi_snark.perm_poly_z_com.x,
    ))
    .unwrap()]);

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
        a_com_fr_poly.clone(),
        b_com_fr_poly.clone(),
        c_com_fr_poly.clone(),
        perm_poly_z_com_fr_poly.clone(),
    ];
    let alfa = arithmetic::hash_transcript(alfa_inputs);
    let alfa_poly = Polynomial::new(vec![alfa]);

    // Changes Fq element to Fr element for transcript hash.
    let t_x_lo_com_fr_poly =
        Polynomial::new(vec![
            Fr::from_bytes(&Fq::to_bytes(&pi_snark.t_x_lo_com.x)).unwrap()
        ]);
    let t_x_mid_com_fr_poly =
        Polynomial::new(vec![
            Fr::from_bytes(&Fq::to_bytes(&pi_snark.t_x_mid_com.x)).unwrap()
        ]);
    let t_x_hi_com_fr_poly =
        Polynomial::new(vec![
            Fr::from_bytes(&Fq::to_bytes(&pi_snark.t_x_hi_com.x)).unwrap()
        ]);

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
        a_com_fr_poly.clone(),
        b_com_fr_poly.clone(),
        c_com_fr_poly.clone(),
        perm_poly_z_com_fr_poly.clone(),
        t_x_lo_com_fr_poly.clone(),
        t_x_mid_com_fr_poly.clone(),
        t_x_hi_com_fr_poly.clone(),
    ];
    let eval_chal_z = arithmetic::hash_transcript(z_inputs);

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
        a_com_fr_poly.clone(),
        b_com_fr_poly.clone(),
        c_com_fr_poly.clone(),
        perm_poly_z_com_fr_poly.clone(),
        t_x_lo_com_fr_poly.clone(),
        t_x_mid_com_fr_poly.clone(),
        t_x_hi_com_fr_poly.clone(),
        pi_snark.a_poly_op.clone(),
        pi_snark.b_poly_op.clone(),
        pi_snark.c_poly_op.clone(),
        pi_snark.sigma_1_op.clone(),
        pi_snark.sigma_2_op.clone(),
        pi_snark.perm_poly_z_w_op.clone(),
    ];
    let open_chal_v = arithmetic::hash_transcript(v_inputs);

    // Changes Fq element to Fr element for transcript hash.
    let w_z_x_com_fr_poly =
        Polynomial::new(vec![
            Fr::from_bytes(&Fq::to_bytes(&pi_snark.w_z_x_com.x)).unwrap()
        ]);
    let w_z_w_x_com_fr_poly =
        Polynomial::new(vec![
            Fr::from_bytes(&Fq::to_bytes(&pi_snark.w_z_w_x_com.x)).unwrap()
        ]);

    // Computes multipoint evaluation challenge u ∈ F.
    let u_inputs = vec![
        common_inputs.q_l_poly.clone(),
        common_inputs.q_r_poly.clone(),
        common_inputs.q_o_poly.clone(),
        common_inputs.q_m_poly.clone(),
        common_inputs.q_c_poly.clone(),
        common_inputs.sigma_1_poly.clone(),
        common_inputs.sigma_2_poly.clone(),
        common_inputs.sigma_3_poly.clone(),
        a_com_fr_poly.clone(),
        b_com_fr_poly.clone(),
        c_com_fr_poly.clone(),
        perm_poly_z_com_fr_poly.clone(),
        t_x_lo_com_fr_poly.clone(),
        t_x_mid_com_fr_poly.clone(),
        t_x_hi_com_fr_poly.clone(),
        pi_snark.a_poly_op.clone(),
        pi_snark.b_poly_op.clone(),
        pi_snark.c_poly_op.clone(),
        pi_snark.sigma_1_op.clone(),
        pi_snark.sigma_2_op.clone(),
        pi_snark.perm_poly_z_w_op.clone(),
        w_z_x_com_fr_poly.clone(),
        w_z_w_x_com_fr_poly.clone(),
    ];
    let multi_eval_chal_u = arithmetic::hash_transcript(u_inputs);

    // Round 5:
    // Computes zero polynomial evaluation.
    let zero_poly = Polynomial::zero_poly(common_inputs.n);
    let zero_poly_eval_z = Polynomial::new(vec![eval(&zero_poly, eval_chal_z)]);

    // Round 6:
    // Computes Lagrange polynomial evaluation L1(z) = ω(zn−1) / n(z−ω).
    let l1_poly_eval_z = Polynomial::new(vec![
        (eval(&Polynomial::l_1_poly(&common_inputs.rou_dom), eval_chal_z)),
    ]);

    // Round 7:
    // We have no public input polynomial.

    // Round 8:
    // Computes r’s constant term.
    let r_0 = l1_poly_eval_z.mul_fr_poly(&alfa.square().neg()).sub_poly(
        &alfa_poly
            .mul_poly(
                &pi_snark
                    .a_poly_op
                    .add_poly(&pi_snark.sigma_1_op.mul_fr_poly(&beta))
                    .add_fr_poly(&gamma),
            )
            .mul_poly(
                &pi_snark
                    .b_poly_op
                    .add_poly(&pi_snark.sigma_2_op.mul_fr_poly(&beta))
                    .add_fr_poly(&gamma),
            )
            .mul_poly(&pi_snark.c_poly_op.add_fr_poly(&gamma))
            .mul_poly(&pi_snark.perm_poly_z_w_op),
    );

    // Round 9:
    // Computes first part of batched polynomial:
    let d_1 = q_m_com * pi_snark.a_poly_op.coeff[0] * pi_snark.b_poly_op.coeff[0]
        + q_l_com * pi_snark.a_poly_op.coeff[0]
        + q_r_com * pi_snark.b_poly_op.coeff[0]
        + q_o_com * pi_snark.c_poly_op.coeff[0]
        + q_c_com;

    let d_2 = pi_snark.perm_poly_z_com
        * pi_snark
            .a_poly_op
            .add_fr_poly(&(beta * eval_chal_z))
            .add_fr_poly(&gamma)
            .mul_poly(
                &pi_snark
                    .b_poly_op
                    .add_fr_poly(&(beta * eval_chal_z * common_inputs.k1))
                    .add_fr_poly(&gamma),
            )
            .mul_poly(
                &pi_snark
                    .c_poly_op
                    .add_fr_poly(&(beta * eval_chal_z * common_inputs.k2))
                    .add_fr_poly(&gamma),
            )
            .mul_poly(&alfa_poly)
            .add_poly(
                &l1_poly_eval_z
                    .mul_fr_poly(&alfa.square())
                    .add_fr_poly(&multi_eval_chal_u),
            )
            .coeff[0];

    let d_3 = sigma_3_com
        * pi_snark
            .a_poly_op
            .add_poly(&pi_snark.sigma_1_op.mul_fr_poly(&beta))
            .add_fr_poly(&gamma)
            .mul_poly(
                &pi_snark
                    .b_poly_op
                    .add_poly(&pi_snark.sigma_2_op.mul_fr_poly(&beta))
                    .add_fr_poly(&gamma),
            )
            .mul_poly(&pi_snark.perm_poly_z_w_op.mul_fr_poly(&(alfa * beta)))
            .coeff[0];

    // Computes quotient polynomial t(X)'s "n". Degree of polynomial is different then index.
    let t_x_n = common_inputs.q_l.len() + 2;

    let d_4 = (pi_snark.t_x_lo_com
        + pi_snark.t_x_mid_com * pow(eval_chal_z, t_x_n)
        + pi_snark.t_x_hi_com * pow(eval_chal_z, 2 * t_x_n))
        * zero_poly_eval_z.coeff[0];

    let batched_d = d_1 + d_2 - d_3 - d_4;

    // Round 10:
    // Computes full batched polynomial commitment [F]1.
    let batched_f = batched_d
        + pi_snark.a_poly_com * open_chal_v
        + pi_snark.b_poly_com * pow(open_chal_v, 2)
        + pi_snark.c_poly_com * pow(open_chal_v, 3)
        + sigma_1_com * pow(open_chal_v, 4)
        + sigma_2_com * pow(open_chal_v, 5);

    // Round 11:
    // Computes group-encoded batch evaluation [E]1.
    let batched_e = G1::generator()
        * &pi_snark
            .a_poly_op
            .mul_fr_poly(&open_chal_v)
            .add_poly(&(pi_snark.b_poly_op).mul_fr_poly(&pow(open_chal_v, 2)))
            .add_poly(&(pi_snark.c_poly_op).mul_fr_poly(&pow(open_chal_v, 3)))
            .add_poly(&(pi_snark.sigma_1_op).mul_fr_poly(&pow(open_chal_v, 4)))
            .add_poly(&(pi_snark.sigma_2_op).mul_fr_poly(&pow(open_chal_v, 5)))
            .add_poly(&pi_snark.perm_poly_z_w_op.mul_fr_poly(&multi_eval_chal_u))
            .add_fr_poly(&r_0.coeff[0].neg())
            .coeff[0];

    // Round 12:
    // Batch validate all evaluations:
    let pair_left = pairing(
        &(pi_snark.w_z_x_com + pi_snark.w_z_w_x_com * multi_eval_chal_u).to_affine(),
        &common_inputs.trusted_setup.1.to_affine(),
    );

    let pair_right = pairing(
        &(pi_snark.w_z_x_com * eval_chal_z
            + pi_snark.w_z_w_x_com * multi_eval_chal_u * eval_chal_z * common_inputs.rou_dom[1]
            + batched_f
            - batched_e)
            .to_affine(),
        &G2::generator().to_affine(),
    );
    pair_left == pair_right
}
