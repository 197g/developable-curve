
#[bench]
fn b_mathru_young_laplace(b: &mut Bencher) {
    b.iter(|| {
        let phi_0 = -0.7653981633974483;
        let x_contact = 0.34641208585537364;
        let z_contact = 0.4145142430976436;

        let h_0 = 0.09026305092159702;
        let fac = 0.9;
        let fac_min = 0.2; // 0.01;
        let fac_max = 10.0; // 2.0;
        let n_max = 1e9 as u32;
        let abs_tol = 1e-6;
        let rel_tol = 1e-3;

        let solver = ProportionalControl::new(n_max, h_0, fac, fac_min, fac_max, abs_tol, rel_tol);
        let ode = YoungLaplace {
            phi_0,
            x_contact,
            z_contact,
            startsign: (x_contact * z_contact - phi_0.sin()).signum(),
        };
        let problem = ExplicitInitialValueProblemBuilder::new(&ode, ode.time_span().0, ode.init_cond())
        .t_end(ode.time_span().1)
        .build();
        solver.solve(&problem, &DormandPrince54::default()).unwrap()
    });
}

#[bench]
fn b_young_laplace(b: &mut Bencher) {
    b.iter(|| {
        let phi_0 = -0.7653981633974483;
        let x_contact = 0.34641208585537364;
        let z_contact = 0.4145142430976436;
        let ode = YoungLaplaceOde {
            phi_0,
            startsign: (x_contact * z_contact - phi_0.sin()).signum(),
        };
        solve_ivp(
            &ode,
            (-phi_0.abs(), 0.),
            Coord([x_contact, z_contact]),
            |_, _| true,
            1e-6,
            1e-3,
        )
    });
}

#[bench]
fn b_mathru_harmonic(b: &mut Bencher) {
    b.iter(|| {
        let h_0 = 0.09026305092159702;
        let fac = 0.9;
        let fac_min = 0.2; // 0.01;
        let fac_max = 10.0; // 2.0;
        let n_max = 1e9 as u32;
        let abs_tol = 1e-6;
        let rel_tol = 1e-3;

        let solver = ProportionalControl::new(n_max, h_0, fac, fac_min, fac_max, abs_tol, rel_tol);
        let ode = HarmonicOde {};
        let problem = ExplicitInitialValueProblemBuilder::new(&ode, 0.0, vector![1.; 0.])
            .t_end(10.0)
            .build();
        solver.solve(&problem, &DormandPrince54::default()).unwrap()
    });
}

#[bench]
fn b_harmonic(b: &mut Bencher) {
    b.iter(|| {
        let ode = HarmonicOde {};
        solve_ivp(&ode, (0., 10.), Coord([1., 0.]), |_, _| true, 1e-6, 1e-3)
    });
}

#[bench]
fn b_coupled_array(b: &mut Bencher) {
    b.iter(|| {
        let ode = CoupledHarmonicArOde {};
        solve_ivp(
            &ode,
            (0., 100.),
            Coord([1., 0., 0., 0.5]),
            |_, _| true,
            1e-6,
            1e-3,
        )
    });
}
#[bench]
fn b_coupled_struct(b: &mut Bencher) {
    b.iter(|| {
        let ode = CoupledHarmonicStructOde {};
        solve_ivp(
            &ode,
            (0., 100.),
            Coord([1., 0., 0., 0.5]),
            |_, _| true,
            1e-6,
            1e-3,
        )
    });
}
