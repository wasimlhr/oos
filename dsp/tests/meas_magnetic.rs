// Copyright (C) 2026 Julian Decker
// SPDX-License-Identifier: GPL-2.0-or-later

//! Inductance, series resistance and reversed current of a simulated inductor.

use oos_dsp::sim::sample;


fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}


fn inductor(l: f64, r: f64, dt: f64, n: usize) -> (Vec<f32>, Vec<f32>) {
    let f = 1000.0;
    let v = sample(
        |t| {
            let x = (t * f).fract();
            let i = if x < 0.5 { -1.0 + 4.0 * x } else { 3.0 - 4.0 * x };
            let di = if x < 0.5 { 4.0 * f } else { -4.0 * f };
            l * di + r * i
        },
        0.0,
        dt,
        n,
    );
    let i = sample(|t| { let x = (t * f).fract(); if x < 0.5 { -1.0 + 4.0 * x } else { 3.0 - 4.0 * x } }, 0.0, dt, n);
    (v, i)
}


#[test]
fn magnetic_inductance_of_a_linear_core_is_the_slope() {
    let l = 1e-3;
    let dt = 2e-6;
    let (v, i) = inductor(l, 0.0, dt, 6000);
    let m = oos_dsp::magnetic::measure(&v, &i, dt, oos_dsp::magnetic::Geom::default());
    assert!(m.kept >= 8, "{:?}", (m.kept, m.dropped, m.note, m.cycles));
    assert!(close(m.l.unwrap(), l, 0.01 * l), "{} {}", m.l.unwrap(), m.note);
    // a sampled pure inductor leaves a fraction of a watt of trapezoid error
    assert!(m.loss.unwrap().abs() < 0.02, "{}", m.loss.unwrap());
    let loop0 = &m.loops[0];
    let (i0, p0) = (loop0[0].1, loop0[0].0);
    for &(psi, cur) in loop0 {
        let est = p0 + l * (cur - i0);
        assert!((psi - est).abs() < 2e-4, "psi {} est {} i {}", psi, est, cur);
    }
}


#[test]
fn magnetic_series_resistance_is_the_loss_and_rcu_removes_it() {
    let dt = 2e-6;
    let (v, i) = inductor(1e-3, 1.0, dt, 6000);
    let m = oos_dsp::magnetic::measure(&v, &i, dt, oos_dsp::magnetic::Geom::default());
    assert!(close(m.loss.unwrap(), 1.0 / 3.0, 0.02), "{} {}", m.loss.unwrap(), m.note);
    let g = oos_dsp::magnetic::Geom { rcu: 1.0, ..oos_dsp::magnetic::Geom::default() };
    let m = oos_dsp::magnetic::measure(&v, &i, dt, g);
    assert!(m.loss.unwrap().abs() < 0.02, "left {} copper {}", m.loss.unwrap(), m.copper);
    assert!(close(m.copper, 1.0 / 3.0, 0.02), "{}", m.copper);
}


#[test]
fn magnetic_short_record_nan_and_reversed_current() {
    let dt = 2e-6;
    let (v, i) = inductor(1e-3, 1.0, dt, 200);
    let m = oos_dsp::magnetic::measure(&v, &i, dt, Default::default());
    assert!(m.l.is_none() && m.loss.is_none() && m.partial, "{}", m.note);
    let nan = vec![f32::NAN; 64];
    let m = oos_dsp::magnetic::measure(&nan, &nan, dt, Default::default());
    assert!(m.l.is_none() && m.loss.is_none(), "{}", m.note);
    let (v, i) = inductor(1e-3, 1.0, dt, 6000);
    let i: Vec<f32> = i.iter().map(|x| -x).collect();
    let m = oos_dsp::magnetic::measure(&v, &i, dt, Default::default());
    assert!(m.loss.unwrap() < 0.0, "{}", m.loss.unwrap());
    assert!(close(m.loss.unwrap(), -1.0 / 3.0, 0.02), "{}", m.loss.unwrap());
}
