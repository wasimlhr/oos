// Copyright (C) 2026 Julian Decker
// SPDX-License-Identifier: GPL-2.0-or-later

//! Double-pulse timings, peak current and recovery charge of simulated records.

use oos_dsp::sim::sample;


fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

/// Piecewise linear function through (t, value) points.
fn pwl(p: &[(f64, f64)], t: f64) -> f64 {
    if t <= p[0].0 {
        return p[0].1;
    }
    for w in p.windows(2) {
        if t <= w[1].0 {
            return w[0].1 + (w[1].1 - w[0].1) * (t - w[0].0) / (w[1].0 - w[0].0);
        }
    }
    p[p.len() - 1].1
}


#[test]
fn double_pulse_timings_peak_and_recovery_charge() {
    let dt = 1e-9;
    let tend = 13e-6;
    let n = (tend / dt) as usize;
    // gate: two pulses, 0..10 V in 20 ns. Current plateau 10 A. Recovery triangle 4 A, 200 ns base.
    let gate = [(0.0, 0.0), (1e-6, 0.0), (1.02e-6, 10.0), (6e-6, 10.0), (6.02e-6, 0.0), (10e-6, 0.0), (10.02e-6, 10.0), (12e-6, 10.0), (12.02e-6, 0.0)];
    let id = [
        (0.0, 0.0),
        (1.1e-6, 0.0),
        (1.2e-6, 10.0),
        (6.05e-6, 10.0),
        (6.15e-6, 0.0),
        (10.05e-6, 0.0),
        (10.15e-6, 10.0),
        (10.25e-6, 14.0),
        (10.35e-6, 10.0),
        (12e-6, 10.0),
    ];
    let vds = [
        (0.0, 400.0),
        (1.05e-6, 400.0),
        (1.15e-6, 1.0),
        (6.1e-6, 1.0),
        (6.2e-6, 400.0),
        (6.25e-6, 480.0),
        (6.4e-6, 400.0),
        (10.05e-6, 400.0),
        (10.15e-6, 1.0),
        (12e-6, 1.0),
    ];
    let g = sample(|t| pwl(&gate, t), 0.0, dt, n);
    let i = sample(|t| pwl(&id, t), 0.0, dt, n);
    let v = sample(|t| pwl(&vds, t), 0.0, dt, n);
    let p = oos_dsp::dpt::measure(&v, &i, Some(&g), oos_dsp::dpt::Fourth::Off, dt, 10.0);
    assert!(close(p.td_on.unwrap(), 58e-9, 2e-9), "tdon {}", p.td_on.unwrap());
    assert!(close(p.tr.unwrap(), 80e-9, 2e-9), "tr {}", p.tr.unwrap());
    assert!(close(p.ton.unwrap(), 138e-9, 3e-9), "ton {}", p.ton.unwrap());
    assert!(close(p.vpeak.unwrap(), 480.0, 1e-6), "vpeak {}", p.vpeak.unwrap());
    assert!(close(p.ipeak.unwrap(), 14.0, 1e-3), "ipeak {}", p.ipeak.unwrap());
    assert!(close(p.irrm.unwrap(), 4.0, 0.05), "irrm {}", p.irrm.unwrap());
    // symmetric triangle, Trr ends at 10 % of Irrm: 190 ns of the 200 ns base
    assert!(close(p.trr.unwrap(), 190e-9, 3e-9), "trr {}", p.trr.unwrap());
    assert!(close(p.qrr.unwrap(), 0.5 * 200e-9 * 4.0, 0.01 * 0.5 * 200e-9 * 4.0), "qrr {}", p.qrr.unwrap());
}


#[test]
fn double_pulse_one_pulse_no_dead_time_and_inverted_current() {
    let dt = 1e-9;
    let n = (8e-6 / dt) as usize;
    let gate = [(0.0, 0.0), (1e-6, 0.0), (1.02e-6, 10.0), (6e-6, 10.0), (6.02e-6, 0.0)];
    let id = [(0.0, 0.0), (1.1e-6, 0.0), (1.2e-6, 10.0), (6.05e-6, 10.0), (6.15e-6, 0.0)];
    let vds = [(0.0, 400.0), (1.05e-6, 400.0), (1.15e-6, 1.0), (6.1e-6, 1.0), (6.2e-6, 400.0), (6.25e-6, 480.0), (6.4e-6, 400.0)];
    let g = sample(|t| pwl(&gate, t), 0.0, dt, n);
    let i = sample(|t| pwl(&id, t), 0.0, dt, n);
    let v = sample(|t| pwl(&vds, t), 0.0, dt, n);
    let pos = oos_dsp::dpt::measure(&v, &i, Some(&g), oos_dsp::dpt::Fourth::Off, dt, 10.0);
    assert!(pos.trr.is_none() && pos.tdt.is_none() && pos.td_on.is_none(), "{}", pos.note);
    assert!(pos.td_off.is_some() && pos.tf.is_some() && pos.toff.is_some(), "{:?}", (pos.td_off, pos.tf, pos.toff));
    let inv: Vec<f32> = i.iter().map(|x| -x).collect();
    let neg = oos_dsp::dpt::measure(&v, &inv, Some(&g), oos_dsp::dpt::Fourth::Off, dt, 10.0);
    assert!(neg.inverted && neg.trr.is_none() && neg.tdt.is_none());
    assert!(close(neg.td_off.unwrap(), pos.td_off.unwrap(), 2e-9), "{} {}", neg.td_off.unwrap(), pos.td_off.unwrap());
    assert!(close(neg.tf.unwrap(), pos.tf.unwrap(), 2e-9));
    assert!(close(neg.toff.unwrap(), pos.toff.unwrap(), 2e-9));
    assert!(neg.td_on.is_none());
}
