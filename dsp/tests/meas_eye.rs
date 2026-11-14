// Copyright (C) 2026 Julian Decker
// SPDX-License-Identifier: GPL-2.0-or-later

//! Eye and jitter additions: strobe-clock folding, BER extrapolation, sub-rate jitter, common mode, SSC.

use oos_dsp::eye::{self, Eye};
use std::f64::consts::PI;


fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

/// PRBS7 (x^7 + x^6 + 1) bit sequence.
fn prbs7(n: usize) -> Vec<bool> {
    let mut s = 0x7fu8;
    (0..n)
        .map(|_| {
            let b = ((s >> 6) ^ (s >> 5)) & 1;
            s = ((s << 1) | b) & 0x7f;
            b == 1
        })
        .collect()
}

/// A source-synchronous strobe (DQS) clocks the data (DQ). The eye is DQ folded on DQS edges, not CDR from DQ.
#[test]
fn eye_folds_data_on_a_strobe_clock() {
    let (rate, spu, nui) = (1e9, 32usize, 2048usize);
    let ui = 1.0 / rate;
    let dt = ui / spu as f64;
    let bits = prbs7(nui);
    let mut dq = Vec::with_capacity(nui * spu);
    let mut dqs = Vec::with_capacity(nui * spu);
    for i in 0..nui * spu {
        let t = i as f64 * dt;
        let k = ((t + 0.5 * ui) / ui).floor() as usize;
        dq.push(if bits[k.min(bits.len() - 1)] { 1.0 } else { 0.0 });
        dqs.push(if (t / ui).floor() as usize % 2 == 0 { 1.0 } else { 0.0 });
    }
    let clk = eye::recover_clock(&dqs, dt, 0.5, 0.1, None).unwrap();
    assert!(close(clk.period, ui, ui * 1e-6), "{}", clk.period);
    let mut e = Eye::new(128, 64, 2, -0.25, 1.25);
    e.accumulate(&dq, dt, &eye::center_strobe(&clk));
    let m = e.measure(clk.period);
    assert!(m.height > 0.7, "{:?}", m);
    assert!(m.width > 0.5 * ui, "{:?}", m);
    let idle = eye::gate_around_clock(&dq, dt, &clk, 1.5);
    assert_eq!(idle.iter().filter(|x| x.is_finite()).count(), dq.len());
}


#[test]
fn eye_ber_subrate_common_mode_and_ssc() {
    use oos_dsp::eye::{self, Clock, EyeMeas, Jitter};
    let m = EyeMeas { one: 1.0, zero: 0.0, sigma_one: 0.01, sigma_zero: 0.01, ..Default::default() };
    let q = eye::norm_inv(1.0 - 1e-12);
    let want = (1.0 - q * 0.01) - q * 0.01;
    assert!(close(m.height_at(1e-12), want, 1e-9), "{}", m.height_at(1e-12));
    let closed = EyeMeas { sigma_one: 1.0, sigma_zero: 1.0, ..m };
    assert_eq!(closed.height_at(1e-12), 0.0);
    let j = Jitter { dj_pp: 1e-12, rj_rms: 1e-12, ..Default::default() };
    let j2 = eye::tj_at(&j, eye::BER_J2);
    let j9 = eye::tj_at(&j, eye::BER_J9);
    assert!(j2 < j9 && close(eye::tj_at(&j, 1e-12), 1e-12 + 2.0 * q * 1e-12, 1e-18), "{} {}", j2, j9);
    assert!(close(eye::eye_width_at(&j, 1e-9, 1e-12), (1.0 - eye::tj_at(&j, 1e-12) / 1e-9).max(0.0), 1e-12));

    let period = 1e-9;
    let trend: Vec<f32> = (0..256).map(|k| (1e-12 * (PI * k as f64).cos()) as f32).collect();
    assert!(close(eye::subrate_pp(&trend, period, 2), 2e-12, 1e-16), "{}", eye::subrate_pp(&trend, period, 2));
    let tone = |div: f64| -> Vec<f32> { (0..512).map(|k| (1e-12 * (PI * k as f64 / (div / 2.0)).cos()) as f32).collect() };
    assert!(close(eye::subrate_pp(&tone(4.0), period, 4), 2e-12, 1e-14), "F/4");
    assert!(close(eye::subrate_pp(&tone(8.0), period, 8), 2e-12, 1e-14), "F/8");
    assert!(close(eye::srj_pp(&trend, period), 2e-12, 1e-15), "{}", eye::srj_pp(&trend, period));
    let wide = tone(2.0);
    let mixed: Vec<f32> = (0..512).map(|k| wide[k] + tone(4.0)[k] + tone(8.0)[k]).collect();
    assert!(close(eye::srj_pp(&mixed, period), 6e-12, 0.05), "{}", eye::srj_pp(&mixed, period));
    assert_eq!(eye::phase_noise(&vec![0.0; 128], period, 1.0, 1.0 / period), 0.0);
    let pn: Vec<f32> = (0..4096).map(|k| (5e-12 * (2.0 * PI * k as f64 / 64.0).sin()) as f32).collect();
    let floor = eye::phase_noise(&pn, period, 1.0 / (period * 256.0), 1.0 / (period * 16.0));
    // the sine's own bin is left out, so this is the window skirt, well under the 5 ps tone
    assert!(floor.is_finite() && floor < 2.5e-12, "phase noise off the line {floor}");
    assert!(eye::npj(&[0.0; 10], &Clock::default(), &trend, &eye::EdgeMeans::default()).is_nan());
    let n_edge = 2048usize;
    let sine: Vec<f64> = (0..n_edge).map(|k| 10e-12 * (2.0 * PI * k as f64 / 32.0).sin()).collect();
    let trend_j: Vec<f32> = sine.iter().map(|x| *x as f32).collect();
    let edges: Vec<f64> = (0..n_edge).map(|k| k as f64 * period).collect();
    let clk = Clock { period, phase: 0.0, edges, ticks: vec![], first_rising: true };
    let left = eye::npj(&sine, &clk, &trend_j, &eye::EdgeMeans::default());
    assert!(left.is_finite() && left < 2e-12, "npj after the sine {left}");
    let ps = eye::psij(&trend_j, &clk, &sine, None, 1.0 / (period * 32.0), 0.4, 0.6 * period, None);
    assert!(close(ps.before.tie_pp, 20e-12, 0.02), "psij before {}", ps.before.tie_pp);
    assert!(ps.after.tie_pp < 2e-12 && ps.after.pj < ps.before.pj, "psij after {} pj {}", ps.after.tie_pp, ps.after.pj);

    assert!(close(eye::dc_common(&[1.0, 1.0], &[3.0, 3.0]), 2.0, 1e-9));
    assert!(close(eye::ac_common_pp(&[0.0, 1.0, 0.0, 1.0], None, 1e-6, false), 1.0, 1e-9));
    let a: Vec<f32> = (0..11).map(|k| k as f32 / 10.0).collect();
    let b: Vec<f32> = (0..11).map(|k| 1.0 - k as f32 / 10.0).collect();
    assert!(close(eye::diff_crossover(&a, &b), 0.5, 1e-6), "{}", eye::diff_crossover(&a, &b));

    let dt = 1e-8;
    let mut v = vec![0.5f32; 500];
    for k in 40..60 {
        v[k] = 1.0;
    }
    let c = Clock { period: 1e-6, phase: 0.0, edges: vec![0.0, 3e-6], ticks: vec![], first_rising: true };
    assert!(close(eye::tnt_ratio(&v, dt, &c), 20.0 * (0.5f64).log10(), 0.05), "{}", eye::tnt_ratio(&v, dt, &c));

    let mut edges = Vec::new();
    let mut t = 0.0;
    for _ in 0..64 {
        let f = 1e6 * (1.0 + 0.05 * (2.0 * PI * 125e3 * t).sin());
        let p = 1.0 / f;
        edges.push(t);
        edges.push(t + 0.5 * p);
        t += p;
    }
    let c = Clock { period: 1e-6, phase: 0.0, edges, ticks: vec![], first_rising: true };
    let (ppm, _) = eye::ssc_deviation(&c);
    assert!(ppm > 50_000.0 && ppm < 150_000.0, "{}", ppm);
    let rate = eye::ssc_rate(&c);
    assert!(close(rate, 125e3, 30e3), "{}", rate);
    let clock = Clock { period: 1e-6, phase: 0.0, edges: (0..32).map(|k| k as f64 * 1e-6).collect(), ticks: vec![], first_rising: true };
    assert!(eye::is_clock(&clock));
    let data = Clock { edges: (0..32).map(|k| k as f64 * 3e-6).collect(), ..clock };
    assert!(!eye::is_clock(&data));
}
