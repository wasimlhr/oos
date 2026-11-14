// Copyright (C) 2026 Julian Decker
// SPDX-License-Identifier: GPL-2.0-or-later

//! DQ0 transform, shaft speed and mechanical power of simulated machine waveforms.

use oos_dsp::sim::sample;
use std::f64::consts::PI;


fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}


#[test]
fn motor_balanced_sines_put_the_peak_on_d() {
    let dt = 1e-4;
    let (f, a, n) = (50.0, 10.0, 8000);
    let w = 2.0 * PI * f;
    let va = sample(|t| a * (w * t).cos(), 0.0, dt, n);
    let vb = sample(|t| a * (w * t - 2.0 * PI / 3.0).cos(), 0.0, dt, n);
    let vc = sample(|t| a * (w * t + 2.0 * PI / 3.0).cos(), 0.0, dt, n);
    let d = oos_dsp::motor::dq0(&va, &vb, &vc, dt, oos_dsp::motor::Theta::Vector).unwrap();
    assert!(close(d.d, a, 0.05), "{:?}", (d.d, d.q, d.zero));
    assert!(close(d.q, 0.0, 0.05) && close(d.zero, 0.0, 0.05), "{:?}", (d.d, d.q, d.zero));
    assert!(oos_dsp::motor::dq0(&va, &vb, &[], dt, oos_dsp::motor::Theta::Vector).is_none());
}


#[test]
fn motor_quadrature_speed_and_swapped_direction() {
    let (hold, cycles, dt, cpr) = (10usize, 50usize, 1e-4, 100.0);
    let n = hold * 4 * cycles;
    let mut a = vec![false; n];
    let mut b = vec![false; n];
    for k in 0..n {
        match (k / hold) % 4 {
            0 => {}
            1 => a[k] = true,
            2 => {
                a[k] = true;
                b[k] = true;
            }
            _ => b[k] = true,
        }
    }
    let s = oos_dsp::motor::shaft_quad(&a, &b, None, dt, 0.0, cpr, false);
    // one count per state change: the first sample establishes the state and does not count
    let counts = (n / hold - 1) as f64;
    let cps = counts / ((n - 1) as f64 * dt);
    let rpm = cps / cpr * 60.0;
    assert!(close(s.rpm.unwrap(), rpm, 1e-6), "{} vs {} ({})", s.rpm.unwrap(), rpm, cps * 360.0 / cpr);
    assert_eq!(s.direction, oos_dsp::motor::Direction::Forward);
    let swapped = oos_dsp::motor::shaft_quad(&a, &b, None, dt, 0.0, cpr, true);
    assert_eq!(swapped.direction, oos_dsp::motor::Direction::Reverse);
    assert!(close(swapped.rpm.unwrap(), -rpm, 1e-6), "{}", swapped.rpm.unwrap());
}


#[test]
fn motor_torque_and_bus_give_mechanical_power_and_efficiency() {
    let (rpm, dt, n) = (600.0, 1e-4, 4000);
    let deg_s = rpm * 6.0;
    let v = sample(|t| deg_s * t, 0.0, dt, n);
    let shaft = oos_dsp::motor::shaft_analog(&v, dt, 0.0, 1.0, 0.0, false, 1);
    assert!(close(shaft.rpm.unwrap(), rpm, 0.05), "{}", shaft.rpm.unwrap());
    let torque = vec![2.0f32; n];
    let bus_v = vec![48.0f32; n];
    let bus_i = vec![4.0f32; n];
    let m = oos_dsp::motor::mechanical(shaft, Some(&torque), 1.0, Some(&bus_v), Some(&bus_i));
    let p = 2.0 * rpm * 2.0 * PI / 60.0;
    assert!(close(m.p_mech.unwrap(), p, 0.01), "{}", m.p_mech.unwrap());
    assert!(close(m.p_elec.unwrap(), 192.0, 1e-6));
    assert!(close(m.eff.unwrap(), p / 192.0, 1e-4), "{}", m.eff.unwrap());
    assert!(!m.generating);
}
