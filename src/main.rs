use std::collections::BTreeSet;

use nalgebra::{DMatrix, DVector};

#[derive(Debug)]
struct Resistor (usize, usize, u32);

#[derive(Debug)]
struct VSource (usize, usize, f64);

#[derive(Debug)]
struct ISource (usize, usize, f64);

fn main() {

    let rs = Vec::from([
        Resistor(1, 2, 5),
        Resistor(2, 0, 3),
        Resistor(3, 0, 10),
    ]);

    let vs = Vec::from([
        VSource(0, 1, 5.0),
        VSource(2, 3, 10.0),
    ]);

    let is = Vec::from([
        ISource(1, 2, 2.0),
    ]);

    let res = solve_dc(&rs, &vs, &is);

    println!("node voltages: {:?}", res);

}

fn solve_dc(resistors: &[Resistor], v_sources: &[VSource], i_sources: &[ISource]) -> Vec<f64> {
    // node population generation
    let mut nodes = BTreeSet::<usize>::new();

    for r in resistors {
        nodes.insert(r.0);
        nodes.insert(r.1);
    }
    
    for vs in v_sources {
        nodes.insert(vs.0);
        nodes.insert(vs.1);
    }

    for is in i_sources {
        nodes.insert(is.0);
        nodes.insert(is.1);
    }

    let ref_node = nodes.first().unwrap();

    let nb_node_rows = nodes.len() - 1;
    let nb_vsources = v_sources.len();

    let mut g_m = DMatrix::<f64>::zeros(nb_node_rows, nb_node_rows);
    let mut b_m = DMatrix::<f64>::zeros(nb_node_rows, nb_vsources);
    let d_m = DMatrix::<f64>::zeros(nb_vsources, nb_vsources);

    let mut z_v = DVector::<f64>::zeros(nb_node_rows + nb_vsources);

    for r in resistors {
        let Resistor(n1, n2, r_val) = r;

        let g: f64 = 1.0 / f64::from(*r_val);

        let n1_idx = nodes.iter().position(|n| n == n1).unwrap();
        let n2_idx = nodes.iter().position(|n| n == n2).unwrap();

        if n1 != ref_node {
            g_m[(n1_idx-1, n1_idx-1)] += g;
        }
        if n2 != ref_node {
            g_m[(n2_idx-1, n2_idx-1)] += g;
        }
        if n1 != ref_node && n2 != ref_node {
            g_m[(n1_idx-1, n2_idx-1)] -= g;
            g_m[(n2_idx-1, n1_idx-1)] -= g;
        }
    }

    for (idx, v) in v_sources.iter().enumerate() {
        let VSource(n_neg, n_pos, v_val) = v;

        let n_neg_idx = nodes.iter().position(|n| n == n_neg).unwrap();
        let n_pos_idx = nodes.iter().position(|n| n == n_pos).unwrap();

        if n_neg != ref_node {
            b_m[(n_neg_idx-1, idx)] = -1.0;
        }

        if n_pos != ref_node {
            b_m[(n_pos_idx-1, idx)] = 1.0;
        }

        z_v[nb_node_rows + idx] = *v_val;
    }

    for i in i_sources {
        let ISource(n_neg, n_pos, i_val) = i;

        let n_neg_idx = nodes.iter().position(|n| n == n_neg).unwrap();
        let n_pos_idx = nodes.iter().position(|n| n == n_pos).unwrap();

        if n_neg != ref_node {
            z_v[n_neg_idx-1] = -i_val;
        }

        if n_pos != ref_node {
            z_v[n_pos_idx-1] = *i_val;
        }
    }
    let c_m = b_m.transpose();

    let a_m = nalgebra::stack![g_m, b_m; c_m, d_m];

    let x_m = a_m.lu().solve(&z_v).unwrap();
    let nodes_v = x_m.as_slice()[0..nb_node_rows].to_owned();

    return nodes_v;
}
