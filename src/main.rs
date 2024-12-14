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

fn create_nodes_set(resistors: &[Resistor], v_sources: &[VSource], i_sources: &[ISource]) -> BTreeSet<usize> {
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

    return nodes;
}

fn solve_dc(resistors: &[Resistor], v_sources: &[VSource], i_sources: &[ISource]) -> Vec<f64> {
    let mut nodes = create_nodes_set(resistors, v_sources, i_sources);
    let ref_node = nodes.pop_first().expect("node list is empty");

    let nb_nodes = nodes.len();
    let nb_vsources = v_sources.len();

    let mut g_m = DMatrix::<f64>::zeros(nb_nodes, nb_nodes);
    let mut b_m = DMatrix::<f64>::zeros(nb_nodes, nb_vsources);
    let d_m = DMatrix::<f64>::zeros(nb_vsources, nb_vsources);

    let mut z_v = DVector::<f64>::zeros(nb_nodes + nb_vsources);

    for r in resistors {
        let Resistor(n1, n2, r_val) = *r;

        let g: f64 = 1.0 / f64::from(r_val);

        let n1_idx = nodes.iter().position(|&n| n == n1);
        let n2_idx = nodes.iter().position(|&n| n == n2);

        if let Some(n1_idx) = n1_idx {
            g_m[(n1_idx, n1_idx)] += g;
        }
        if let Some(n2_idx) = n2_idx {
            g_m[(n2_idx, n2_idx)] += g;
        }
        if n1 != ref_node && n2 != ref_node {
            let n1_idx = n1_idx.unwrap();
            let n2_idx = n2_idx.unwrap();

            g_m[(n1_idx, n2_idx)] -= g;
            g_m[(n2_idx, n1_idx)] -= g;
        }
    }

    for (idx, v) in v_sources.iter().enumerate() {
        let VSource(n_neg, n_pos, v_val) = *v;

        let n_neg_idx = nodes.iter().position(|&n| n == n_neg);
        let n_pos_idx = nodes.iter().position(|&n| n == n_pos);

        if let Some(n_neg_idx) = n_neg_idx {
            b_m[(n_neg_idx, idx)] = -1.0;
        }

        if let Some(n_pos_idx) = n_pos_idx {
            b_m[(n_pos_idx, idx)] = 1.0;
        }

        z_v[nb_nodes + idx] = v_val;
    }

    for i in i_sources {
        let ISource(n_neg, n_pos, i_val) = *i;

        let n_neg_idx = nodes.iter().position(|&n| n == n_neg);
        let n_pos_idx = nodes.iter().position(|&n| n == n_pos);

        if let Some(n_neg_idx) = n_neg_idx {
            z_v[n_neg_idx] = -i_val;
        }

        if let Some(n_pos_idx) = n_pos_idx {
            z_v[n_pos_idx] = i_val;
        }
    }
    let c_m = b_m.transpose();

    let a_m = nalgebra::stack![g_m, b_m; c_m, d_m];

    let x_m = a_m.lu().solve(&z_v).unwrap();
    let nodes_v = x_m.as_slice()[0..nb_nodes].to_owned();

    return nodes_v;
}
