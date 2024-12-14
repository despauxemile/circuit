use std::collections::BTreeSet;

use nalgebra::{DMatrix, DVector};

#[derive(Debug)]
struct Resistor (usize, usize, u32);

#[derive(Debug)]
struct VSource (usize, usize, f64);

#[derive(Debug)]
struct ISource (usize, usize, f64);

#[derive(Debug)]
enum Element {
    Resistor(Resistor),
    VSource(VSource),
    ISource(ISource),
}

impl Element {
    fn get_nodes(&self) -> Vec<usize> {
        match self {
            Element::Resistor(resistor) => vec!(resistor.0, resistor.1),
            Element::VSource(vsource) => vec!(vsource.0, vsource.1),
            Element::ISource(isource) => vec!(isource.0, isource.1),
        }
    }
}

impl From<Resistor> for Element {
    fn from(value: Resistor) -> Self {
        Self::Resistor(value)
    }
}
impl From<VSource> for Element {
    fn from(value: VSource) -> Self {
        Self::VSource(value)
    }
}
impl From<ISource> for Element {
    fn from(value: ISource) -> Self {
        Self::ISource(value)
    }
}

fn main() {

    let elts: Vec<Element> = Vec::from([
        Resistor(1, 2, 5).into(),
        Resistor(2, 0, 3).into(),
        Resistor(3, 0, 10).into(),
        VSource(0, 1, 5.0).into(),
        VSource(2, 3, 10.0).into(),
        ISource(1, 2, 2.0).into(),
    ]);

    let res = solve_dc(&elts);

    println!("node voltages: {:?}", res);

}

fn create_nodes_set(elements: &[Element]) -> BTreeSet<usize> {
    // node population generation
    let mut nodes = BTreeSet::<usize>::new();

    elements.iter().for_each(|elt| {
        elt.get_nodes().iter().for_each(|&n| {
            nodes.insert(n);
        })
    });

    return nodes;
}

fn solve_dc(elements: &[Element]) -> Vec<f64> {
    let mut nodes = create_nodes_set(elements);
    let ref_node = nodes.pop_first().expect("node list is empty");

    let nb_elt_nodes = elements.iter().filter(|&elt| {
        match elt {
            Element::Resistor(_) => true,
            _ => false,
        }
    }).collect::<Vec<_>>().len();

    let nb_vsources = elements.iter().filter(|&elt| {
        match elt {
            Element::VSource(_) => true,
            _ => false,
        }
    }).collect::<Vec<_>>().len();

    let mut g_m = DMatrix::<f64>::zeros(nb_elt_nodes, nb_elt_nodes);
    let mut b_m = DMatrix::<f64>::zeros(nb_elt_nodes, nb_vsources);
    let d_m = DMatrix::<f64>::zeros(nb_vsources, nb_vsources);

    let mut z_v = DVector::<f64>::zeros(nb_elt_nodes + nb_vsources);

    let mut vsource_idx = 0usize;
    for elt in elements {

        match elt {
            Element::Resistor(resistor) => {
                let Resistor(n1, n2, r_val) = *resistor;

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
            },
            Element::VSource(vsource) => {
                let VSource(n_neg, n_pos, v_val) = *vsource;

                let n_neg_idx = nodes.iter().position(|&n| n == n_neg);
                let n_pos_idx = nodes.iter().position(|&n| n == n_pos);

                if let Some(n_neg_idx) = n_neg_idx {
                    b_m[(n_neg_idx, vsource_idx)] = -1.0;
                }

                if let Some(n_pos_idx) = n_pos_idx {
                    b_m[(n_pos_idx, vsource_idx)] = 1.0;
                }

                z_v[nb_elt_nodes + vsource_idx] = v_val;

                vsource_idx += 1;
            },
            Element::ISource(isource) => {
                let ISource(n_neg, n_pos, i_val) = *isource;

                let n_neg_idx = nodes.iter().position(|&n| n == n_neg);
                let n_pos_idx = nodes.iter().position(|&n| n == n_pos);

                if let Some(n_neg_idx) = n_neg_idx {
                    z_v[n_neg_idx] = -i_val;
                }

                if let Some(n_pos_idx) = n_pos_idx {
                    z_v[n_pos_idx] = i_val;
                }
            },
            //_ => todo!()
        }

    }

    let c_m = b_m.transpose();

    let a_m = nalgebra::stack![g_m, b_m; c_m, d_m];

    let x_m = a_m.lu().solve(&z_v).unwrap();
    let nodes_v = x_m.as_slice()[0..nodes.len()].to_owned();

    return nodes_v;
}
