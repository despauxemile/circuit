import numpy as np

def solve_DC(
    resistors: list[tuple[int, int, int]],
    v_sources: list[tuple[int, int, int]],
    i_sources: list[tuple[int, int, int]],
) -> np.ndarray:
    # generate nodes list
    nodes = set()

    for r in resistors:
        nodes.update([ r[0], r[1] ])

    for v in v_sources:
        nodes.update([ v[0], v[1] ])

    for i in i_sources:
        nodes.update([ i[0], i[1] ])

    ref_node = next(iter(nodes))
    nodes.remove(ref_node)

    node_idx = lambda n: list(nodes).index(n) if n != ref_node else -1

    # MNA
    # A@x = z
    # [ G B ] * [ X ] = [ z ]
    # [ C D ]   [ X ] = [ z ]
    # G: conductance matrix
    # B: voltage sources connection
    # C: B.T
    # D: zero
    # X: unknowns
    # z: sources

    n_n = len(nodes)
    n_vs = len(v_sources)

    G_m = np.zeros((n_n, n_n))
    B_m = np.zeros((n_n, n_vs))
    D_m = np.zeros((n_vs, n_vs))

    z_m = np.zeros(n_n + n_vs)

    # generate G_m (R)
    for n1, n2, R in resistors:
        G = 1/R

        n1_idx = node_idx(n1)
        n2_idx = node_idx(n2)

        if n1 != ref_node:
            G_m[n1_idx][n1_idx] += G

        if n2 != ref_node:
            G_m[n2_idx][n2_idx] += G

        if n1 != ref_node and n2 != ref_node:
            G_m[n1_idx][n2_idx] -= G
            G_m[n2_idx][n1_idx] -= G

    # generate B_m + populate z_m (V)
    for v_idx, (n_neg, n_pos, V) in enumerate(v_sources):
        n_neg_idx = node_idx(n_neg)
        n_pos_idx = node_idx(n_pos)

        if n_neg != ref_node:
            B_m[n_neg_idx][v_idx] = -1

        if n_pos != ref_node:
            B_m[n_pos_idx][v_idx] = 1

        z_m[n_n + v_idx] = V

    # populate z_m (I)
    for n_neg, n_pos, I in i_sources:
        n_neg_idx = node_idx(n_neg)
        n_pos_idx = node_idx(n_pos)

        if n_neg != ref_node:
            z_m[n_neg_idx] = -I

        if n_pos != ref_node:
            z_m[n_pos_idx] = I

    C_m = B_m.T

    A_m = np.vstack([
        np.hstack([G_m, B_m]),
        np.hstack([C_m, D_m])
    ])

    print(A_m)
    print(z_m)

    x = np.linalg.solve(A_m, z_m)

    return x

def main():
    # R : (n1, n2, R)
    R1 = (1, 2, 5)
    R2 = (2, 0, 3)
    R3 = (3, 0, 10)
    Rs = [R1, R2, R3]

    # V : (n-, n+, V)
    V1 = (0, 1, 5)
    V2 = (2, 3, 10)
    Vs = [V1, V2]

    I1 = (1, 2, 2)
    Is = [I1]

    out = solve_DC(Rs, Vs, Is)
    print(out)
        
if __name__ == "__main__":
    main()
