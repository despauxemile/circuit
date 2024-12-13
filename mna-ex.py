import numpy as np

# R * I = V
# G * V = I

def main():
    R1 = 2
    R2 = 4
    R3 = 1
    Va = 3
    Vb = 2

    G1 = 1/R1
    G2 = 1/R2
    G3 = 1/R3

    A = np.array([
        [G1, -G1, 0, 1, 0],
        [-G1, G1+G2+G3, -G2, 0, 0],
        [0, -G2, G2, 0, 1],
        [1, 0, 0, 0, 0],
        [0, 0, 1, 0, 0]
    ])

    b = np.array([0, 0, 0, Va, Vb])

    x = np.linalg.solve(A, b)

    print(x)
    print(A)
        
if __name__ == "__main__":
    main()
