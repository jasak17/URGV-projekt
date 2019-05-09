import numpy as np


def trilateration(P1, P2, P3, r1, r2, r3):
    p1 = np.array([0, 0, 0])
    p2 = np.array([P2[0] - P1[0], P2[1] - P1[1], P2[2] - P1[2]])
    p3 = np.array([P3[0] - P1[0], P3[1] - P1[1], P3[2] - P1[2]])

    v1 = p2 - p1
    v2 = p3 - p1

    Xn = (v1)/np.linalg.norm(v1)

    tmp = np.cross(v1, v2)

    Zn = (tmp)/np.linalg.norm(tmp)

    Yn = np.cross(Xn, Zn)

    i = np.dot(Xn, v2)
    d = np.dot(Xn, v1)
    j = np.dot(Yn, v2)

    X = ((r1**2)-(r2**2)+(d**2))/(2*d)
    Y = (((r1**2)-(r3**2)+(i**2)+(j**2))/(2*j))-((i/j)*(X))
    Z1 = np.sqrt(max(0, r1 ** 2 - X ** 2 - Y ** 2))
    Z2 = -Z1

    K1 = P1 + X*Xn + Y * Yn + Z1 * Zn
    K2 = p1 + X * Xn + Y * Yn - Z2 * Zn
    K1[2] = 0
    K2[2] = 0
    return np.round(K1, 0), np.round(K2, 0)


if __name__ == "__main__":
    # PRIMER 1
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument('number', type=float, nargs='+')
    args = p.parse_args()


    vector1, vector2, vector3 = np.array([args.number[0], args.number[1], args.number[2]]), np.array([args.number[4], args.number[5], args.number[6]]), np.array([args.number[8], args.number[9], args.number[10]])
    R1 = args.number[3]
    R2 = args.number[7]
    R3 = args.number[11]

    print(trilateration(vector1, vector2, vector3, R1, R2, R3))


