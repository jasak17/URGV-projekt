import numpy as np
import math



def trilateration(P1, P2, P3, R1, R2, R3):
    p1 = np.array([0, 0, 0])
    p2 = np.array([P2[0] - P1[0], P2[1] - P1[1], P2[2] - P1[2]])
    p3 = np.array([P3[0] - P1[0], P3[1] - P1[1], P3[2] - P1[2]])

    v1 = p2 - p1
    v2 = p3 - p1

    Xn = v1 / np.sqrt(np.sum(v1**2))

    Z = np.cross(v1, v2)
    Zn = Z / np.sqrt(np.sum(Z ** 2))

    Yn = np.cross(Xn, Zn)

    d = np.dot(Xn, v1)
    i = np.dot(Xn, v2)
    j = np.dot(Yn, v2)

    X = (((R1 ** 2) - (R2 ** 2) + (d ** 2))) / (2 * d)
    Y = ((R1 ** 2 - R3 ** 2 + i ** 2+ j**2)) / ((2 * j)-(i/j)*X)

    a = abs(R1**2 - X ** 2 - Y ** 2)

    Z1 = math.sqrt(a)
    Z2 = -math.sqrt(abs(R1 ** 2 - X ** 2 - Y ** 2))

    K1 = P1 + X*Xn + Y * Yn + Z1 * Zn
    K2 = P1 + X * Xn + Y * Yn - Z2 * Zn



if __name__ == "__main__":
    # PRIMER 1
    vector1, vector2, vector3 = np.array([2, 1, 0]), np.array([4, 3, 0]), np.array([4, 4, 1])
    r1 = 2
    r2 = 2
    r3 = 2.449

    trilateration(vector1, vector2, vector3, r1, r2, r3)
