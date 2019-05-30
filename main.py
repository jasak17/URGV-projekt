import numpy as np
import argparse
import math

def trilateration3(P1, P2, P3, r1, r2, r3):
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
    print("\nprimer za 3 postaje: ")
    return K1, 0



def trilateration4(v1, v2, v3, v4):
    p1 = np.asarray(v1[0])
    p2 = np.asarray(v2[0])
    p3 = np.asarray(v3[0])
    p4 = np.asarray(v4[0])
    r1 = np.asarray(v1[1])
    r2 = np.asarray(v2[1])
    r3 = np.asarray(v3[1])
    r4 = np.asarray(v4[1])
    e_x = (p2 - p1) / np.linalg.norm(p2 - p1)
    i = np.dot(e_x, (p3 - p1))
    e_y = (p3 - p1 - (i * e_x)) / (np.linalg.norm(p3 - p1 - (i * e_x)))
    e_z = np.cross(e_x, e_y)
    d = np.linalg.norm(p2 - p1)
    j = np.dot(e_y, (p3 - p1))
    x = ((r1 ** 2) - (r2 ** 2) + (d ** 2)) / (2 * d)
    y = (((r1 ** 2) - (r3 ** 2) + (i ** 2) + (j ** 2)) / (2 * j)) - ((i / j) * (x))
    z1 = np.sqrt(r1 ** 2 - x ** 2 - y ** 2)
    z2 = np.sqrt(r1 ** 2 - x ** 2 - y ** 2) * (-1)
    ans1 = p1 + (x * e_x) + (y * e_y) + (z1 * e_z)
    ans2 = p1 + (x * e_x) + (y * e_y) + (z2 * e_z)
    dist1 = np.linalg.norm(p4 - ans1)
    dist2 = np.linalg.norm(p4 - ans2)
    if np.abs(r4 - dist1) < np.abs(r4 - dist2):
        return ans1
    else:
        return ans2


if __name__ == "__main__":
    r1 = 0.9661
    r2 = 2.2039
    r3 = 1.7243
    r4 = 4.0825
    r5 = 4.7003
    r6 = 6.7794

    vector1 = np.array([[2, 1, 0], [r1]])
    vector2 = np.array([[4, 3, 0], [r2]])
    vector3 = np.array([[3, 3, 1], [r3]])
    vector4 = np.array([[2, 3, 4], [r4]])
    vector5 = np.array([[5, 5, 2], [r5]])
    vector6 = np.array([[3, 8, 3], [r6]])

    # for i in range(1,4):

    rez = trilateration4(vector1, vector2, vector3, vector4)
    rez1 = trilateration4(vector2, vector3, vector4, vector5)
    rez2 = trilateration4(vector3, vector1, vector5, vector6)
    rez3 = trilateration4(vector4, vector5, vector6, vector1)
    rez4 = trilateration4(vector5, vector6, vector1, vector2)
    rez5 = trilateration4(vector6, vector1, vector2, vector3)

    avg = (rez + rez1 + rez2 + rez3 + rez4 + rez5) / 6
    print(avg)
    err = ((rez - avg) ** 2 + (1 - avg) ** 2 + (rez2 - avg) ** 2 + (rez3 - avg) ** 2 + (rez4 - avg) ** 2 + (
                rez5 - avg) ** 2)/6

    print(math.sqrt(err[0]))
    print(math.sqrt(err[1]))
    print(math.sqrt(err[2]))

    p = argparse.ArgumentParser()
    p.add_argument('number', type=float, nargs='+')
    args = p.parse_args()

    vector1, vector2, vector3 = np.array([args.number[0], args.number[1], args.number[2]]), np.array(
        [args.number[4], args.number[5], args.number[6]]), np.array([args.number[8], args.number[9], args.number[10]])
    R1 = args.number[3]
    R2 = args.number[7]
    R3 = args.number[11]

    print(trilateration3(vector1, vector2, vector3, R1, R2, R3))