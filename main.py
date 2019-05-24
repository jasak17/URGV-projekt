import numpy as np


def trilateration(P1, P2, P3, P4):
    p1 = np.array([0, 0, 0])
    p2 = np.array([P2[0][0] - P1[0][0], P2[0][1] - P1[0][1], P2[0][2] - P1[0][2]])
    p3 = np.array([P3[0][0] - P1[0][0], P3[0][1] - P1[0][1], P3[0][2] - P1[0][2]])

    r1 = np.asarray(P1[1])
    r2 = np.asarray(P2[1])
    r3 = np.asarray(P3[1])
    r4 = np.asarray(P4[1])



    v1 = p2 - p1
    v2 = p3 - p1

    Xn = (v1)/np.linalg.norm(v1)

    tmp = np.cross(v1, v2)

    Zn = (tmp)/np.linalg.norm(tmp)

    Yn = np.cross(Xn, Zn)

    a = np.dot((np.asarray(P4[0]) - np.asarray(P1[0])), Xn)
    b = np.dot((np.asarray(P4[0]) - np.asarray(P1[0])), Yn)
    c = np.dot((np.asarray(P4[0]) - np.asarray(P1[0])), Zn)

    p4 = np.array([a, b, c])

    i = np.dot(Xn, v2)
    d = np.dot(Xn, v1)
    j = np.dot(Yn, v2)

    X = ((r1**2)-(r2**2)+(d**2))/(2*d)
    Y = (((r1**2)-(r3**2)+(i**2)+(j**2))/(2*j))-((i/j)*(X))

    Z = (r1**2 - r4**2 + a**2 + b**2 + c**2) / 2*c - (a/c)*X - (b/c)*Y

    # Z1 = np.sqrt(max(0, r1 ** 2 - X ** 2 - Y ** 2))
    Z2 = -Z

    K1 = P1[0] + X*Xn + Y * Yn + Z * Zn
    K2 = p1 + X * Xn + Y * Yn - Z2 * Zn

    # K2[2] = 0
    return K1


if __name__ == "__main__":
    # PRIMER 1
    # import argparse
    # #
    # # p = argparse.ArgumentParser()
    # # p.add_argument('number', type=float, nargs='+')
    # # args = p.parse_args()

    # vector1, vector2, vector3, vector4 = np.array([args.number[0], args.number[1], args.number[2]]), np.array([args.number[4], args.number[5], args.number[6]]), np.array([args.number[8], args.number[9], args.number[10]]), np.array([args.number[12], args.number[13], args.number[14]])
    # R1 = args.number[3]
    # R2 = args.number[7]
    # R3 = args.number[11]
    # R4 = args.number[15]
    # R5 = args.number[19]
    # R6 = args.number[23]
    # R7 = args.number[27]
    # R8 = args.number[31]
    # R9 = args.number[35]

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

    print(trilateration(vector1, vector2, vector3, vector4))





# def trilateration(p1, p2, p3, p4, r1, r2, r3, r4):
#     e_x = (p2 - p1) / np.linalg.norm(p2 - p1)
#     i = np.dot(e_x, (p3 - p1))
#     e_y = (p3 - p1 - (i * e_x)) / (np.linalg.norm(p3 - p1 - (i * e_x)))
#     e_z = np.cross(e_x, e_y)
#     d = np.linalg.norm(p2 - p1)
#     j = np.dot(e_y, (p3 - p1))
#     x = ((r1 ** 2) - (r2 ** 2) + (d ** 2)) / (2 * d)
#     y = (((r1 ** 2) - (r3 ** 2) + (i ** 2) + (j ** 2)) / (2 * j)) - ((i / j) * (x))
#     z1 = np.sqrt(r1 ** 2 - x ** 2 - y ** 2)
#     z2 = np.sqrt(r1 ** 2 - x ** 2 - y ** 2) * (-1)
#     ans1 = p1 + (x * e_x) + (y * e_y) + (z1 * e_z)
#     ans2 = p1 + (x * e_x) + (y * e_y) + (z2 * e_z)
#     dist1 = np.linalg.norm(p4 - ans1)
#     dist2 = np.linalg.norm(p4 - ans2)
#     if np.abs(r4 - dist1) < np.abs(r4 - dist2):
#         return ans1
#     else:
#         return ans2
