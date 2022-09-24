from ova import *
import sympy as sym


if __name__ == '__main__':
    r, theta, z = sym.symbols('r, theta, z')
    a, b, c, d, e, f, g, h, i = sym.symbols('a, b, c, d, e, f, g, h, i')
    x, y, z = sym.symbols('x, y, z')
    cylind = CoordinateSystem(coordinate='Cylindrical', \
                              r=r, t=theta, p=z)
    Ar = sym.Function('A_r')(r, theta, z)
    At = sym.Function('A_t')(r, theta, z)
    Az = sym.Function('A_z')(r, theta, z)
    A = [Ar, At, Az]
    Br = sym.Function('B_r')(r, theta, z)
    Bt = sym.Function('B_t')(r, theta, z)
    Bz = sym.Function('B_z')(r, theta, z)
    B = [Br, Bt, Bz]
    res = cylind.dot(A, B)
    sym.pprint(res)
    T = [[a, b, c], [d, e, f], [g, h, i]]
    C = [x, y, z]
    # benchmarked with mathematica
    sym.pprint(cylind.dot(C, T))
    sym.pprint(cylind.dot(T, C))

    pass
