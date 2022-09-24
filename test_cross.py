from ova import *
import sympy as sym


if __name__ == '__main__':
    r, theta, z = sym.symbols('r, theta, z')
    a, b, c, d, e, f, g, h, i = sym.symbols('a, b, c, d, e, f, g, h, i')
    A = [a, b, c]
    B = [d, e, f]
    cylind = CoordinateSystem(coordinate='Cylindrical', r=r, t=theta, p=z)
    res = cylind.cross(A, B)
    sym.pprint(res)

    pass
