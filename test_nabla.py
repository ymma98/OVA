from ova import *
import sympy as sym


if __name__ == '__main__':
    # test grad
    r = sym.symbols('r', positive=True)
    theta, z = sym.symbols('theta, z')

    cylind = CoordinateSystem(coordinate='Cylindrical', r=r, t=theta, p=z)
    Ur = sym.Function('Ur')(r, theta, z)
    Ut = sym.Function('Ut')(r, theta, z)
    Uz = sym.Function('Uz')(r, theta, z)
    U = [Ur, Ut, Uz]
    phi = sym.Function('phi')(r, theta, z)
    print("**************************************")
    # grad of phi
    print(" grad of phi, phi is a scalar function")
    sym.pprint(cylind.grad(phi))
    print("**************************************")
    # grad u
    print(" grad of vector u, u is a vector")
    sym.pprint(cylind.grad(U))
    print("**************************************")
    # u dot grad(u)
    print(" u dot grad u, u is a vector")
    sym.pprint(cylind.dot(U, cylind.grad(U)))
    print("**************************************")
    # curl of u
    print("curl of u, u is a vector")
    sym.pprint(cylind.curl(U))
    print("**************************************")
    sym.pprint(cylind.div(U))
    Trr = sym.Function('Trr')(r, theta, z)
    Trt = sym.Function('Trt')(r, theta, z)
    Trz = sym.Function('Trz')(r, theta, z)
    Ttr = sym.Function('Ttr')(r, theta, z)
    Ttt = sym.Function('Ttt')(r, theta, z)
    Ttz = sym.Function('Ttz')(r, theta, z)
    Tzr = sym.Function('Tzr')(r, theta, z)
    Tzt = sym.Function('Tzt')(r, theta, z)
    Tzz = sym.Function('Tzz')(r, theta, z)
    T = sym.Matrix([[Trr, Trt, Trz], [Ttr, Ttt, Ttz], [Tzr, Tzt, Tzz]])
    print("div of tensor T, T is a 2-rank tensor")
    sym.pprint(cylind.div(T))
    print("**************************************")
    print("Laplacian of phi, phi is a scalar function")
    sym.pprint(cylind.Laplacian(phi))
