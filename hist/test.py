from ova import *
import sympy as sym


if __name__ == '__main__':
    # test calc co_metirc matrix
    r, theta, z = sym.symbols('r, theta, z')
    x_expr = r * sym.cos(theta)
    y_expr = r * sym.sin(theta)
    z_expr = z
    co_metric = get_co_metric(x_expr, y_expr, z_expr, r, theta, z)
    sym.pprint(co_metric)

    # test calc contra_metric matrix
    x, y, z = sym.symbols('x, y, z')
    r_expr = sym.sqrt(x**2 + y**2)
    theta_expr = sym.asin(y / sym.sqrt(x**2 + y**2))
    z_expr = z
    contra_metric = get_contra_metric(r_expr, theta_expr, z_expr, r, theta, z)
    sym.pprint(contra_metric)

    pass
