"""
Orthogonal curvilinear coordinates Vector Analysis (OVA)
ymma98@qq.com
"""

import sympy as sym


def def_co_metric(r_expr, t_expr, p_expr):
    """
    assuming the new coordinate has 
    :param r_expr: expression 
    """


class CoordinateSystem:
    """
    defines a orthogonal curvilinear coordinate
    system representation.
    """
    def __init__(self, coordinate='Cartesian',
                 r = sym.Symbol('x'),
                 t = sym.Symbol('y'),
                 p = sym.Symbol('z'),
                 para1=sym.Symbol('R0'),
                 para2=sym.Symbol('epsilon'),
                 new_contra_metric = sym.Matrix(
                                           [
                                               [1, 0, 0],
                                               [0, 1, 0],
                                               [0, 0, 1]
                                           ]
                                       )
                 ):
        """
        :param coordinate: (str) name of the coordinate system
                         inherent coordinate systems:
                         'Cartesian'
                         'Cylinder'
                         'Sphere'
                         'Toroidal'
        :param r: (sympy.Symbol) name of the coordinate, xi^1
        :param t: (sympy.Symbol) name of the coordinate, xi^2
        :param p: (sympy.Symbol) name of the coordinate, xi^3
        :param para1: (sympy.Symbol) name of one of parameters of the coordinate system,
                      only used for special coordinate systems such as the 'Toroidal'
        :param para2: (sympy.Symbol) name of one of parameters of the coordinate system,
                      only used for special coordinate systems such as the 'Toroidal'
        :param new_contra_metric: (sympy.Matrix) user defined metrix
        """
            pass









