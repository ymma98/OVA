"""
Orthogonal curvilinear coordinates Vector Analysis (OVA)
ymma98@qq.com
"""

import sympy as sym


def get_co_metric(x_expr, y_expr, z_expr,
                  r, t, p):
    """
    assuming the new coordinate (r, t, p)
    g_{ij} = x.diff(xi^i) * x.diff(xi^j)
    :param x_expr (sympy expr): expression x = x_expr(r, t, p)
    :param y_expr (sympy expr): expression y = y_expr(r, t, p)
    :param z_expr (sympy expr): expression z = z_expr(r, t, p)
    :param r (sympy symbol): 1st coordinate symbol, like 'r'
    :param t (sympy symbol): 2nd coordinate symbol, like 'theta'
    :param p (sympy symbol): 3rd coordinate symbol, like 'z'
    """
    coord = [r, t, p]
    co_metric = sym.Matrix(3, 3, lambda i, j: sym.simplify(
                           x_expr.diff(coord[i]) * x_expr.diff(coord[j]) + \
                           y_expr.diff(coord[i]) * y_expr.diff(coord[j]) + \
                           z_expr.diff(coord[i]) * z_expr.diff(coord[j]) ) )
    return co_metric

def get_contra_metric(r_expr, t_expr, p_expr,
                      r, t, p,
                      x=sym.Symbol('x'),
                      y=sym.Symbol('y'),
                      z=sym.Symbol('z')):
    """
    assuming the new coordinate (r, t, p)
    g^{ij} = xi^i.diff(x^i) * xi^j.diff(x^j)
    :param r_expr (sympy expr): expression r = r_expr(x, y, z)
    :param t_expr (sympy expr): expression t = t_expr(x, y, z)
    :param p_expr (sympy expr): expression p = p_expr(x, y, z)
    :param r (sympy symbol): 1st coordinate symbol, like 'r'
    :param t (sympy symbol): 2nd coordinate symbol, like 'theta'
    :param p (sympy symbol): 3rd coordinate symbol, like 'z'
    :param x (sympy symbol): 1st coordinate symbol in Cartesian, like 'x'
    :param y (sympy symbol): 2nd coordinate symbol in Cartesian, like 'y'
    :param z (sympy symbol): 3rd coordinate symbol in Cartesian, like 'z'
    """
    xi_expr = [r_expr, t_expr, p_expr]
    contra_metric = sym.Matrix(3, 3, lambda i, j: sym.simplify(
                           xi_expr[i].diff(x) * xi_expr[j].diff(x) + \
                           xi_expr[i].diff(y) * xi_expr[j].diff(y) + \
                           xi_expr[i].diff(z) * xi_expr[j].diff(z)
                           ) )
    contra_metric = contra_metric.subs(r_expr, r)
    contra_metric = contra_metric.subs(t_expr, t)
    contra_metric = contra_metric.subs(p_expr, p)
    return contra_metric




class CoordinateSystem:
    """
    defines a orthogonal curvilinear coordinate
    system representation.
    """
    def __init__(self, coordinate='Cartesian',
                 r = sym.Symbol('x'),
                 t = sym.Symbol('y'),
                 p = sym.Symbol('z'),
                 param1=sym.Symbol('R0'),
                 param2=sym.Symbol('epsilon'),
                 new_co_metric = None
                 ):
        """
        :param coordinate: (str) name of the coordinate system
                         inherent coordinate systems:
                         'Cartesian'
                         'Cylindrical'
                         'Spherical'
                         'Toroidal' (not completed)
        :param r: (sympy.Symbol) name of the coordinate, xi^1
        :param t: (sympy.Symbol) name of the coordinate, xi^2
        :param p: (sympy.Symbol) name of the coordinate, xi^3
        :param param1: (sympy.Symbol) name of one of parameters of the coordinate system,
                      only used for special coordinate systems such as the toroidal system
        :param param2: (sympy.Symbol) name of one of parameters of the coordinate system,
                      only used for special coordinate systems such as the toroidal system
        :param new_co_metric: (sympy.Matrix) user defined contra-variant metrix, 3x3
        """
        self.coord = coordinate
        self.r = r
        self.t = t
        self.p = p
        self.param1 = param1
        self.param2 = param2

        if not new_co_metric:  # use default coordinate system
            self.use_default_contra_metric()
        else:  # user defined coordinate system
            self.co_metric = new_co_metric
            self.contra_metric = new_co_metric.inv()

        self.g = self.co_metric.det()
        self.V = sym.sqrt(self.g)
        self.J = self.V
        self.h1 = sym.sqrt(self.co_metric[0, 0])
        self.h2 = sym.sqrt(self.co_metric[1, 1])
        self.h3 = sym.sqrt(self.co_metric[2, 2])
        self.h = [self.h1, self.h2, self.h3]

    def dot(self, A, B):
        """
        :param A: (list) 1x3 list (vector) or 3x3 list (tensor)
        :param B: (list) 1x3 list (vector) or 3x3 list (tensor)
        """
        case_a = 'v'
        case_b = 'v'
        if isinstance(A, list) and isinstance(A[0], list):
            # A is a tensor
            case_a = 't'
        if isinstance(A, sym.Matrix):
            if A.shape == (3, 3):
                case_a = 't'
            else:
                raise ValueError("sym.Matrix shape is not (3, 3)")

        if isinstance(B, list) and isinstance(B[0], list):
            # B is a tensor
            case_b = 't'
        if isinstance(B, sym.Matrix):
            if B.shape == (3, 3):
                case_b = 't'
            else:
                raise ValueError("sym.Matrix shape is not (3, 3)")

        case = case_a + case_b  # 'vv' means vector dot vector
                     # 'vt' means vector dot tensor
                     # 'tv' means tensor dot vector

        if case == 'vv':
            res = sym.S(0)
            for i in [0, 1, 2]:
                res += A[i] * B[i]
            # simplify res
            res = sym.simplify(res)
        elif case == 'vt':
            # a dot T = Ai * Tij e_j
            res = [sym.S(0), sym.S(0), sym.S(0)]
            for i in [0, 1, 2]:
                for j in [0, 1, 2]:
                    if not isinstance(B, sym.Matrix):
                        res[j] += A[i] * B[i][j]
                        res[j] = sym.simplify(res[j])
                    else:
                        res[j] += A[i] * B[i, j]
                        res[j] = sym.simplify(res[j])

        elif case == 'tv':
            # T dot b = bj * Tij e_i
            res = [sym.S(0), sym.S(0), sym.S(0)]
            for i in [0, 1, 2]:
                for j in [0, 1, 2]:
                    if not isinstance(A, sym.Matrix):
                        res[i] += A[i][j] * B[j]
                        res[i] = sym.simplify(res[i])
                    else:
                        res[i] += A[i, j] * B[j]
                        res[i] = sym.simplify(res[i])
        else:
            raise ValueError("wrong input.Shape!")

        return res

    def cross(self, A, B):
        """
        :param A: (list) 1x3 list (vector)
        :param B: (list) 1x3 list (vector)
        """
        res = [sym.S(0), sym.S(0), sym.S(0)]
        for i in [0, 1, 2]:
            for j in [0, 1, 2]:
                for k in [0, 1, 2]:
                    res[k] += self.levi_civita(i,j,k) * \
                              A[i] * B[j]
        for i in range(len(res)):
            res[i] = sym.simplify(res[i])
        return res

    def grad(self, u):
        """
        u could be a scalar function or a vector
        :param u: (sympy.Function or expression or list/vector)
        """
        xi_list = [self.r, self.t, self.p]
        h_list = [self.h1, self.h2, self.h3]
        if isinstance(u, list):
            # u is a vector
            res = sym.Matrix([
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0]] )

            for i in [0, 1, 2]:
                for j in [0, 1, 2]:
                    res[i, j] += u[j].diff(xi_list[i]) / h_list[i]
                    res[i, i] += u[j] / (h_list[i] * h_list[j]) \
                                   * h_list[i].diff(xi_list[j])
            for i in [0, 1, 2]:
                for k in [0, 1, 2]:
                    res[i, k] -= u[i] / (h_list[i] * h_list[k]) * \
                                 h_list[i].diff(xi_list[k])
            # simplify
            for i in [0, 1, 2]:
                for j in [0, 1, 2]:
                    res[i, j] = sym.simplify(res[i, j])
        else:
            # u is a scalar function or expression
            res = [sym.S(0), sym.S(0), sym.S(0)]
            for i in [0, 1, 2]:
                res[i] = sym.simplify( u.diff(xi_list[i]) / h_list[i] )
        return res


    def curl(self, u):
        """
        :param u: (list) 3x1 vector
        """
        xi_list = [self.r, self.t, self.p]
        h_list = [self.h1, self.h2, self.h3]
        res = [sym.S(0), sym.S(0), sym.S(0)]
        for k in [0, 1, 2]:
            for i in [0, 1, 2]:
                for j in [0, 1, 2]:
                    temp_expr = h_list[i] * u[i]
                    res[k] += self.levi_civita(k, j, i) * \
                              h_list[k] / self.V * \
                              temp_expr.diff(xi_list[j])
        for i in [0, 1, 2]:
            res[i] = sym.simplify(res[i])

        return res

    def div(self, u):
        """
        u could be a scalar function or a vector
        :param u: (list) 1x3 list (vector) or 3x3 list (tensor)
        """
        xi_list = [self.r, self.t, self.p]
        h_list = [self.h1, self.h2, self.h3]
        case_u = 'v'  # default: u is a 1x3 vector
        if isinstance(u, list) and isinstance(u[0], list):
            # u is a tensor in list format
            case_u = 't'
        if isinstance(u, sym.Matrix):
            # B is a tensor in sym.Matrix format
            if u.shape == (3, 3):
                case_u = 't'
            else:
                raise ValueError("sym.Matrix shape is not (3, 3)")

        if case_u == 'v':
            res = sym.S(0)
            for i in [0, 1, 2]:
                temp_expr = self.V * u[i] / h_list[i]
                res += sym.S(1) / self.V * temp_expr.diff(xi_list[i])
        elif case_u == 't':
            res = [sym.S(0), sym.S(0), sym.S(0)]
            if not isinstance(u, sym.Matrix):
                for i in [0, 1, 2]:
                    for k in [0, 1, 2]:
                        res[k] += sym.S(1) / h_list[i] * u[i][k].diff(xi_list[i])
                        temp_expr = self.V / h_list[i]
                        res[k] += u[i][k] / self.V * temp_expr.diff(xi_list[i])
                for i in [0, 1, 2]:
                    for k in [0, 1, 2]:
                        res[i] += u[i][k] / (h_list[i] * h_list[k]) \
                                  * h_list[i].diff(xi_list[k])
                for i in [0, 1, 2]:
                    for j in [0, 1, 2]:
                        res[j] -= u[i][i] / (h_list[i] * h_list[j]) \
                                  * h_list[i].diff(xi_list[j])
            else:
                for i in [0, 1, 2]:
                    for k in [0, 1, 2]:
                        res[k] += sym.S(1) / h_list[i] * u[i, k].diff(xi_list[i])
                        temp_expr = self.V / h_list[i]
                        res[k] += u[i, k] / self.V * temp_expr.diff(xi_list[i])
                for i in [0, 1, 2]:
                    for k in [0, 1, 2]:
                        res[i] += u[i, k] / (h_list[i] * h_list[k]) \
                                  * h_list[i].diff(xi_list[k])
                for i in [0, 1, 2]:
                    for j in [0, 1, 2]:
                        res[j] -= u[i, i] / (h_list[i] * h_list[j]) \
                                  * h_list[i].diff(xi_list[j])
        else:
            raise ValueError(r"wrong case_u value, neither 'v' nor 't'")
        return res

    def Laplacian(self, u):
        """
        :param u: (sym.Function or expression)
        """
        res = sym.S(0)
        xi_list = [self.r, self.t, self.p]
        h_list = [self.h1, self.h2, self.h3]
        for i in [0, 1, 2]:
            temp_expr = self.V / (h_list[i]**2) * u.diff(xi_list[i])
            res += sym.S(1) / self.V * temp_expr.diff(xi_list[i])
        return res





    def levi_civita(self, i, j, k):
        """
        Levi-Civita symbol, i.e. \epsilon_{ijk} or \epsilon&{ijk}
        :param i: (int), 0 - 2
        :param j: (int)  0 - 2
        :param k: (int)  0 - 2
        """
        res = sym.S(0)
        if (i==j or i==k or j==k):
            res = sym.S(0)
        elif (
             ([i,j,k] == [0,1,2]) or  # (1,2,3)
             ([i,j,k] == [1,2,0]) or  # (2,3,1)
             ([i,j,k] == [2,0,1])     # (3,1,2)
        ):
          res = sym.S(1)
        elif (
             [i,j,k] == [0,2,1] or  # (1,3,2)
             [i,j,k] == [1,0,2] or  # (2,1,3)
             [i,j,k] == [2,1,0]     # (3,2,1)
        ):
          res = sym.S(-1)
        else:
            raise ValueError("wrong input for Levi-Civita symbol")

        return res





    def use_default_contra_metric(self):
        if (self.coord == 'Cartesian'):
            self.contra_metric = sym.Matrix(
                [
                    [1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1]
                ] )
        elif (self.coord == 'Cylindrical'):
            # add assumption: r >= 0
            self.contra_metric = sym.Matrix(
                [
                    [1, 0,             0],
                    [0, 1 / self.r**2, 0],
                    [0, 0,             1]
                ] )
        elif (self.coord == 'Spherical'):
            self.contra_metric = sym.Matrix(
                [
                    [1, 0,             0                                   ],
                    [0, 1 / self.r**2, 0                                   ],
                    [0, 0,             1 / (self.r**2 * sym.sin(self.t)**2)]
                ] )
        else:
            raise ValueError("coordinate {0} does not exists!".format(self.coord))
        self.co_metric = self.contra_metric.inv()









