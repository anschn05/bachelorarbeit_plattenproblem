


# import numpy as np
# import sympy as sp
# x = sp.Symbol('x')
# y = sp.Symbol('y')

# E = 70e9
# t = 0.02
# nu = 0.23

# def f(x,y):
#     return x**2 * (1-x)**2 * y**2 * (1-y)**2
# def f_x(x,y):
#     return 2*x*(1-x)**2*y**2*(1-y)**2 - 2*x**2*(1-x)*y**2*(1-y)**2
# def f_y(x,y):
#     return 2*y*(1-y)**2*x**2*(1-x)**2 - 2*y**2*(1-y)*x**2*(1-x)**2

# def f_xy(x,y):
#     return (2*x*(1-x)**2 - 2*x**2*(1-x))*(2*y*(1-y)**2 - 2*y**2*(1-y))

# def f_xx(x,y):
#     return (y**2*(1-y)**2)(2*(1-x)**2-8*x*(1-x)-2*x**2)

# def f_yy(x,y):
#     return x**2*(1-x)**2(2*(1-y)**2-8*y*(1-y)+2*y**2)

# def D():
#     return t**3/12*E/(1+nu)

# def C_dach(A):
#     return D()*(A + nu/(1-nu)*np.trace(A)*np.eye(2))



# if __name__=="__main__":
#     nabla2_w = np.matrix([[f_xx(x,y),f_xy(x,y)],[f_xy(x,y),f_yy(x,y)]])
#     M = C_dach(nabla2_w)

