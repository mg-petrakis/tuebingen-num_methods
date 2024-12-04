# -*- coding: utf-8 -*-

import numpy as np
"""
Problem 2: calculate integral let's say from  0 to 2π
"""
def f(x):
    if x == 0:
        return 1
    else:
        temp = (2 ** x)*np.sin(x)/x
        return temp

#first we make codes for each method

def simpsons_third(funct, init_x, final_x, n: int):
    # where we use a function as an argument, the integration limits and the
    # number of steps to the final point: h = (b-a)/n
    # simpson's 1/3 rule uses an even number of steps
    xpts = np.linspace(init_x, final_x, num = (n+1))
    h = (final_x - init_x)/n
    temp = 0
    # here we choose to ignore the sub interval if n odd
    #so we need to choose large n
    for i in range(int(n/2)):
        temp += (h/3)*(funct(xpts[2*i]) + 4*funct(xpts[2*i + 1]) + funct(xpts[2*i + 2]))
    return temp

def romberg_simpson_third(funct, init_x, final_x, n: int):
    xpts = np.linspace(init_x, final_x, num = (n+1))
    #for the improved method n should be a multiple of four
    if n%4 != 0:
        print("Please use as n a multiple of 4.")
        return 1
    else:
        h = 2*(final_x - init_x)/n
        temp = 0
        for i in range(int(n/4)):
            temp += (h/45)*(7*funct(xpts[4*i]) + 7*funct(xpts[4*i + 4]) + 32*funct(xpts[4*i + 1]) + 32*funct(xpts[4*i + 3]) + 12*funct(xpts[4*i + 2]))
        return temp

def four_point_gauss_legendre(funct, init_x, final_x):
    xpoints = [-0.861136116, -0.3394810436, 0.3394810436, 0.8611363116]
    Acoeffs = [0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451]
    #see pdf for the substitution
    mu = (init_x + final_x)/2
    lamda = (final_x - init_x)/2
    temp = 0
    for i in range(4):
        temp += lamda*Acoeffs[i]*funct(mu + lamda*xpoints[i])
    return temp

"""
#Test: should give ~2
print(four_point_gauss_legendre(np.sin, 0, np.pi))
"""

#So, for the given integral:
print("For the integral of (2^x)sinc(x) from 0 to 2π, we calculate:\n")
prob2_simps = simpsons_third(f, 0, 2*np.pi,4)
print("With the 1/3 Simpson's rule, 4 points: ", prob2_simps)
prob2_romb = romberg_simpson_third(f, 0, 2*np.pi, 4)
print("With the Romberg - 1/3 Simpson's rule, 4 points: ", prob2_romb)
prob2_gauss = four_point_gauss_legendre(f, 0, 2*np.pi)
print("With the 4-point Gauss-Legendre method: ", prob2_gauss)
print("\nThe most accurate result, using 1028 steps with the Romberg-Simpson rule:\n")
print(romberg_simpson_third(f, 0, 2*np.pi, 1028))

"""
Finally, Problem 3: 2D integrals
"""
print("\nProblem 3:")
# 1st integral: it's composed of seperable polynomials, so has higher derivatives
# zero and should be exact
def f1x(x):
    return x
def f1y(x):
    return x**2
res = romberg_simpson_third(f1x, 0 , 2, 4)*romberg_simpson_third(f1y, 0 , 1, 4)
print("\nFor the first integral:", res)
print("while true value: 2/3.")


#2nd integral: mixed limits of integration
# Let's make a variation:

def twod_simpsons_third(func2d, a1, b1, a2, b2, n: int):
    # where (a1,b1) and (a2,b2) respective limits of integration
    x1_pts = np.linspace(a1, b1, num = (n+1))
    h = (b1 - a1)/n
    temp = 0
    def funct(x2):
        return simpsons_third(lambda x1: func2d(x1,x2), a2(x2), b2(x2), n)
    for i in range(int(n/2)):
        temp += (h/3)*(funct(x1_pts[2*i]) + 4*funct(x1_pts[2*i + 1]) + funct(x1_pts[2*i + 2]))
    return temp

def f(x,y):
    return (x*y*y)

def a2(x):
    return 2*x

def b2(x):
    return 2

res2 = twod_simpsons_third(f, 0, 1, a2, b2, 8)
print("\nFor the second integral:", res2)
print("while true value: 4/15 = 0.26666666666.")

#3rd integral: similarly but with different order

def a3(x):
    return 0

def b3(x):
    return (x/2)

def g(x,y):
    return f(y,x)

res3 = twod_simpsons_third(g, 0, 2, a3, b3, 32)
print("\nFor the third integral:", res3)
print("while true value: 4/15 = 0.26666666666.")
#as verified in Mathematica