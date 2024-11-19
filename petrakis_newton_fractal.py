# -*- coding: utf-8 -*-
# Author: Marios Gavriil Petrakis
#we import the complex numbers module cmath, where
#the imaginary unit is used as "1j"
import matplotlib.pyplot as plt
import math
import cmath

#we define f(z) and f'(z) regularly, as it is an algebraic function
def f(z: complex):
    return (z**3 -1)
#we define the complex derivative for the Newton method
def gradf(z: complex):
    return (3*z**2)

#z0 initial guess, epsilon tolerance, M max steps

def solve_cnewton(z0: complex, epsilon: float, M: int):
# we initialize the iteration number to 1
    n_iter = 1
# z_help is used to check convergence according to the tolerance
    z_help = z0 + 0.5
# a list to store the consecutive differences Δz_n and test convergence
    conv_test = []
    while abs(z_help-z0) > epsilon:
        z_help = z0
# condition to terminate the loop if the iterations exceed the limit M        
        if n_iter > M:
            break
# break the loop if we happen to land on a root
        elif f(z0) == 0 +0j:
            break
# proper Newton method
        elif gradf(z0) != 0 + 0j:
            z0 -=f(z0)/gradf(z0)
            conv_test.append(abs(z0-z_help))
# escape stationary points (z=0)
        else:
            z0 += 0.0001
        n_iter += 1
# to output mulitple different data types I'll use a dictionary
    temp = {"root": z0, "iterations": n_iter, "quality": conv_test}
    return temp

# fractal plotter

# exact roots
root1 = 1 +0j
root2 = -0.5 + cmath.sin((2/3)*cmath.pi)*1j
root3 = -0.5 - cmath.sin((2/3)*cmath.pi)*1j

# data function
def gridsol(N: int, z1: complex, z2: complex, eps: float, M: int):
# store the initial point - result data in different lists
    startpoints = []
    rootfound = []
    iterlog = []
    # positive real part: moving rightwards
    a = (z2-z1).real
    # negative imaginary part: moving downwards
    b = (z2-z1).imag
    for k in range(N+1):
        for l in range(N+1):
            z0 = z1 + k*(a/N) + l*(b/N)*1j
            # record starting points
            startpoints.append(z0)
            newt = solve_cnewton(z0, eps, M)
            # record roots
            rootfound.append(newt["root"])
            # record iterations number in log10 corresponding to the roots
            iterlog.append(math.log10(newt["iterations"]))
    temp = {"start": startpoints, "roots": rootfound, "log_itns": iterlog}
    return temp

# subprogramme take the real, imaginary parts of list values
def listimag(lista):
    return [z.imag for z in lista]

def listreal(lista):
    return [z.real for z in lista]

# plotting as a heat map

def fractalplotter(N: int, z1: complex, z2: complex, eps: float, M: int):
    grid = gridsol(N, z1, z2, eps, M)
    zpoints = grid["start"]
# we use the imaginary part because it is distinct for the 3 roots
    values = listimag(grid["roots"])
    #scatter plot and mapping color coded points on the complex plane
    plt.scatter(x=listreal(zpoints), y=listimag(zpoints), s=3, c=values, cmap='viridis')
    #normalize color by the min and max imag parts of the roots
    plt.clim(-math.sin((2/3) * math.pi), math.sin((2/3) * math.pi))
    plt.colorbar(label="Imaginary part of converged root")
    # same scale for both axes so the image is not distorted
    plt.axis('scaled')
    plt.title("Basins of attraction with Newton method")
    plt.xlabel("Re(z)")
    plt.ylabel("Im(z)")
    plt.xlim((z1.real, z2.real))
    plt.ylim((z2.imag, z1.imag))
    plt.savefig("3rd_zoom.png", dpi=1200)

def iterplotter(N: int, z1: complex, z2: complex, eps: float, M: int):
    grid = gridsol(N, z1, z2, eps, M)
    zpoints = grid["start"]
# we use the imaginary part because it is distinct for the 3 roots
    values = grid["log_itns"]
    plt.scatter(x=listreal(zpoints), y=listimag(zpoints), s=0.01, c=values, cmap='plasma')
    plt.clim(0, math.log10(M))
    plt.colorbar(label="$\log_{10}(n_{iter})$")
    plt.axis('scaled')
    plt.title("Iterations for convergence")
    plt.xlabel("Re(z)")
    plt.ylabel("Im(z)")
    plt.xlim((z1.real, z2.real))
    plt.ylim((z2.imag, z1.imag))
    plt.savefig("whole_t.png", dpi=1200)
    
#fractalplotter(1000,-0.80 + 0.01j, -0.78 - 0.01j, 1e-6, 200)
iterplotter(1000, -2 + 2j, 2-2j, 1e-6, 100)
#print(gridsol(150, -2 + 2j, 2-2j, 1e-6, 100)["log_itns"])