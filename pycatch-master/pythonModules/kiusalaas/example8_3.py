## example8_3
import numpy
import run_kut5
import linInterp
import printSoln

def initCond(u):  # Initial values of [y,y',y"];
                  # use 'u' if unknown
    return numpy.array([0.0, 0.0, u])  

def r(u):  # Boundary condition residual--see Eq. (8.3)
    X,Y = run_kut5.integrate(F,xStart,initCond(u),xStop,h)
    y = Y[len(Y) - 1]
    r = y[0] - 2.0  
    return r

def F(x,y):  # First-order differential equations                     
    F = numpy.zeros((3),dtype=float)
    F[0] = y[1]
    F[1] = y[2]
    F[2] = 2.0*y[2] + 6.0*x*y[0]
    return F

xStart = 5.0        # Start of integration
xStop = 0.0         # End of integration
u1 = 1.0            # 1st trial value of unknown init. cond.
u2 = 2.0            # 2nd trial value of unknown init. cond.
h = -0.1            # initial step size
freq = 2            # printout frequency
u = linInterp.linInterp(r,u1,u2)
X,Y = run_kut5.integrate(F,xStart,initCond(u),xStop,h)
printSoln.printSoln(X,Y,freq)
raw_input("\nPress return to exit")
