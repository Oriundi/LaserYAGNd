#!/usr/bin/env python

'''
Documentation, License etc.
Xiao-Bass speed equations for
YAG:Nd laser
@package laserYAGNd
'''

import numpy
import matplotlib.pyplot as plt
from scipy.integrate import ode
from config import *


def equation(t, x):
    # Xiao-Bass speed equations
    return numpy.array([R - Cg * x[0] * x[2] - x[0]/tau_g,
                        Ca * x[1] * x[2] + (Na_total - x[1])/tau_a,
                        (x[0]*g - x[1]*a1 - (Na_total-x[1])*a2 - 2*gamma)
                        * x[2]/tau_r + (x[0]+Ng_total)*Cg_eps])

r = ode(equation).set_integrator('dopri5', atol=1e-12, rtol=1e-12)
r.set_initial_value(x0, t0)


def laserYAGNd(x_out, tau):
    while r.successful() and r.t < t1:
        r.integrate(r.t+dt)
        x_out = numpy.vstack([x_out, numpy.abs(r.y)])
        tau = numpy.vstack([tau, r.t])
    nd = x_out[:, 0]
    na = x_out[:, 1]
    q = x_out[:, 2]
    p = h*nu_p/tau_r * math.log(1/R2) * q * 1e6
    return nd, na, q, p, tau


def main():
    Nd, Na, q, P, tau = laserYAGNd(x_out_initial, tau_initial)
    plt.figure(1)
    plt.subplot(311)
    plt.plot(tau, Nd)
    plt.subplot(312)
    plt.plot(tau, Na)
    plt.subplot(313)
    plt.plot(tau, P)
    plt.show()
    # sys.exit(app.exec_())

if __name__ == '__main__':
    main()
