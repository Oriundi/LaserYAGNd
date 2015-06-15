#!/usr/bin/env python3

'''
Documentation, License etc.
Xiao-Bass speed equations for
YAG:Nd laser
@package laserYAGNd
'''

import matplotlib.pyplot as plt
from scipy.integrate import ode, simps
from scipy.signal import argrelmax
from config import *
from optimization.optim_R2 import optimizeR2
# import logging
from optimization.optim_lg import optimizelg
from optimization.optim_la import optimizela
from optimization.optim_T0 import optimizeT0
import time
# import multiprocessing
# import os
# global R

optimize_events = {
    "Optimize_R2": optimizeR2,
    "Optimize_lg": optimizelg,
    "Optimize_la": optimizela,
    "Optimize_T0": optimizeT0
}

interval_events = {
    "Optimize_R2": numpy.linspace(0.7, 1, num=5, retstep=True),
    "Optimize_lg": numpy.linspace(500, 1500, num=5, retstep=True),
    "Optimize_la": numpy.linspace(20, 150, num=5, retstep=True),
    "Optimize_T0": numpy.linspace(0.7, 1, num=5, retstep=True)
}


def equation(t, x):
    # Xiao-Bass speed equations
    return numpy.array([R - Cg * x[0] * x[2] - x[0]/tau_g,
                        Ca * x[1] * x[2] + (Na_total - x[1])/tau_a,
                        (x[0]*g - x[1]*a1 - (Na_total-x[1])*a2 - 2*gamma)
                        * x[2]/tau_r + (x[0]+Ng_total)*Cg_eps])

r = ode(equation).set_integrator('dop853', atol=1e-12, rtol=1e-12)
r.set_initial_value(x0, time_start)


def laserYAGNd(x_out, tau):
    # ii = 0
    # global R
    while r.successful() and r.t < time_end:
        r.integrate(r.t+dt)
        x_out = numpy.vstack([x_out, numpy.abs(r.y)])
        tau = numpy.vstack([tau, r.t])
        # print(x_out[ii, 0])
        # R11 = Pi * (1 - math.exp(-2 * sigma_p * x_out[ii, 0] * lg))  # um^-3 * us^-1
        # R22 = A * lg * h * nu_p
        # R = R11 / R22 * 1e-3
        # ii += 1
    nd = x_out[:, 0]
    na = x_out[:, 1]
    q = x_out[:, 2]
    p = h*nu_l/tau_r * math.log(1/R2) * q * 1e6
    tau[0] = 0
    return nd, na, q, p, tau


def main():
    # pool_size = multiprocessing.cpu_count()
    # os.system('taskset -cp 0-%d %s' % (pool_size, os.getpid()))
    # print("Pool size:", pool_size)
    # pool = multiprocessing.Pool(processes=pool_size)
    # print("Ng total %f " % Ng_total)
    # p # rint("Na total %f" % Na_total)

    optimize = False
    # optimize = "Optimize_R2"
    # optimize = "Optimize_lg"
    # optimize = "Optimize_la"
    # optimize = "Optimize_T0"
    print('Optimize event - %s' % optimize)

    # print("R11 = %f" % R11, "R22 = %f" % R22, "R = %f" % R, "R_ = %f" % R_)
    # print("Ng_total = %f" % Ng_total, "Na_total = %f" % Na_total)

    if optimize is False:
        _time = time.time()
        Nd, Na, q, P, tau = laserYAGNd(x_out_initial, tau_initial)
        print("Elapsed time %f" % (time.time() - _time))
        freq = 0
        imp_width = 0
        avg_P = 0
        avg_imp_energy = 0

        #print("local max %s" % local_max)
        # print((numpy.where(P[local_max] > 1)))
        # print(numpy.size(numpy.where(P[local_max] > 1)))

        # P = h*nu_p/tau_r * math.log(1/R2) * q * 1e6

        local_max = argrelmax(P, order=100)
        for jj in local_max[0]:
            if P[jj] < 2:
                local_max = numpy.delete(local_max, numpy.where(local_max == jj))
        try:
            local_max = numpy.delete(local_max, 0)
            local_max = numpy.delete(local_max, len(local_max)-1)
        except IndexError:
            print("No maximum over 2 W")
            pass
        # print("local max %s" % local_max)

        tau_m = tau[local_max]
        # Find laser frequency
        for jj in range(0, len(tau_m)-1):
            freq += 1 / (tau_m[jj+1] - tau_m[jj]) * 1e3
        freq /= (numpy.size(local_max) - 1)
        print("Laser frequency = %2.4f kHz" % freq)
        # Find average impulse width
        width_cut = 2
        for jj in local_max:
            more = jj
            less = jj
            while P[more] > P[jj]/width_cut:
                more += 1
            while P[less] > P[jj]/width_cut:
                less -= 1
            imp_width += ((tau[more] + tau[more-1])/2 - (tau[less]+tau[less+1])/2)
            # print("Local max %f" % tau[jj], "Left halfcut %f" % tau[less], "Right halfcut %f" % tau[more])
        imp_width /= numpy.size(local_max) - 1
        print("Impulse average width %2.4f ns" % (imp_width * 1e3))
        # Find average impulse energy
        energy_cut = 2
        for jj in local_max:
            more = jj
            less = jj
            while P[more] > P[jj]/energy_cut:
                more += 1
            while P[less] > P[jj]/energy_cut:
                less -= 1
            x = P[less:more]
            imp_energy = simps(x, dx=dt)
            avg_imp_energy += imp_energy
            # print("Impulse energy %2.4f" % imp_energy)
        avg_imp_energy /= numpy.size(local_max) - 1
        print("Impulse average energy %2.4f nJ" % avg_imp_energy)
        # Find average impulse max power
        avg_P += numpy.sum(P[local_max])
        avg_P /= numpy.size(local_max)
        print("Impulse average power %2.4f W" % avg_P)

        # print('tau0 = %f' % tau[0])
        # print('tau_last = %f' % tau[len(tau)-1])

        fig = plt.figure(1)
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)

        ax1.set_xlabel('Time, us')
        ax1.set_ylabel('Concentration, Nd - blue, Cr - red ')

        ax2.set_xlabel('Time, us')
        ax2.set_ylabel('Power, W')

        ax1.plot(tau, Nd * 1e12, 'b')
        ax1.hold()
        plt_ax1 = ax1.plot(tau, Na * 1e12, color='r')
        ax1.hold()

        ax2.plot(tau, P, color='b')
        ax2.hold()

        # plt.savefig('out.jpg', format='jpg', dpi=600)
        plt.show()
    else:
        interval = interval_events[optimize]
        out = optimize_events[optimize](interval)


if __name__ == '__main__':
    main()
