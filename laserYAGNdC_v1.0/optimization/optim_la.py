from scipy.signal import argrelmax, find_peaks_cwt
from scipy.integrate import ode
from config import *
from time import time
import matplotlib.pyplot as plt
import logging
from scipy.integrate import simps

numpy.set_printoptions(precision=8)


def solve_equation(x00, t00, t11, tau_r0, Cg0, Ca0, a10, a20, Na_total0, Cg_eps0, x_out, tau):
    def equation(t, x):
        # Xiao-Bass speed equations
        return numpy.array([R - Cg0 * x[0] * x[2] - x[0]/tau_g,
                           Ca0 * x[1] * x[2] + (Na_total0 - x[1])/tau_a,
                           (x[0]*g - x[1]*a10 - (Na_total-x[1])*a20 - 2*gamma)
                           * x[2]/tau_r0 + (x[0]+Ng_total)*Cg_eps0])

    r = ode(equation).set_integrator('dop853', atol=1e-12, rtol=1e-12)
    r.set_initial_value(x00, t00)

    while r.successful() and r.t < t11:
        r.integrate(r.t+dt)
        x_out = numpy.vstack([x_out, numpy.abs(r.y)])
        tau = numpy.vstack([tau, r.t])

    q = x_out[:, 2]
    P = h*nu_l/tau_r * math.log(1/R2) * q * 1e6
    local_max = argrelmax(P)
    print("Found %i local maxima" % numpy.size(numpy.where(P[local_max] > 1)))
    if numpy.size(numpy.where(P[local_max] > 1)) < 5:
        return "Not success", t11
    else:
        return "Success", x_out, tau


def optimizela(interval):
    _time = time()
    freq = numpy.zeros([numpy.size(interval[0])])
    # print(numpy.size(freq))
    imp_width = numpy.zeros([numpy.size(interval[0])])
    avg_P = numpy.zeros([numpy.size(interval[0])])
    avg_imp_energy = numpy.zeros([numpy.size(interval[0])])
    max_calc_time = 500
    for ii in range(0, numpy.size(interval[0])):
        tau = time_start
        la = interval[0][ii]   # reflectivity of output mirror
        lr = n * (lg + la)
        Veff = lr / la * Vg         # um^3, effective mode volume
        tau_r = 2 * lr / c0
        Cg = sigma_g * c0 / Veff    # 1/us
        Ca = -1 * sigma_a1 * c0 / Veff  # 1/us
        a1 = 2 * sigma_a1 * la  # um^3
        a2 = 2 * sigma_a2 * la  # um^3
        Na_total = math.log(T0) / (-1 * sigma_a1 * la)
        Cg_eps = eps * c0 * sigma_g * lg / lr  # um^3/us.
        x0 = numpy.array([Ng_total, Na_total, 0])
        x_out = x0
        # print(R2)
        # gamma_2 = -math.log(R2)  # losses on output mirror
        # gamma = gamma_i + 0.5 * (gamma_1 + gamma_2)  # losses
        result = [0]
        calc_time = 50
        while result[0] is not "Success" and calc_time <= max_calc_time:
            result = solve_equation(x0, time_start, calc_time, tau_r, Cg, Ca, a1, a2, Na_total, Cg_eps, x_out, tau)
            if result[0] is not "Success" and calc_time <= max_calc_time:
                logging.warning("No enough local maxima. Increase calculation time to %i us" % (calc_time+50))
            calc_time += 50

        try:
            # print(result[0])
            x_out = result[1]
            tau = result[2]

            nd = x_out[:, 0]
            na = x_out[:, 1]
            q = x_out[:, 2]
            P = h*nu_l/tau_r * math.log(1/R2) * q * 1e6
        except IndexError:
            logging.warning("There is no impulses in first %i us. Go to next value." % calc_time)

        try:
            local_max = argrelmax(P, order=100)
            for jj in local_max[0]:
                if P[jj] < 5:
                    local_max = numpy.delete(local_max, numpy.where(local_max == jj))
            local_max = numpy.delete(local_max, 0)
            local_max = numpy.delete(local_max, len(local_max)-1)
            # print("local max %s" % local_max)

            tau_m = tau[local_max]
            # Find laser frequency
            for jj in range(0, len(tau_m)-1):
                freq[ii] += 1 / (tau_m[jj+1] - tau_m[jj]) * 1e3
            freq[ii] /= (numpy.size(local_max) - 1)
            print("Laser frequency = %2.4f kHz" % freq[ii])
            # Find average impulse width
            width_cut = 2
            for jj in local_max:
                more = jj
                less = jj
                while P[more] > P[jj]/width_cut:
                    more += 1
                while P[less] > P[jj]/width_cut:
                    less -= 1
                imp_width[ii] += ((tau[more] + tau[more-1])/2 - (tau[less]+tau[less+1])/2)
                # print("Local max %f" % tau[jj], "Left halfcut %f" % tau[less], "Right halfcut %f" % tau[more])
            imp_width[ii] /= numpy.size(local_max) - 1
            print("Impulse average width %2.4f ns" % (imp_width[ii] * 1e3))
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
                avg_imp_energy[ii] += imp_energy
                # print("Impulse energy %2.4f" % imp_energy)
            avg_imp_energy[ii] /= numpy.size(local_max) - 1
            print("Impulse average energy %2.4f uJ" % avg_imp_energy[ii])
            # Find average impulse max power
            avg_P[ii] += numpy.sum(P[local_max])
            avg_P[ii] /= numpy.size(local_max)
            print("Impulse average max power %2.4f W" % avg_P[ii])

            # plt.figure(1)
            # plt.subplot(311)
            # plt.plot(tau, nd)
            # plt.subplot(312)
            # plt.plot(tau, na)
            # plt.subplot(313)
            # plt.plot(tau, P)
            # plt.show()
        except Exception:
            logging.warning("Cannot make calculations")

    plt.figure(2)
    plt.subplot(221)
    plt.xlabel("Absorber thickness")
    plt.ylabel("Laser frequency, kHz")
    plt.plot(interval[0], freq)
    plt.subplot(222)
    plt.xlabel("Absorber thickness")
    try:
        plt.ylabel("Impulse width on level P/%i" % width_cut)
    except Exception:
        print("width_cut is not set")
        pass
    plt.plot(interval[0], imp_width)
    plt.subplot(223)
    plt.xlabel("Absorber thickness")
    plt.ylabel("Average impulse power, W")
    plt.plot(interval[0], avg_P)
    plt.subplot(224)
    plt.xlabel("Absorber thickness")
    plt.ylabel("Average impulse energy, uJ")
    plt.plot(interval[0], avg_imp_energy)
    plt.show()

    print("Elapsed time on optimization %f" % (time() - _time))
    # return nd, na, q, p, tau
