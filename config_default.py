__author__ = 'oriundi'

import math
import numpy


tau_g = 230  # um, lifetime of high laser level (0.23ms)
tau_a = 3.5  # um, lifetime of excited level 3T2 of ions Cr(4+)(d) (3.5us)

# Parameters for the first equation
sigma_g = 3.5 * 1e-11  # um^2 efficient cross-section of laser transition
c0 = 3 * 1e8  # um^2/us, speed of light
rp = 300    # um, pumping beam radius
rl = 300    # um, laser beam radius
Pi = 0.2  # W, incident power
alpha_p = 5*1e-7  # 1/um, absorption coefficient on the wavelength
h = 6.626069 * 1e-34 * 1e6  # J*us, Plank constant
wavelength = 0.808  # um, radiation excitation wavelength
n = 1.8169  # YAG:Nd refractive index at 1.064 um
nu_p = c0 / wavelength  # 1/us, radiation excitation frequency

# Parameters for the second equation
sigma_a1 = 2.4 * 1e-10    # um^2, cross-section transition between Cr(4+)(d) levels in absorber
Na_total = 1*1e5  # um^-3, Cr(4+) ion concentration in tetraidal state

# Parameters for the third equation
lg = 1000   # um, thickness of active medium
la = 50  # um, absorber thickness (10-250 um)
sigma_a2 = 2.8 * 1e-11    # um^2, cross-section transition between Cr(4+)(d) levels in absorber
#gamma_i = alpha_p * lg  # absorption losses in active media
R1 = 1  # reflectivity of input mirror
#gamma_1 = -math.log(R1)  # losses on input mirror
R2 = 0.99   # reflectivity of output mirror
#gamma_2 = -math.log(R2)  # losses on output mirror
#gamma = gamma_i + 0.5 * (gamma_1 + gamma_2)  # losses
eps = 10**-13   # describes relative power of spontaneous emission, dimensionless
Ng_total = 1.66 * 1e5   # um^-3, Activator concentration (1.2 at %)
#lr = n * (lg + la)  # um, resonator optical length, lr = n(l_g + l_a)
#tau_r = 2 * lr / c0  # us, time of photon full pass in resonator (2l'/c0), t_r = 2*l_opt/c

# sigma_p = 0.495 * 10**-11		# um^2, effective cross-section absorption on excitation wavelength

#V = math.pi * rp**2 * lg  	# um^3, pumping beam volume in generating media
#Vg = 0.5 * math.pi * rl**2 * lg
#Veff = lr / la * Vg         # um^3, effective mode volume

# R - excitation speed
# Ae = pi * rp^2 / 4; % um^2
# Pa = Pi * (1 - exp(-2 * alpha_p * lg));
# R1 = Pi / (Ae * lg * h * nu_p / n * 8.2 *1e6);
# R1 = Pi * (1 - exp(-2 * lg * sigma_p * Ng_total)) / (Ae * lg * h * nu_p * 1e6);  % um^-3 * us^-1
# R1 = Pi * (Ng_total - exp(-2 * lg * sigma_p * Ng_total)) / (Ae * lg * h * nu_p * 1e6 * 0.75e6);  % um^-3 * us^-1
#eta = rl/rp
#R = eta * Pi / (V * h * nu_p * 1e6)

# Coefs
#Cg = sigma_g * c0 / Veff   # 1/us
#Ca = -1 * sigma_a1 * c0 / Veff   # 1/us
#g = 2 * sigma_g * lg    # um^3
#a1 = 2 * sigma_a1 * la  # um^3
#a2 = 2 * sigma_a2 * la  # um^3
#Cg_eps = eps * c0 * sigma_g * lg / lr  # um^3/us.

# Initial values
#t0 = 0
#tau_m = 50
#dt = 0.01
#x0 = numpy.array([Ng_total, Na_total, 0])
#x_out_initial = x0
#tau_initial = t0