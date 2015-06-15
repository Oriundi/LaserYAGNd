import math
import numpy

# global R

tau_g = 230  # um, lifetime of high laser level (0.23ms)
tau_a = 3.5  # um, lifetime of excited level 3T2 of ions Cr(4+)(d) (3.5us)

# Parameters for the first equation
sigma_g = 3.5 * 1e-11  # um^2 efficient cross-section of laser transition
c0 = 3 * 1e8  # um^2/us, speed of light
rp = 300    # um, mode radius
rl = 300    # um, laser beam radius
Pi = 0.2  # W, pumping power
alpha_p = 5 * 1e-7  # 1/um, absorption coefficient on the wavelength
h = 6.626069 * 1e-34 * 1e6  # J*us, Plank constant
h_ = 1.054571 * 1e-34 * 1e6
wavelength_p = 0.808  # um, radiation excitation wavelength
n = 1.8169  # YAG:Nd refractive index at 1.064 um
nu_p = c0 / (n * wavelength_p)  # 1/us, radiation excitation frequency
wavelength_l = 1.064
nu_l = c0/wavelength_l

# Parameters for the second equation
sigma_a1 = 2.4 * 1e-10    # um^2, cross-section transition between Cr(4+)(d) levels in absorber
# Na_total = 1*1e5  # um^-3, Cr(4+) ion concentration in tetraidal state

T0 = 0.9988   # Initial absorber transmission ( 0 ... 1)
Ng_total_perc = 0.001208    # 0.001208   # Initial Activator concentration

# Parameters for the third equation
lg = 1000  # um, thickness of active medium
la = 50  # um, absorber thickness (10-250 um)
sigma_a2 = 2.8 * 1e-11    # um^2, cross-section transition between Cr(4+)(d) levels in absorber
gamma_i = alpha_p * lg  # absorption losses in active media
R1 = 1  # reflectivity of input mirror
gamma_1 = -math.log(R1)  # losses on input mirror
R2 = 0.99   # reflectivity of output mirror

gamma_2 = -math.log(R2)  # losses on output mirror
gamma = gamma_i + 0.5 * (gamma_1 + gamma_2)  # losses
eps = 10**-13   # describes relative power of spontaneous emission, dimensionless
# Ng_total = 1.66*1e5   # um^-3, Activator concentration (1.2 at %)
lr = n * (lg + la)  # um, resonator optical length, lr = n(l_g + l_a)
tau_r = 2 * lr / c0  # us, time of photon full pass in resonator (2l'/c0), t_r = 2*l_opt/c

sigma_p = 0.495 * 10**-11		# um^2, effective cross-section absorption on excitation wavelength

V = math.pi * rp**2 * lg  	# um^3, pumping beam volume in generating media
Vg = 1/4 * math.pi * rl**2 * lg   # mode volume in generating medium
Veff = lr / lg * Vg         # um^3, effective mode volume
# print(lr, lg, Vg)

Ng_total = (3 * 4.55) / ((3 * 88.9 + 5 * 27.0 + 12 * 16.0) * 1.6726 * 1e-24) * Ng_total_perc / 100 * 1e-12
# print('T0 = %f, sigma_a1 = %f, la = %f' % (T0, sigma_a1 * 1e10, la))
Na_total = math.log(T0) / (-1 * sigma_a1 * la)

# R - excitation speed
A = math.pi * rp**2   # um^2
# Pa = Pi * (1 - exp(-2 * alpha_p * lg));
# R1 = Pi / (Ae * lg * h * nu_p / n * 8.2 *1e6);
Pi_losses = 1 - math.exp(-2 * sigma_p * Ng_total * lg)
# print("Pi losses %f" % Pi_losses)
R11 = Pi * Pi_losses  # W
R22 = A * lg * h * nu_p
R_ = R11 / R22 * 1e-6
# R1 = Pi * (Ng_total - exp(-2 * lg * sigma_p * Ng_total)) / (Ae * lg * h * nu_p * 1e6 * 0.75e6);  # um^-3 * us^-1
eta = 1
R = eta * Pi / (V * h * nu_p) * 1e-6
# print(R)

# Coefs
Cg = sigma_g * c0 / Veff  # 1/us
Ca = -1 * sigma_a1 * c0 / Veff   # 1/us
g = 2 * sigma_g * lg  # um^3
a1 = 2 * sigma_a1 * la    # um^3
a2 = 2 * sigma_a2 * la    # um^3
Cg_eps = eps * c0 * sigma_g * lg / lr  # um^3/us.

# print(Cg, Ca, g, a1, a2, Cg_eps)

# Initial values
time_start = 0
time_end = 50
dt = 0.01

x0 = numpy.array([Ng_total, Na_total, 0])
x_out_initial = x0
tau_initial = time_end

