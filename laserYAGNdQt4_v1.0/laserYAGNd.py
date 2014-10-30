#!/usr/bin/env python3

'''
Documentation, License etc.
Xiao-Bass speed equations for
YAG:Nd laser
@package laserYAGNd
'''

import sys
import os
import math
import configparser
from scipy.integrate import ode
from laserYAGNd_ui import *
import numpy
from PyQt4 import QtCore, QtGui
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar

version_file = open(os.path.join(os.path.curdir, 'VERSION'))
version = version_file.read().strip()

class fig_etmpty(FigureCanvas):
	def __init__(self):
		# Standard Matplotlib code to generate the plot
		self.fig0 = Figure()
		self.axes0 = self.fig0.add_subplot(111)
		FigureCanvas.__init__(self, self.fig0)


class fig_NdNa(FigureCanvas):
	def __init__(self, parent, tau, Nd, Na):
		# Standard Matplotlib code to generate the plot
		self.fig1 = Figure()
		self.axes1 = self.fig1.add_subplot(111)
		self.axes1.hold(True)

		self.axes1.plot(tau, Nd)
		self.axes1.plot(tau, Na)
		# initialize the canvas where the Figure renders into
		FigureCanvas.__init__(self, self.fig1)
		self.setParent(parent)
		FigureCanvas.setSizePolicy(self, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
		FigureCanvas.updateGeometry(self)


class fig_P(FigureCanvas):
	def __init__(self, parent, tau, P):
		# Standard Matplotlib code to generate the plot
		self.fig2 = Figure()
		self.axes2 = self.fig2.add_subplot(111)
		self.axes2.hold(False)

		# Nd, Na, q, P, tau = laserYAGNd(x_out_initial, tau_initial)

		self.axes2.plot(tau, P)
		# initialize the canvas where the Figure renders into
		FigureCanvas.__init__(self, self.fig2)
		self.setParent(parent)
		FigureCanvas.setSizePolicy(self, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
		FigureCanvas.updateGeometry(self)


class NavToolbar(NavigationToolbar):
	# only display the buttons we need
	toolitems = [t for t in NavigationToolbar.toolitems if t[0] in ('Home', 'Pan', 'Zoom', 'Save')]


class MainApp(QtGui.QMainWindow):
	def __init__(self, parent=None):
		QtGui.QMainWindow.__init__(self, parent)
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self)
		self.ui.input_parameter.setCurrentIndex(0)

		self.initialize_figures()
		self.open_default_config()

		self.connect(self.ui.button_calculate, QtCore.SIGNAL("clicked()"), self.app_run)
		self.connect(self.ui.actionOpen_config, QtCore.SIGNAL("activated()"), self.open_config)
		self.connect(self.ui.actionSave_config, QtCore.SIGNAL("activated()"), self.save_config)

	def open_default_config(self):
		cfgfile = configparser.ConfigParser()
		cfgfile.read('.config_default.cfg')

		self.ui.value_tau_m.setProperty("value", cfgfile.get('Times', 'tau_m'))
		self.ui.value_dt.setProperty("value", cfgfile.get('Times', 'dt'))
		self.ui.value_tau_g.setProperty("value", cfgfile.get('Times', 'tau_g'))
		self.ui.value_tau_a.setProperty("value", cfgfile.get('Times', 'tau_a'))

		self.ui.value_lg.setProperty("value", cfgfile.get('Resonator', 'lg'))
		self.ui.value_la.setProperty("value", cfgfile.get('Resonator', 'la'))
		self.ui.value_rl.setProperty("value", cfgfile.get('Resonator', 'rl'))
		self.ui.value_n.setProperty("value", cfgfile.get('Resonator', 'n'))

		self.ui.value_sigma_g.setProperty("value", cfgfile.get('Cross-sections', 'sigma_g'))
		self.ui.value_sigma_a1.setProperty("value", cfgfile.get('Cross-sections', 'sigma_a1'))
		self.ui.value_sigma_a2.setProperty("value", cfgfile.get('Cross-sections', 'sigma_a2'))

		self.ui.value_Ng_total.setProperty("value", cfgfile.get('Concentrations', 'Ng_total'))
		self.ui.value_Ng_total_perc.setProperty("value", cfgfile.get('Concentrations', 'Ng_total_perc'))
		self.ui.value_Na_total.setProperty("value", cfgfile.get('Concentrations', 'Na_total'))

		self.ui.value_alpha_p.setProperty("value", cfgfile.get('Losses', 'alpha_p'))
		self.ui.value_R1.setProperty("value", cfgfile.get('Losses', 'R1'))
		self.ui.value_R2.setProperty("value", cfgfile.get('Losses', 'R2'))

		self.ui.value_Pi.setProperty("value", cfgfile.get('Pumping', 'Pi'))
		self.ui.value_wavelength.setProperty("value", cfgfile.get('Pumping', 'wavelength'))
		self.ui.value_rp.setProperty("value", cfgfile.get('Pumping', 'rp'))

	def open_config(self):
		FileDialog = QtGui.QFileDialog()
		file_filter = 'Config Files (*.cfg)'
		default_dir = os.path.expanduser("~")
		filename = QtGui.QFileDialog.getOpenFileName(FileDialog, "Open Configuration File", default_dir, file_filter)

		cfgfile = configparser.ConfigParser()
		cfgfile.read(filename)

		self.ui.value_tau_m.setProperty("value", cfgfile.get('Times', 'tau_m'))
		self.ui.value_dt.setProperty("value", cfgfile.get('Times', 'dt'))
		self.ui.value_tau_g.setProperty("value", cfgfile.get('Times', 'tau_g'))
		self.ui.value_tau_a.setProperty("value", cfgfile.get('Times', 'tau_a'))

		self.ui.value_lg.setProperty("value", cfgfile.get('Resonator', 'lg'))
		self.ui.value_la.setProperty("value", cfgfile.get('Resonator', 'la'))
		self.ui.value_rl.setProperty("value", cfgfile.get('Resonator', 'rl'))
		self.ui.value_n.setProperty("value", cfgfile.get('Resonator', 'n'))

		self.ui.value_sigma_g.setProperty("value", cfgfile.get('Cross-sections', 'sigma_g'))
		self.ui.value_sigma_a1.setProperty("value", cfgfile.get('Cross-sections', 'sigma_a1'))
		self.ui.value_sigma_a2.setProperty("value", cfgfile.get('Cross-sections', 'sigma_a2'))

		self.ui.value_Ng_total.setProperty("value", cfgfile.get('Concentrations', 'Ng_total'))
		self.ui.value_Ng_total_perc.setProperty("value", cfgfile.get('Concentrations', 'Ng_total_perc'))
		self.ui.value_Na_total.setProperty("value", cfgfile.get('Concentrations', 'Na_total'))

		self.ui.value_alpha_p.setProperty("value", cfgfile.get('Losses', 'alpha_p'))
		self.ui.value_R1.setProperty("value", cfgfile.get('Losses', 'R1'))
		self.ui.value_R2.setProperty("value", cfgfile.get('Losses', 'R2'))

		self.ui.value_Pi.setProperty("value", cfgfile.get('Pumping', 'Pi'))
		self.ui.value_wavelength.setProperty("value", cfgfile.get('Pumping', 'wavelength'))
		self.ui.value_rp.setProperty("value", cfgfile.get('Pumping', 'rp'))

	def save_config(self):
		FileDialog = QtGui.QFileDialog()
		file_filter = 'Config Files (*.cfg)'
		default_dir = os.path.expanduser("~")
		filename = QtGui.QFileDialog.getSaveFileName(FileDialog, "Save Configuration File", default_dir, file_filter)

		cfgfile = configparser.ConfigParser()
		cfgfile1 = open(filename, 'w')

		cfgfile.add_section('Times')
		cfgfile.set('Times', 'tau_m', self.ui.value_tau_m.text())
		cfgfile.set('Times', 'dt', self.ui.value_dt.text())
		cfgfile.set('Times', 'tau_g', self.ui.value_tau_g.text())
		cfgfile.set('Times', 'tau_a', self.ui.value_tau_a.text())

		cfgfile.add_section('Resonator')
		cfgfile.set('Resonator', 'lg', self.ui.value_lg.text())
		cfgfile.set('Resonator', 'la', self.ui.value_la.text())
		cfgfile.set('Resonator', 'rl', self.ui.value_rl.text())
		cfgfile.set('Resonator', 'n', self.ui.value_n.text())

		cfgfile.add_section('Cross-sections')
		cfgfile.set('Cross-sections', 'sigma_g', self.ui.value_sigma_g.text())
		cfgfile.set('Cross-sections', 'sigma_a1', self.ui.value_sigma_a1.text())
		cfgfile.set('Cross-sections', 'sigma_a2', self.ui.value_sigma_a2.text())

		cfgfile.add_section('Concentrations')
		cfgfile.set('Concentrations', 'Ng_total', self.ui.value_Ng_total.text())
		cfgfile.set('Concentrations', 'Ng_total_perc', self.ui.value_Ng_total_perc.text())
		cfgfile.set('Concentrations', 'Na_total', self.ui.value_Na_total.text())

		cfgfile.add_section('Losses')
		cfgfile.set('Losses', 'alpha_p', self.ui.value_alpha_p.text())
		cfgfile.set('Losses', 'R1', self.ui.value_R1.text())
		cfgfile.set('Losses', 'R2', self.ui.value_R2.text())

		cfgfile.add_section('Pumping')
		cfgfile.set('Pumping', 'pi', self.ui.value_Pi.text())
		cfgfile.set('Pumping', 'wavelength', self.ui.value_wavelength.text())
		cfgfile.set('Pumping', 'rp', self.ui.value_rp.text())

		cfgfile.write(cfgfile1)

	def equation(self, t, x):
		tau_g = self.ui.value_tau_g.valueFromText(self.ui.value_tau_g.text())
		tau_a = self.ui.value_tau_a.valueFromText(self.ui.value_tau_a.text())
		Ng_total = self.ui.value_Ng_total.valueFromText(self.ui.value_Ng_total.text()) * 1e-12
		Na_total = self.ui.value_Na_total.valueFromText(self.ui.value_Na_total.text()) * 1e-12
		gamma, tau_r, R, Cg, Ca, g, a1, a2, Cg_eps, nu_p = self.update_input_data()

		# Xiao-Bass speed equations
		eq = numpy.array([R - Cg * x[0] * x[2] - x[0]/tau_g,
					Ca * x[1] * x[2] + (Na_total - x[1])/tau_a,
					(x[0]*g - x[1]*a1 - (Na_total-x[1])*a2 - 2*gamma)
					* x[2]/tau_r + (x[0]+Ng_total)*Cg_eps])
		return eq

	def laserYAGNd(self, x_out, tau):
		tau_m = self.ui.value_tau_m.valueFromText(self.ui.value_tau_m.text())
		dt = self.ui.value_dt.valueFromText(self.ui.value_dt.text())

		progress = QtGui.QProgressDialog("Calculating...", "Cancel", 0, tau_m, self)
		progress.setWindowModality(QtCore.Qt.WindowModal)
		progress.setMinimumDuration(0)

		r = ode(self.equation).set_integrator('dopri5', atol=1e-12, rtol=1e-12)
		r.set_initial_value(x_out, tau)
		while r.successful() and r.t < tau_m:
			r.integrate(r.t+dt)
			x_out = numpy.vstack([x_out, numpy.abs(r.y)])
			tau = numpy.vstack([tau, r.t])
			progress.setValue(r.t)
			if progress.wasCanceled():
				break
		progress.setValue(tau_m)
		nd = x_out[:, 0]
		na = x_out[:, 1]
		q = x_out[:, 2]
		return nd, na, q, tau

	def update_input_data(self):
		c0 = 3 * 1e8  # um^2/us, speed of light
		h = 6.626069 * 1e-34 * 1e6  # J*us, Plank constant
		eps = 10**-13   # describes relative power of spontaneous emission, Resonatorless

		lg = self.ui.value_lg.valueFromText(self.ui.value_lg.text())
		la = self.ui.value_la.valueFromText(self.ui.value_la.text())
		rl = self.ui.value_rl.valueFromText(self.ui.value_rl.text())
		n = self.ui.value_n.valueFromText(self.ui.value_n.text())

		sigma_g = self.ui.value_sigma_g.valueFromText(self.ui.value_sigma_g.text()) * 1e8
		sigma_a1 = self.ui.value_sigma_a1.valueFromText(self.ui.value_sigma_a1.text()) * 1e8
		sigma_a2 = self.ui.value_sigma_a2.valueFromText(self.ui.value_sigma_a2.text()) * 1e8

		# Ng_total_perc = self.ui.value_Ng_total_perc.valueFromText(self.ui.value_Ng_total_perc.text())

		alpha_p = self.ui.value_alpha_p.valueFromText(self.ui.value_alpha_p.text()) * 1e-4
		R1 = self.ui.value_R1.valueFromText(self.ui.value_R1.text())
		R2 = self.ui.value_R2.valueFromText(self.ui.value_R2.text())

		Pi = self.ui.value_Pi.valueFromText(self.ui.value_Pi.text())
		wavelength = self.ui.value_wavelength.valueFromText(self.ui.value_wavelength.text())
		rp = self.ui.value_rp.valueFromText(self.ui.value_rp.text())

		gamma_i = alpha_p * lg  # absorption losses in active media
		gamma_1 = -math.log(R1)  # losses on input mirror
		gamma_2 = -math.log(R2)  # losses on output mirror
		gamma = gamma_i + 0.5 * (gamma_1 + gamma_2)  # losses
		lr = n * (lg + la)  # um, resonator optical length, lr = n(l_g + l_a)
		tau_r = 2 * lr / c0  # us, time of photon full pass in resonator (2l'/c0), t_r = 2*l_opt/c
		nu_p = c0 / wavelength  # 1/us, radiation excitation frequency

		V = math.pi * rp**2 * lg  	# um^3, pumping beam volume in generating media
		Vg = 0.5 * math.pi * rl**2 * lg
		Veff = lr / la * Vg         # um^3, effective mode volume

		eta = rl/rp
		R = eta * Pi / (V * h * nu_p * 1e6)

		# Coefs
		Cg = sigma_g * c0 / Veff   # 1/us
		Ca = -1 * sigma_a1 * c0 / Veff   # 1/us
		g = 2 * sigma_g * lg    # um^3
		a1 = 2 * sigma_a1 * la  # um^3
		a2 = 2 * sigma_a2 * la  # um^3
		Cg_eps = eps * c0 * sigma_g * lg / lr  # um^3/us.

		return gamma, tau_r, R, Cg, Ca, g, a1, a2, Cg_eps, nu_p

	def initialize_figures(self):
		global vbl1, vbl2
		global qmc1, qmc2, ntb1, ntb2

		vbl1 = QtGui.QVBoxLayout(self.ui.plot_ndna)
		qmc1 = fig_etmpty()
		ntb1 = NavToolbar(qmc1, self.ui.plot_ndna)
		ntb1.setParent(None)
		vbl1.addWidget(qmc1)

		vbl2 = QtGui.QVBoxLayout(self.ui.plot_P)
		qmc2 = fig_etmpty()
		ntb2 = NavToolbar(qmc2, self.ui.plot_P)
		ntb2.setParent(None)
		vbl2.addWidget(qmc2)

	def app_run(self):
		global vbl1, vbl2
		global qmc1, qmc2, ntb1, ntb2

		h = 6.626069 * 1e-34 * 1e6  # J*us, Plank constant
		gamma, tau_r, R, Cg, Ca, g, a1, a2, Cg_eps, nu_p = self.update_input_data()
		R2 = self.ui.value_R2.valueFromText(self.ui.value_R2.text())

		Ng_total = self.ui.value_Ng_total.valueFromText(self.ui.value_Ng_total.text()) * 1e-12
		Na_total = self.ui.value_Na_total.valueFromText(self.ui.value_Na_total.text()) * 1e-12
		x_out_initial = numpy.array([Ng_total, Na_total, 0])
		tau_initial = 0

		Nd, Na, q, tau = self.laserYAGNd(x_out_initial, tau_initial)
		P = h * nu_p / tau_r * math.log(1/R2) * q * 1e6
		print(max(P))

		qmc1.setParent(None)
		ntb1.setParent(None)
		qmc1 = fig_NdNa(self.ui.plot_ndna, tau, Nd, Na)
		ntb1 = NavToolbar(qmc1, self.ui.plot_ndna)
		vbl1.addWidget(qmc1)
		vbl1.addWidget(ntb1)
		ntb1.setFocus()

		qmc2.setParent(None)
		ntb2.setParent(None)
		qmc2 = fig_P(self.ui.plot_P, tau, P)
		ntb2 = NavToolbar(qmc2, self.ui.plot_P)
		vbl2.addWidget(qmc2)
		vbl2.addWidget(ntb2)
		ntb2.setFocus()


def main():
	app = QtGui.QApplication(sys.argv)
	#app.setStyle('Oxigen')
	window = MainApp()
	window.setWindowTitle("Laser YAG:Nd")
	window.show()
	sys.exit(app.exec_())


if __name__ == '__main__':
	main()
