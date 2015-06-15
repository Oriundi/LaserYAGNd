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

#from scipy import seterr
from scipy import integrate, optimize
#from scipy.signal import argrelmax
from laserYAGNd_ui_800x600 import *
import numpy
from PyQt4 import QtCore, QtGui
import pyqtgraph as pqg
#import sympy
#from multiprocessing import Pool
# import sqlalchemy


try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8

    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

version_file = open(os.path.join(os.path.curdir, _fromUtf8('VERSION')))
version = version_file.read().strip()


class MainApp(QtGui.QMainWindow):
    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.input_parameter.setCurrentIndex(0)
        self.open_default_config()

        self.connect(self.ui.button_calculate, QtCore.SIGNAL(_fromUtf8("clicked()")), self.app_run)
        self.connect(self.ui.actionOpen_config, QtCore.SIGNAL(_fromUtf8("activated()")), self.open_config)
        self.connect(self.ui.actionSave_config, QtCore.SIGNAL(_fromUtf8("activated()")), self.save_config)
        self.connect(self.ui.actionShow_Energy, QtCore.SIGNAL(_fromUtf8("activated()")), self.show_energy)

        font = QtGui.QFont()
        font.setFamily("Liberation Sans")
        font.setPointSize(11)
        font.setBold(False)
        font.setWeight(50)
        self.ui.plot_NdNa.setBackground('w')
        self.ui.plot_P.setBackground('w')
        self.ui.plot_NdNa.setTitle(_fromUtf8('Crystal properties'), font=font)
        self.ui.plot_P.setTitle(_fromUtf8('Power'), font='k')
        self.ui.plot_NdNa.setLabel('left', text=_fromUtf8('Concentration, cm<sup>-3</sup>,<br>'
                                                          'red - N<sub>g</sub>, green - N<sub>a</sub>'), font='k')
        self.ui.plot_NdNa.setLabel('bottom',  text=_fromUtf8('Time, μs'), font='k')
        self.ui.plot_P.setLabel('left', text=_fromUtf8('Power, W'), font='k')
        self.ui.plot_P.setLabel('bottom', text=_fromUtf8('Time, μs'), font='k')

    def open_default_config(self):
        cfgfile = configparser.ConfigParser()
        cfgfile.read(_fromUtf8('.config_default.cfg'))

        self.ui.value_tau_m.setProperty("value", cfgfile.get('Times', 'tau_m'))
        self.ui.value_dt.setProperty("value", cfgfile.get('Times', 'dt'))
        self.ui.value_tau_g.setProperty("value", cfgfile.get('Times', 'tau_g'))
        self.ui.value_tau_a.setProperty("value", cfgfile.get('Times', 'tau_a'))

        self.ui.value_lg.setProperty("value", cfgfile.get('Resonator', 'lg'))
        self.ui.value_la.setProperty("value", cfgfile.get('Resonator', 'la'))
        self.ui.value_rl.setProperty("value", cfgfile.get('Resonator', 'rl'))
        self.ui.value_n.setProperty("value", cfgfile.get('Resonator', 'n'))

        self.ui.value_sigma_g.setProperty("value", float(cfgfile.get('Cross-sections', 'sigma_g')) / 1e-19)
        self.ui.value_sigma_a1.setProperty("value", float(cfgfile.get('Cross-sections', 'sigma_a1')) / 1e-18)
        self.ui.value_sigma_a2.setProperty("value", float(cfgfile.get('Cross-sections', 'sigma_a2')) / 1e-19)

        #self.ui.value_Ng_total.setProperty("value", float(cfgfile.get('Crystal properties', 'Ng_total')) / 1e17)
        self.ui.value_Ng_total_perc.setProperty("value", cfgfile.get('Crystal properties', 'Ng_total_perc'))
        self.ui.value_T0.setProperty("text", float(cfgfile.get('Crystal properties', 'T0')))

        self.ui.value_alpha.setProperty("value", float(cfgfile.get('Losses', 'alpha')))
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

        self.ui.value_sigma_g.setProperty("value", float(cfgfile.get('Cross-sections', 'sigma_g')) / 1e-19)
        self.ui.value_sigma_a1.setProperty("value", float(cfgfile.get('Cross-sections', 'sigma_a1')) / 1e-18)
        self.ui.value_sigma_a2.setProperty("value", float(cfgfile.get('Cross-sections', 'sigma_a2')) / 1e-19)

        self.ui.value_Ng_total_perc.setProperty("value", cfgfile.get('Crystal properties', 'Ng_total_perc'))
        self.ui.value_T0.setProperty("text", float(cfgfile.get('Crystal properties', 'Na_total')))

        self.ui.value_alpha.setProperty("value", float(cfgfile.get('Losses', 'alpha')))
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
        cfgfile.set('Cross-sections', 'sigma_g', str(float(self.ui.value_sigma_g.text()) * 1e-19))
        cfgfile.set('Cross-sections', 'sigma_a1', str(float(self.ui.value_sigma_a1.text()) * 1e-18))
        cfgfile.set('Cross-sections', 'sigma_a2', str(float(self.ui.value_sigma_a2.text()) * 1e-19))

        cfgfile.add_section('Crystal properties')
        cfgfile.set('Crystal properties', 'Ng_total_perc', self.ui.value_Ng_total_perc.text())
        cfgfile.set('Crystal properties', 'Na_total', str(float(self.ui.value_T0.text())))

        cfgfile.add_section('Losses')
        cfgfile.set('Losses', 'alpha', str(float(self.ui.value_alpha.text())))
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
        T0 = self.ui.value_T0.valueFromText(self.ui.value_T0.text())
        Ng_total_perc = self.ui.value_Ng_total_perc.valueFromText(self.ui.value_Ng_total_perc.text())
        gamma, tau_r, R, Ng_total, Na_total, Cg, Ca, g, a1, a2, Cg_eps, nu_p = self.update_input_data()
        eq = []

        #Ng_total = 1.66e7
        #Na_total = 1e6

        # Xiao-Bass speed equations
        # eq[0] = R - Cg * x[0] * x[2] - x[0]/tau_g
        # eq[1] = Ca * x[1] * x[2] + (Na_total - x[1])/tau_a
        # eq[2] = (x[0]*g - x[1]*a1 - (Na_total-x[1])*a2 - 2*gamma) * x[2]/tau_r + (x[0]+Ng_total) * Cg_eps
        eq = numpy.array([R - Cg * x[0] * x[2] - x[0]/tau_g,
                           Ca * x[1] * x[2] + (Na_total - x[1])/tau_a,
                           (x[0]*g - x[1]*a1 - (Na_total-x[1])*a2 - 2*gamma)
                          * x[2]/tau_r + (x[0]+Ng_total)*Cg_eps])
        return eq

    def solve_equation(self, x_out, tau):
        tau_m = self.ui.value_tau_m.valueFromText(self.ui.value_tau_m.text())
        dt = self.ui.value_dt.valueFromText(self.ui.value_dt.text())

        progress = QtGui.QProgressDialog("Calculating...", "Cancel", 0, tau_m, self)
        progress.setWindowModality(QtCore.Qt.WindowModal)
        progress.setMinimumDuration(0)

        # ts = numpy.linspace(tau, dt, tau_m)
        # r = integrate.odeint(self.equation, x_out, ts)

        r = integrate.ode(self.equation).set_integrator('dop853', atol=1e-12, rtol=1e-12)
        r.set_initial_value(x_out, tau)
        while r.successful() and r.t < tau_m:
            r.integrate(r.t+dt)
            x_out = numpy.vstack([x_out, numpy.abs(r.y)])
            tau = numpy.vstack([tau, r.t])
            progress.setValue(r.t)
            if progress.wasCanceled():
                break
        progress.setValue(tau_m)
        if not r.successful():
            error_msg = QtGui.QMessageBox()
            error_msg.setWindowModality(QtCore.Qt.WindowModal)
            font = QtGui.QFont()
            font.setPointSize(14)
            font.setBold(True)
            font.setFamily("Liberation Sans")
            error_msg.setText('Please make time step smaller')
            error_msg.setFont(font)
            error_msg.setFocus()
            error_msg.exec_()
        nd = x_out[:, 0]
        na = x_out[:, 1]
        q = x_out[:, 2]
        return nd, na, q, tau

    def update_input_data(self):
        c0 = 3e4   #c0 = 3 * 1e8  # ucm/us, speed of light
        h = 6.626069 * 1e-34 * 1e6  # J*us, Plank constant
        eps = 10**-13   # describes relative power of spontaneous emission, Resonatorless

        lg = self.ui.value_lg.valueFromText(self.ui.value_lg.text())
        la = self.ui.value_la.valueFromText(self.ui.value_la.text())
        rl = self.ui.value_rl.valueFromText(self.ui.value_rl.text())
        n = self.ui.value_n.valueFromText(self.ui.value_n.text())

        sigma_g = self.ui.value_sigma_g.valueFromText(self.ui.value_sigma_g.text()) * 1e-19  # * 1e8
        sigma_a1 = self.ui.value_sigma_a1.valueFromText(self.ui.value_sigma_a1.text()) * 1e-18   # 1e8
        sigma_a2 = self.ui.value_sigma_a2.valueFromText(self.ui.value_sigma_a2.text()) * 1e-19  # 1e8

        Ng_total_perc = self.ui.value_Ng_total_perc.valueFromText(self.ui.value_Ng_total_perc.text())
        T0 = self.ui.value_T0.valueFromText(self.ui.value_T0.text())

        alpha = self.ui.value_alpha.valueFromText(self.ui.value_alpha.text())   # 1e-4
        R1 = self.ui.value_R1.valueFromText(self.ui.value_R1.text())
        R2 = self.ui.value_R2.valueFromText(self.ui.value_R2.text())

        Pi = self.ui.value_Pi.valueFromText(self.ui.value_Pi.text())
        wavelength = self.ui.value_wavelength.valueFromText(self.ui.value_wavelength.text())
        rp = self.ui.value_rp.valueFromText(self.ui.value_rp.text())

        gamma_i = alpha * lg  # absorption losses in active media
        gamma_1 = -math.log(R1)  # losses on input mirror
        gamma_2 = -math.log(R2)  # losses on output mirror
        gamma = gamma_i + 0.5 * (gamma_1 + gamma_2)  # losses
        lr = n * (lg + la)  # um, resonator optical length, lr = n(l_g + l_a)
        tau_r = 2 * lr / c0  # us, time of photon full pass in resonator (2l'/c0), t_r = 2*l_opt/c
        nu_p = c0 / (n * wavelength)  # 1/us, radiation excitation frequency

        V = math.pi * rp**2 * lg  	# um^3, pumping beam volume in generating media
        Vg = 1 / 4 * math.pi * rl**2 * lg
        Veff = lr / lg * Vg          # um^3, effective mode volume
        # print(lr, lg, Vg)

        eta = rl/rp
        R = eta * Pi / (V * h * nu_p) * 1e-6
        # print(R)

        Ng_total = 0
        Na_total = 0
        # print('T0 = %f, sigma_a1 = %f, la = %f' % (T0, sigma_a1 * 1e18, la))
        Ng_total = (3 * 4.55) / ((3 * 88.9 + 5 * 27.0 + 12 * 16.0) * 1.6726 * 1e-24) * Ng_total_perc / 100
        Na_total = numpy.log(T0) / (-1 * sigma_a1 * la)
        # print("Ng_total out = %f" % Ng_total)
        # print("Na_total out = %f" % Na_total)
        # Ng_total = 1.66e18
        # Na_total = 1e18

        # Coefs
        Cg = sigma_g * c0 / Veff   # 1/us
        Ca = -1 * sigma_a1 * c0 / Veff   # 1/us
        g = 2 * sigma_g * lg    # um^3
        a1 = 2 * sigma_a1 * la  # um^3
        a2 = 2 * sigma_a2 * la  # um^3
        Cg_eps = eps * c0 * sigma_g * lg / lr  # um^3/us.
        # print(Cg, Ca, g, a1, a2, Cg_eps)

        return gamma, tau_r, R, Ng_total, Na_total, Cg, Ca, g, a1, a2, Cg_eps, nu_p

    def show_energy(self):
        imp_energy = QtGui.QMessageBox()
        imp_energy.setWindowModality(QtCore.Qt.WindowModal)
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setBold(False)
        font.setFamily("Liberation Sans")
        text = ''
        try:
            for ii in xrange(0, imp):
                text += 'Energy of {} impulse is {:.2f} μJ \n'.format(ii+1, energy[ii]*1e6)

            imp_energy.setText(text)
            imp_energy.setFont(font)
            imp_energy.setFocus()
            imp_energy.exec_()
        except ValueError:
            pass

    def app_run(self):
        self.ui.plot_NdNa.clear()
        self.ui.plot_P.clear()
        error_msg = None
        global imp, energy

        font_error = QtGui.QFont()
        font_error.setPointSize(14)
        font_error.setBold(True)
        font_error.setFamily("Liberation Sans")

        h = 6.626069 * 1e-34 * 1e6  # J*us, Plank constant
        gamma, tau_r, R, Ng_total, Na_total, Cg, Ca, g, a1, a2, Cg_eps, nu_p = self.update_input_data()
        #print(Ng_total, Na_total)

        tau_g = self.ui.value_tau_g.valueFromText(self.ui.value_tau_g.text())
        dt = self.ui.value_dt.valueFromText(self.ui.value_dt.text())
        sigma_g = self.ui.value_sigma_g.valueFromText(self.ui.value_sigma_g.text()) * 1e-19   # * 1e8
        sigma_a1 = self.ui.value_sigma_a1.valueFromText(self.ui.value_sigma_a1.text()) * 1e-18   # * 1e8
        sigma_a2 = self.ui.value_sigma_a2.valueFromText(self.ui.value_sigma_a2.text()) * 1e-19   # * 1e8
        lg = self.ui.value_lg.valueFromText(self.ui.value_lg.text())
        la = self.ui.value_la.valueFromText(self.ui.value_la.text())
        alpha = self.ui.value_alpha.valueFromText(self.ui.value_alpha.text()) * 1e-3   # * 1e-4
        R1 = self.ui.value_R1.valueFromText(self.ui.value_R1.text())
        R2 = self.ui.value_R2.valueFromText(self.ui.value_R2.text())
        n = self.ui.value_n.valueFromText(self.ui.value_n.text())
        Ng_total_perc = self.ui.value_Ng_total_perc.valueFromText(self.ui.value_Ng_total_perc.text())
        T0 = self.ui.value_T0.valueFromText(self.ui.value_T0.text())
        lr = n * (lg + la)

        #Ngi = 1 / (sigma_g * lg) * (sigma_a1 * la * Na_total + alpha * lr - 1 / 2 * math.log(R1 * R2))
        #Ngi = 1.66e7
        #Ngmax = R * tau_g
        #Na_total = 1e6
        
        # print("tau_m = %f" % (tau_m))
        # print("tau_a = %f" % (tau_a))
        # print("tau_g = %f" % (tau_g))
        # print("tau_a = %f" % (tau_a))

        # print("Ng_total = %f" % (Ng_total * 1e-17))
        # print("Na_total = %f" % (Na_total * 1e-17))
        # print("dt = %f" % (dt))

        # print("x_out_initial = %f" % (x_out_initial))
        # print("tau_initial = %f" % (tau_initial))

        # print("sigma_a2 = %f" % (sigma_a2))
        # print("lg = %f" % (lg))
        # print("la = %f" % (la))
        # print("rl = %f" % (rl))
        # print("alpha = %f" % (alpha))
        # print("R2 = %f" % (R2))
        # print("n = %f" % (n))
        # print("Ng_total_perc = %f" % (Ng_total_perc))
        # print("T0 = %f" % (T0))
        # print("R1 = %f" % (R1))
        
        
        x_out_initial = numpy.array([Ng_total, Na_total, 0])
        tau_initial = 0
        Ng, Na, q, tau = self.solve_equation(x_out_initial, tau_initial)
        P = numpy.array(h * (3e4 / 0.0001064) / tau_r * math.log(1/R2) * q * 1e6)

        penNd = pqg.mkPen(width=1, color='r')
        penNa = pqg.mkPen(width=1, color='g')
        penP = pqg.mkPen(width=1, color='b')
        self.ui.plot_P.plot(tau[:, 0],  P, pen=penP)
        #self.ui.plot_NdNa.plot(tau[:, 0], Ng * 1e12, pen=penNd)
        #self.ui.plot_NdNa.plot(tau[:, 0], Na * 1e12, pen=penNa)
        self.ui.plot_NdNa.plot(tau[:, 0], Ng, pen=penNd)
        self.ui.plot_NdNa.plot(tau[:, 0], Na, pen=penNa)
'''
        # Analisys
        #T0 = math.exp(-sigma_a1 * la * Na_total)
        #T = numpy.exp(-sigma_a1 * la * Na)
        #self.ui.plot_P.plot(tau[:, 0], T * 100, pen=penNd)

        #T00 = 0.9
        #Na0 = - math.log(T00) / (sigma_a1 * la)
        #Na0 = 1e8
        #T000 = numpy.exp(-sigma_a1 * la * Na0)
        #print(Na0)
        #print(T000)

        local_max = numpy.array(argrelmax(P))
        for ii in xrange(numpy.size(local_max)):
            index = local_max[0, ii]
            if P[index] < 1e-3:
                local_max[0, ii] = 0
        mask = numpy.where(local_max == 0)
        local_max = numpy.delete(local_max, mask[1])
        imp = numpy.size(local_max)
        imp_length = 0.1
        cut = int(imp_length/dt)
        max_P = []
        energy_cut = math.exp(1) ** 2
        for ii in xrange(imp):
            value = local_max[ii]
            index_max = numpy.zeros([2*cut])
            P_local_max = P[value]
            if value-cut < 0:
                jj0 = 0
            else:
                jj0 = local_max[ii] - cut
            if value+cut > numpy.size(P)-1:
                jjn = numpy.size(P) - 1
            else:
                jjn = local_max[ii] + cut
            for jj in xrange(jj0, jjn):
                index_max[jj-jj0] = jj
            index_max = numpy.delete(index_max, numpy.where(P[jj0:jjn] < P_local_max/energy_cut))
            if numpy.size(index_max) != 0:
                max_P.append(index_max)
        for ii in xrange(imp):
            index_max_zeros = numpy.array(numpy.where(max_P[ii] == 0))
            max_P[ii] = numpy.delete(max_P[ii], index_max_zeros)
        #print(max_P)

        dt = self.ui.value_dt.valueFromText(self.ui.value_dt.text())
        energy = numpy.zeros(imp)
        try:
            for ii in xrange(imp):
                x = P[max_P[ii][0]:max_P[ii][numpy.size(max_P[ii]) - 1] + 1]
                if numpy.size(x) < 3:
                    if not error_msg:
                        error_msg = QtGui.QMessageBox()
                        error_msg.setText('Not enough point to draw impulse properly. Please make time step smaller.')
                        error_msg.setFont(font_error)
                        error_msg.exec_()
                        break
                energy[ii] = integrate.simps(x, dx=dt) * 1e-6
            #print(energy)
        except IndexError:
            if max_P[:][0] != max_P[:][1]-1 or max_P[:][numpy.size(max_P)-2] != max_P[:][numpy.size(max_P)-1]-1:
                if not error_msg:
                    error_msg = QtGui.QMessageBox()
                    error_msg.setText('Not enough point to draw impulse properly. Please make time step smaller.')
                    error_msg.setFont(font_error)
                    error_msg.exec_()
            for ii in xrange(1, numpy.size(max_P)-2):
                if not error_msg:
                    if max_P[:][ii] != max_P[:][ii-1] + 1 and max_P[:][ii] != max_P[:][ii+1] - 1:
                        error_msg = QtGui.QMessageBox()
                        error_msg.setText('Not enough point to draw impulse properly. Please make time step smaller')
                        error_msg.setFont(font_error)
                        error_msg.exec_()
            if not error_msg:
                error_msg = QtGui.QMessageBox()
                error_msg.setText('No impulses')
                error_msg.setFont(font_error)
                error_msg.exec_()
'''


def main():
    app = QtGui.QApplication(sys.argv)
    #app.setStyle('Plastique')
    window = MainApp()
    window.setWindowTitle("Laser YAG:Nd")
    window.show()
    sys.exit(app.exec_())
if __name__ == '__main__':
    main()
