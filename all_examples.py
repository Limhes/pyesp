#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 11:05:22 2019

@author: limhes
"""

import PyESP
import enum, math
import matplotlib.pyplot as plt
plt.style.use('ggplot')

class Technique(enum.Enum):
    mode_CV = 0
    mode_SWV = 1
    mode_CA = 2
    mode_SDC = 3
    def __int__(self):
        return self.value

class IR_mode(enum.Enum):
    ir_none = 0
    ir_ru = 1
    ir_dl = 2
    def __int__(self):
        return self.value
    

# default experimental variables are loaded on setup (necessary to force reload of module --> needs some work):
PyESP.setup()
# experimental variables can be changed afterwards:
PyESP.set_params(IP=1000, FP=1000) # change initial and final potential
step_potential = 0.002 # V
PyESP.set_params(SI=int(step_potential*1000.))
#PyESP.set_params(CP=-200, CT=10) # conditioning out-of-equilibrium
#PyESP.set_params(V2=500, FP=0, ncycle=3) # cycling back over the re-oxidation wave
#PyESP.set_params(NC=5, SC=5) # take last cycle
#PyESP.set_params(AM=0) # sample at all four points of the step
#PyESP.set_params(IR=int(IR_mode.ir_dl), DL=3e-6, RU=10e5) # works
#PyESP.set_params(TE=150.) # very cold chemistry
#PyESP.set_params(Mode=int(Technique.mode_SWV), FR=50., PH=25) # square wave voltammetry --> works
#PyESP.set_params(Mode=int(Technique.mode_CA), V1=500, V2=-500, VD=10, T2=200) #  chronoamperometry --> works
#PyESP.set_params(Mode=int(Technique.mode_SDC), WE=1, AR=1.0) # polarography --> works

# addSpecies(conc [mM], diff_const [cm2/s])
A = PyESP.addSpecies(1.0, 1.0e-5)
B = PyESP.addSpecies(0.0, 1.0e-5)
C = PyESP.addSpecies(0.0, 1.0e-5)
D = PyESP.addSpecies(0.0, 1.0e-5)
CO2 = PyESP.addSpecies(10.0, 1.0e-5)

# addRedox(ox, red, n, E [mV], k_e [cm/s], alpha)
# addChemical(a, b, c, d, k_f, k_b)
# where reaction is: a + b <-> c + d and k_x is in 1/s (1st order) or 1/(mM*s) (2nd order)

# ErCi mechanism:
redox1 = PyESP.addRedox(A, B, 1, -200, 1e3, 0.5)
chem1 = PyESP.addChemical(B, 0, C, 0, 1.0e3, 0.)

# square scheme mechanism with homogeneous reactant:
#redox1 = PyESP.addRedox(A, B, 1, -200, 1e3, 0.5)
#redox2 = PyESP.addRedox(D, C, 1, 200, 1e3, 0.5)
#chem1 = PyESP.addChemical(B, CO2, C, 0, 1.0e2, 0.1)
#chem2 = PyESP.addChemical(A, CO2, D, 0, 1.0, 1.0e3)

# plot CVs at three different scan rates (for CA, this needs to be adapted):
for scan_rate in [0.1, 0.3, 1.0, 3.0, 10.]:
    # change scan rate
    PyESP.set_params(ST=step_potential/scan_rate)
    # run simulation
    [potential, current] = PyESP.simulate()
    # plot results
    plt.plot(potential, [i*1.0e6/math.sqrt(scan_rate) for i in current], label='Scan rate {:.1f} V/s'.format(scan_rate))

plt.xlabel('Potential [V]')
plt.ylabel('Current [uA]/sqrt(scan rate [V/s])')
plt.xlim(-0.55, 0.2)
plt.legend()
plt.show()

# free memory --> temporary approach (might not even be necessary...)
PyESP.destroy()
