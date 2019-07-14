#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 12:43:43 2019

@author: limhes
"""

import PyESP
import math
import matplotlib.pyplot as plt
plt.style.use('ggplot')

PyESP.setup()
step_potential = 0.002 # V
PyESP.set_params(SI=int(step_potential*1000.),AM=1,AP=0)
A = PyESP.addSpecies(1.0, 1.0e-5)
B = PyESP.addSpecies(0.0, 1.0e-5)
C = PyESP.addSpecies(0.0, 1.0e-5)
redox1 = PyESP.addRedox(A, B, 1, -200, 1e3, 0.5)
chem1 = PyESP.addChemical(B, 0, C, 0, 10., 0.)

#for scan_rate in [0.1, 0.3, 1.0, 3.0, 10.]:
for scan_rate in [1.0]:
    PyESP.set_params(ST=step_potential/scan_rate)
    [potential, current] = PyESP.simulate()
    plt.plot(potential, [i*1.0e6/math.sqrt(scan_rate) for i in current], label='Scan rate {:.1f} V/s'.format(scan_rate))

plt.xlabel('Potential [V]')
plt.ylabel('Current [uA]/sqrt(scan rate [V/s])')
plt.xlim(-0.55, 0.2)
plt.legend()
plt.show()

PyESP.destroy()
