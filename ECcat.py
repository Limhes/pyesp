#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 11:28:36 2019

@author: limhes
"""

import PyESP
import math
import matplotlib.pyplot as plt
plt.style.use('ggplot')

max_i = []
max_v = []
concs_sub = [0, 1, 3, 10, 30]

for conc_sub in concs_sub:
    PyESP.setup()
    
    step_potential = 0.002 # V
    scan_rate = 0.1 # V/s
    
    PyESP.set_params(SI=int(step_potential*1000.))
    PyESP.set_params(ST=step_potential/scan_rate)
    
    A = PyESP.addSpecies(1.0, 1.0e-5)
    B = PyESP.addSpecies(0.0, 1.0e-5)
    substrate = PyESP.addSpecies(conc_sub, 1.0e-5)
    product = PyESP.addSpecies(0.0, 1.0e-5)
    redox1 = PyESP.addRedox(A, B, 1, 0, 1e3, 0.5)
    chem1 = PyESP.addChemical(B, substrate, A, product, 1.0e2, 0.)
    
    [potential, current] = PyESP.simulate()
    print(conc_sub)
    if conc_sub == 0: conc_div = 1.0
    else: conc_div = conc_sub
    plt.plot(potential, [i*1.0e6/math.sqrt(scan_rate)/conc_div for i in current], label='[substrate] = {:.1f} mM'.format(conc_sub))
    
    max_i.append(max(current))
    max_v.append(potential[current.index(max(current))])
    
    PyESP.destroy()
    
    
plt.xlabel('Potential [V]')
plt.ylabel('Current [uA]/sqrt(scan rate [V/s])/([substrate] [mM])')
#plt.xlim(-0.55, 0.55)
plt.legend()
plt.show()

plt.figure()
plt.plot(concs_sub[1:], max_i[1:])
plt.xlabel('[substrate] [mM]')
plt.ylabel('i_p [uA]')
plt.show()




