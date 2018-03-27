#!/usr/bin/python
import sys
sys.path.insert(0, r'./Influenza')
import numpy as np
import Simulation 



s = Simulation.run_Simulation(paramValues = {"vacEfficacy":[0.2, 0.8]})

vacsUsed = s.simulateWithVaccine([0]*32, [0.2, 0.8])

#if not s.options.quiet:
#    print 'Doses:\t\t\t %s' % ', '.join('%g' % v for v in vacsUsed)
#    print

s.outputInfo()
