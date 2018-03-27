#!/usr/bin/python

### note: run as  python optimize.py  -o Infections  100000000 0.2  40000000  0.8
### 100, 40 mil : # of vaccine doses
### 0.2 0.8: vaccine efficacies

import sys
sys.path.insert(0, r'./Influenza')
import Optimization
import csv


#import Simulation 
#s = Simulation.run_Simulation()
#vacsUsed = s.simulateWithVaccine([0],[0.7012, 0.5514, 0.3264, 0.4528, 0.6532], 0.1)
#print s.outputInfo()


o = Optimization.Optimization()

#if not o.options.quiet:
#    from miscellany import formatVacSchedule
#    print 'Vaccine schedule:\t %s' % formatVacSchedule(o.options.vacSchedule)
#    print 'Objective:\t\t %s' % o.options.objective

o.optimize()

#print o.outputInfo()
print o.short_output()


