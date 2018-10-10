#!/usr/bin/python
import sys
sys.path.insert(0, r'../Influenza')
import numpy as np
np.warnings.filterwarnings('ignore')
import pandas as pd
import Simulation
import csv
import multiprocessing as mp
from multiprocessing import Pool
import itertools
import scipy.optimize
from scipy.optimize import fmin_cobyla
from scipy.optimize import basinhopping
import time
##########################################
result_list = []
def log_results(result):
    result_list.append(result)
################################################################
def universal_worker(input_pair):
    function, args0 = input_pair
    args = list(args0[0]) + [args0[1]]
    return function(*args)

def pool_args(function, *args):
    return zip(itertools.repeat(function), zip(itertools.repeat(args[:-1]) ,args[-1]))
   
###########################################################################

def run_efficacy_simulation_incidence(paramList, efficacy, doses, sub_iter):
    
    vacDoses = [doses,0]
    vacEfficacy = [efficacy,0]
    s = Simulation.run_Simulation(paramValues = {"vacEfficacy": vacEfficacy, "betaList": paramList[0:3], "susceptibility_H1": paramList[3:7], "susceptibility_H3": paramList[7:11], "susceptibility_B": paramList[11:15]}, index = sub_iter, calibration = True)
    vaccineCoverage = s.compute_typical_vaccination(vacDoses)
    vacsUsedTypical, vacsUsedUniversal = s.simulateWithVaccine(vaccineCoverage, vacEfficacy, vacDoses)
    incidenceL, incidenceH, infections_H1, infections_H3, infections_B, perc_H1, perc_H3, perc_B,hospitalizationsL, hospitalizationsH, deathsL, deathsH = s.calibration_output()
    incidence = sum(incidenceL)/1e6 + sum(incidenceH)/1e6

    incidence_H1_0 = (100*sum(list(infections_H1)[0:2]))/(1.*sum(list(infections_H1)))
    incidence_H1_5 = (100*sum(list(infections_H1)[2:6]))/(1.*sum(list(infections_H1)))
    incidence_H1_25 = (100*sum(list(infections_H1)[6:14]))/(1.*sum(list(infections_H1)))
    incidence_H1_65 = (100*sum(list(infections_H1)[14:]))/(1.*sum(list(infections_H1)))
    
    incidence_H3_0 = (100*sum(list(infections_H3)[0:2]))/(1.*sum(list(infections_H3)))
    incidence_H3_5 = (100*sum(list(infections_H3)[2:6]))/(1.*sum(list(infections_H3)))
    incidence_H3_25 = (100*sum(list(infections_H3)[6:14]))/(1.*sum(list(infections_H3)))
    incidence_H3_65 = (100*sum(list(infections_H3)[14:]))/(1.*sum(list(infections_H3)))
    
    incidence_B_0 = (100*sum(list(infections_B)[0:2]))/(1.*sum(list(infections_B)))
    incidence_B_5 = (100*sum(list(infections_B)[2:6]))/(1.*sum(list(infections_B)))
    incidence_B_25 = (100*sum(list(infections_B)[6:14]))/(1.*sum(list(infections_B)))
    incidence_B_65 = (100*sum(list(infections_B)[14:]))/(1.*sum(list(infections_B)))
    return incidence, perc_H1, perc_H3, perc_B, incidence_H1_0, incidence_H1_5, incidence_H1_25, incidence_H1_65, incidence_H3_0, incidence_H3_5, incidence_H3_25, incidence_H3_65, incidence_B_0, incidence_B_5, incidence_B_25, incidence_B_65

#####################################################################
def optimize(cons1, cons2, data_incidence, beta_dict_opt, sub_range):
    
    beta_H1 = []
    beta_H3 = []
    beta_B=[]
    H1_0 = []
    H1_5 = []
    H1_25 = []
    H1_65 = []
    H3_0 = []
    H3_5 = []
    H3_25 = []
    H3_65 = []
    B_0= []
    B_5 = []
    B_25 = []
    B_65 = []
    for sub_iter in range(sub_range[0], sub_range[1]):
        time1 = time.time()
        min_deviation = 1e7
        beta_opt = None
        for trial in xrange(1):
            beta_init = list(np.random.uniform(low=0.0015, high=0.0025, size=3))
            susc_init = list(np.random.uniform(low=0.4, high=1, size=12))
            BL0_incidence = beta_init +susc_init
            bounds = [(0.00001,1)]*15
            beta_opt_raw =  scipy.optimize.minimize(evaluateObjective_incidence, BL0_incidence, args=(data_incidence, sub_iter), method='TNC', jac=None, bounds=bounds, constraints = {cons1, cons2}, tol=None, callback=None, options={'disp': None, 'maxls': 20, 'iprint': -1, 'gtol': 1e-05, 'eps': 1e-08, 'maxiter': 1500000, 'ftol': 2.220446049250313e-09, 'maxcor': 10, 'maxfun': 15000})
            print beta_opt_raw
            deviation = beta_opt_raw['fun']
            print sub_iter, trial, deviation
            if deviation < min_deviation: beta_opt = [num for num in beta_opt_raw['x']]
        print ("run time"), time.time() - time1, deviation, beta_opt
        beta_H1.append(beta_opt[0])
        beta_H3.append(beta_opt[1])
        beta_B.append(beta_opt[2])
        H1_0.append(beta_opt[3])
        H1_5.append(beta_opt[4])
        H1_25.append(beta_opt[5])
        H1_65.append(beta_opt[6])
        H3_0.append(beta_opt[7])
        H3_5.append(beta_opt[8])
        H3_25.append(beta_opt[9])
        H3_65.append(beta_opt[10])
        B_0.append(beta_opt[11])
        B_5.append(beta_opt[12])
        B_25.append(beta_opt[13])
        B_65.append(beta_opt[14])

    return beta_H1, beta_H3, beta_B, H1_0, H1_5, H1_25, H1_65, H3_0, H3_5, H3_25, H3_65, B_0, B_5, B_25, B_65

###########################################################################
def evaluateObjective_incidence(paramList, data_incidence, sub_iter):
	
    time2 = time.time()
    [doses, efficacy, data_incidence, data_perc_H1, data_perc_H3, data_perc_B, data_incidence_H1_0, data_incidence_H1_5, data_incidence_H1_25, data_incidence_H1_65, data_incidence_H3_0, data_incidence_H3_5, data_incidence_H3_25, data_incidence_H3_65, data_incidence_B_0, data_incidence_B_5, data_incidence_B_25, data_incidence_B_65] = data_incidence
	
    vax_incidence, vax_perc_H1, vax_perc_H3, vax_perc_B, vax_incidence_H1_0, vax_incidence_H1_5, vax_incidence_H1_25, vax_incidence_H1_65, vax_incidence_H3_0, vax_incidence_H3_5, vax_incidence_H3_25, vax_incidence_H3_65, vax_incidence_B_0, vax_incidence_B_5, vax_incidence_B_25, vax_incidence_B_65 =run_efficacy_simulation_incidence(paramList, efficacy[sub_iter], doses[sub_iter], sub_iter)
		
	
	## only incidence, H1 and H3 perc and not B because perc B = 100 - (perc_H1+H3)
    deviation = ((vax_incidence - data_incidence[sub_iter])**2)**2 + \
    ((vax_perc_H1 - data_perc_H1[sub_iter])**2)**2 + \
    (vax_perc_H3 - data_perc_H3[sub_iter])**2 + \
    (vax_incidence_H1_0 - data_incidence_H1_0[sub_iter])**2  + \
    (vax_incidence_H1_5 - data_incidence_H1_5[sub_iter])**2 + \
    (vax_incidence_H1_25 - data_incidence_H1_25[sub_iter])**2 +\
    (vax_incidence_H1_65 - data_incidence_H1_65[sub_iter])**2 +\
    (vax_incidence_H3_0 - data_incidence_H3_0[sub_iter])**2  +\
    (vax_incidence_H3_5 - data_incidence_H3_5[sub_iter])**2 +\
    (vax_incidence_H3_25 - data_incidence_H3_25[sub_iter])**2 +\
    (vax_incidence_H3_65 - data_incidence_H3_65[sub_iter])**2 +\
    (vax_incidence_B_0 - data_incidence_B_0[sub_iter])**2  +\
    (vax_incidence_B_5 - data_incidence_B_5[sub_iter])**2 +\
    (vax_incidence_B_25 - data_incidence_B_25[sub_iter])**2 +\
    (vax_incidence_B_65 - data_incidence_B_65[sub_iter])**2 
	
    #print ("betalist"),  deviation, 
    return deviation


######################################################################33

	