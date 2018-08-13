#!/usr/bin/python
import sys
sys.path.insert(0, r'../../Influenza')
import numpy as np
np.warnings.filterwarnings('ignore')
import pandas as pd
import Simulation
import csv
import mpi4py.MPI
from multiprocessing import Pool
import itertools
import scipy.optimize
from scipy.optimize import fmin_cobyla
from scipy.optimize import basinhopping
import time
from copy import deepcopy
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

def run_efficacy_simulation_burden(paramList, betalist, susc_list, efficacy, doses, sub_iter):
	
	
    vacDoses = [doses,0]
    vacEfficacy = [efficacy,0]

    s = Simulation.run_Simulation(paramValues = {"vacEfficacy":vacEfficacy, "betaList": betalist, "susceptibility_H1": susc_list[0:4], "susceptibility_H3": susc_list[4:8], "susceptibility_B": susc_list[8:12], "vac_eff_hospitalization": paramList[0], "vac_eff_mortality": paramList[1], "prob_hosp_scaling": paramList[2], "prob_death_scaling": paramList[3]}, index = sub_iter, calibration = True)
    vaccineCoverage = s.compute_typical_vaccination(vacDoses)
    vacsUsedTypical, vacsUsedUniversal = s.simulateWithVaccine(vaccineCoverage, vacEfficacy, vacDoses)
    incidenceL, incidenceH, infections_H1, infections_H3, infections_B, perc_H1, perc_H3, perc_B,hospitalizationsL, hospitalizationsH, deathsL, deathsH = s.calibration_output()
    
    hosp  = sum(hospitalizationsL)/1e3 + sum(hospitalizationsH)/1e3
    deaths = sum(deathsL)/1e3 + sum(deathsH)/1e3
	
	
    return hosp, deaths

#####################################################################
def optimize(cons1, cons2, data_list,  sub_iter):
    
    [iteration, year, doses, efficacy, data_hospitalization, data_mortality,beta_H1, beta_H3, beta_B, susceptibility_H1_0, susceptibility_H1_5,susceptibility_H1_25,susceptibility_H1_65,susceptibility_H3_0, susceptibility_H3_5,susceptibility_H3_25,susceptibility_H3_65,susceptibility_B_0, susceptibility_B_5,susceptibility_B_25,susceptibility_B_65]= data_list

    eff_hosp = []
    eff_death=[]
    hosp_scaling  =[]
    death_scaling = []
    ## optimize efficacy parameters (paramList, beta_list, susc_list, vacEfficacy, vacDoses, data_hospitalization, data_mortality)
    bounds = [(0.00000001,1)]*2 + [(0.00000001,1)]*2
    min_eff_deviation = 1e10
    for trial in xrange(1):
	BL0_burden = list(np.random.uniform(low=0.00000001, high=1, size=2)) + list(np.random.uniform(low=0.00000001, high=0.1, size=2))
	beta_opt =[beta_H1[sub_iter], beta_H3[sub_iter], beta_B[sub_iter], susceptibility_H1_0[sub_iter], susceptibility_H1_5[sub_iter],susceptibility_H1_25[sub_iter],susceptibility_H1_65[sub_iter],susceptibility_H3_0[sub_iter], susceptibility_H3_5[sub_iter],susceptibility_H3_25[sub_iter],susceptibility_H3_65[sub_iter],susceptibility_B_0[sub_iter], susceptibility_B_5[sub_iter],susceptibility_B_25[sub_iter],susceptibility_B_65[sub_iter]]
	eff_opt_raw =  scipy.optimize.minimize(evaluateObjective_burden, BL0_burden, args=(beta_opt[0:3], beta_opt[3:], efficacy[sub_iter], doses[sub_iter], data_hospitalization[sub_iter], data_mortality[sub_iter], sub_iter), method='TNC', jac=None, bounds=bounds, constraints = {cons1, cons2}, tol=None, callback=None, options={'disp': None, 'maxls': 20, 'iprint': -1, 'gtol': 1e-08, 'eps': 1e-08, 'maxiter': 150000000, 'ftol': 2.220446049250313e-09, 'maxcor': 10, 'maxfun': 1500000})
	eff_deviation = eff_opt_raw['fun']
	#print ("optimized  efficacy = "), sub_iter, eff_deviation, trial
	if eff_deviation < min_eff_deviation:
	    min_eff_deviation = deepcopy(eff_deviation)
	    eff_opt = [num for num in eff_opt_raw['x']]
		
    
    print ("final eff =="), sub_iter, min_eff_deviation, eff_opt
    
    eff_hosp.append(eff_opt[0])
    eff_death.append(eff_opt[1])
    hosp_scaling.append(eff_opt[2])
    death_scaling.append(eff_opt[3])
    
    vax_hosp, vax_death = run_efficacy_simulation_burden(eff_opt, beta_opt[0:3], beta_opt[3:], efficacy[sub_iter], doses[sub_iter], sub_iter)
    print ("data vs calibration = "), sub_iter, data_hospitalization[sub_iter], data_mortality[sub_iter], vax_hosp, vax_death
     
    return eff_hosp, eff_death, hosp_scaling, death_scaling


    
###########################################################################
def evaluateObjective_burden(paramList, beta_list, susc_list, vacEfficacy, vacDoses, data_hosp, data_mort, sub_iter):
	
	
	vax_hosp, vax_mortality = run_efficacy_simulation_burden(paramList, beta_list, susc_list, vacEfficacy, vacDoses, sub_iter)
	
	deviation = (vax_hosp - data_hosp)**2 + (vax_mortality - data_mort)**2
	#print ("efficacy param list"), paramList, deviation
	#print ("param fit --"), np.mean(vax_hosp), np.mean(vax_mortality)
	#print ("data"), data_hosp, data_mort
	
	return deviation

###################################################

def cons1(x):
    return 1 - max(x)-1e-08

def cons2(x):
    return min(x)-1e-08


#######################################################################
def create_files(num_files):

	writer = {}
	header = ["iter", "incidence(millions)", "vac_eff_hospitalization", "vac_eff_mortality", "prob_hosp_scaling", "prob_death_scaling"]
	
	for num in xrange(num_files):
		writer[num] = csv.writer(open('efficacy_parameters_num_'+str(num)+'.csv','wb'))
		writer[num].writerow(header)
	
	return writer

##########################################################################
def start_pool(index, mpi_index,mpi_chunksize,data_list, data_infections):
	
	for num in range(mpi_index[0], mpi_index[1]):
		print ("starting pool"), num
		eff_hosp, eff_death,hosp_scaling, death_scaling = optimize(cons1, cons2, data_list,  num)
		elem1 = [num, data_infections[num],eff_hosp, eff_death,hosp_scaling, death_scaling]
		writer[index].writerow(elem1)
	
	
	
######################################################################33
if __name__ == "__main__":
	
	elements = {}
	
	df  = pd.read_csv("sampled_parameter_5000_set_with_data_13July2018.csv")
	year = df['year'].tolist()
	iteration = df['iter'].tolist()
	total_doses = df['data_vacDoses'].tolist()
	efficacy = df['data_vacEfficacy'].tolist()
	efficacy = [num/100. for num in efficacy]
	data_infections = df['data_incidence'].tolist()
	data_infections = [num/1e6 for num in data_infections]
	data_hospitalization = df['data_hospitalizations'].tolist()
	data_hospitalization = [num/1e3 for num in data_hospitalization]
	data_mortality = df['data_mortality'].tolist()
	data_mortality = [num/1e3 for num in data_mortality]
	
	df1 = pd.read_csv("results_calibrated_parameters_COMBINED.csv")
	beta_H1 = df1['beta_H1'].tolist()
	beta_H3 = df1['beta_H3'].tolist()
	beta_B = df1['beta_B'].tolist()
	susceptibility_H1_0 = df1['susceptibility_H1_0'].tolist()
	susceptibility_H1_5 = df1['susceptibility_H1_5'].tolist()
	susceptibility_H1_25 = df1['susceptibility_H1_25'].tolist()
	susceptibility_H1_65= df1['susceptibility_H1_65'].tolist()
	susceptibility_H3_0= df1['susceptibility_H3_0'].tolist()
	susceptibility_H3_5= df1['susceptibility_H3_5'].tolist()
	susceptibility_H3_25 = df1['susceptibility_H3_25'].tolist()
	susceptibility_H3_65 = df1['susceptibility_H3_65'].tolist()
	susceptibility_B_0 = df1['susceptibility_B_0'].tolist()
	susceptibility_B_5 = df1['susceptibility_B_5'].tolist()
	susceptibility_B_25= df1['susceptibility_B_25'].tolist()
	susceptibility_B_65 = df1['susceptibility_B_65'].tolist()

	data_list = [iteration, year, total_doses, efficacy, data_hospitalization, data_mortality,
		     beta_H1, beta_H3, beta_B, susceptibility_H1_0, susceptibility_H1_5,susceptibility_H1_25,susceptibility_H1_65,
		     susceptibility_H3_0, susceptibility_H3_5,susceptibility_H3_25,susceptibility_H3_65,
		     susceptibility_B_0, susceptibility_B_5,susceptibility_B_25,susceptibility_B_65]
	
	rank = mpi4py.MPI.COMM_WORLD.Get_rank()
	size = mpi4py.MPI.COMM_WORLD.Get_size()
	mpi_chunks = 50
	mpi_chunksize = 5000/mpi_chunks
	mpi_indexes = [((num*mpi_chunksize, num*mpi_chunksize+ mpi_chunksize)) for num in xrange(mpi_chunks)]
	writer = create_files(mpi_chunks)
	for (index, mpi_index) in enumerate(mpi_indexes):
		if index%size!=rank: continue
		print "Task number %d being done by processor %d of %d" % (index, rank, size)
		start_pool(index, mpi_index, mpi_chunksize, data_list, data_infections)

	

			
		
		
	