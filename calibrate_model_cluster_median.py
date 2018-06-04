#!/usr/bin/python
import sys
sys.path.insert(0, r'./Influenza')
import numpy as np
np.warnings.filterwarnings('ignore')
import pandas as pd
import Simulation
import csv
import multiprocessing as mp
from multiprocessing import Pool
import itertools
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

def run_efficacy_simulation_incidence(paramList, vacEfficacy, vacDoses, vaccination_type, sub_range):
	
	incidence_list = []
	perc_h1_list =[]
	perc_h3_list = []
	perc_b_list = []

	for sub_iter in range(sub_range[0], sub_range[1]):
		#print sub_iter
		s = Simulation.run_Simulation(paramValues = {"vacEfficacy":vacEfficacy, "betaList": paramList[0:3]}, index = sub_iter, calibration = True)
		vaccineCoverage = s.compute_typical_vaccination(vacDoses)
		vacsUsedTypical, vacsUsedUniversal = s.simulateWithVaccine(vaccineCoverage, vacEfficacy, vacDoses)
		incidenceL, incidenceH, infections_H1, infections_H3, infections_B, perc_H1, perc_H3, perc_B,hospitalizationsL, hospitalizationsH, deathsL, deathsH = s.calibration_output()
	
		incidence_list.append(sum(incidenceL)/1e6 + sum(incidenceH)/1e6)
		perc_h1_list.append(perc_H1)
		perc_h3_list.append(perc_H3)
		perc_b_list.append(perc_B)
		
	
	return incidence_list, perc_h1_list, perc_h3_list, perc_b_list


	
###########################################################################

def run_efficacy_simulation_burden(paramList, betalist, vacEfficacy, vacDoses, vaccination_type, sub_range):
	
	hosp_list = []
	death_list =[]
	for sub_iter in range(sub_range[0], sub_range[1]):
		s = Simulation.run_Simulation(paramValues = {"vacEfficacy":vacEfficacy, "betaList": betalist, "vac_eff_hospitalization": paramList[0], "vac_eff_mortality": paramList[1]}, index = sub_iter, calibration = True)
		vaccineCoverage = s.compute_typical_vaccination(vacDoses)
		vacsUsedTypical, vacsUsedUniversal = s.simulateWithVaccine(vaccineCoverage, vacEfficacy, vacDoses)
		incidenceL, incidenceH, infections_H1, infections_H3, infections_B, perc_H1, perc_H3, perc_B,hospitalizationsL, hospitalizationsH, deathsL, deathsH = s.calibration_output()
		hosp_list.append(sum(hospitalizationsL)/1e3 + sum(hospitalizationsH)/1e3)
		death_list.append(sum(deathsL)/1e3 + sum(deathsH)/1e3)
	
	
	return hosp_list, death_list


###########################################################################
def evaluateObjective_incidence(paramList, vacEfficacy, vacDoses, data):
	
	time1 = time.time()
	processes = 15
	pool = Pool()
	chunksize = 1000/processes
	sub_indexes = [((num*chunksize, num*chunksize+ chunksize)) for num in xrange(processes)]
	res = pool.map_async(universal_worker, pool_args(run_efficacy_simulation_incidence,paramList, vacEfficacy, vacDoses, "typical", sub_indexes))
	pool.close()
	pool.join()
	results = res.get()
	
	vax_incidence = [result[0] for result in results]
	vax_perc_H1 = [result[1] for result in results]
	vax_perc_H3 = [result[2] for result in results]
	vax_perc_B = [result[3] for result in results]
	
	vax_incidence = [item for sublist in vax_incidence for item in sublist]
	vax_perc_H1 = [item for sublist in vax_perc_H1 for item in sublist]
	vax_perc_H3 = [item for sublist in vax_perc_H3 for item in sublist]
	vax_perc_B = [item for sublist in vax_perc_B for item in sublist]
	print ("run time"), time.time() - time1
		
	[data_incidence, data_perc_H1, data_perc_H3, data_perc_B] = data
	## only incidence, H1 and H3 perc and not B because perc B = 100 - (perc_H1+H3)
	deviation = (np.median(vax_incidence) - data_incidence)**2 + (np.median(vax_perc_H1) - data_perc_H1)**2 + (np.median(vax_perc_H3) - data_perc_H3)**2
	
	if (np.median(vax_incidence) <5): deviation+= 5000
	if (np.median(vax_perc_H1) < 5): deviation+= 5000
	if (np.median(vax_perc_H3) < 5):deviation+=5000
	if (np.median(vax_perc_B) < 5): deviation+=5000
	
	print ("betalist"), paramList, deviation
	print ("param fit --"), np.median(vax_incidence), np.median(vax_perc_H1), np.median(vax_perc_H3), 100 - np.median(vax_perc_H1) - np.median(vax_perc_H3)
	print ("data=="), data
	return deviation

###########################################################################
def evaluateObjective_burden(paramList, vacEfficacy, vacDoses, data):
	
	betalist = data[0:3]
	processes = 15
	pool = Pool()
	chunksize = 1000/processes
	sub_indexes = [((num*chunksize, num*chunksize+ chunksize)) for num in xrange(processes)]
	time1 = time.time()
	res = pool.map_async(universal_worker, pool_args(run_efficacy_simulation_burden, paramList, betalist, vacEfficacy, vacDoses, "typical", sub_indexes))
	pool.close()
	pool.join()
	results = res.get()
	vax_hosp = [result[0] for result in results]
	vax_mortality = [result[1] for result in results]
	
	vax_hosp = [item for sublist in vax_hosp for item in sublist]
	vax_mortality = [item for sublist in vax_mortality for item in sublist]
		
	[data_hopitalization, data_mortality] = data[3:]
	deviation = (np.median(vax_hosp) - data_hospitalization)**2 + (np.median(vax_mortality) - data_mortality)**2
		
	return deviation
######################################################################33
class MyBounds_incidence(object):
	def __init__(self, xmax=[1  , 1, 1], xmin=[0.00001,0.00001,0.00001] ):
		self.xmax = np.array(xmax)
		self.xmin = np.array(xmin)
        def __call__(self, **kwargs):
		x = kwargs["x_new"]
		tmax = bool(np.all(x < self.xmax))
		tmin = bool(np.all(x > self.xmin))
		return tmax and tmin
######################################################################33
class MyBounds_burden(object):
	def __init__(self, xmax=[1,1], xmin=[0,0] ):
		self.xmax = np.array(xmax)
		self.xmin = np.array(xmin)
        def __call__(self, **kwargs):
		x = kwargs["x_new"]
		tmax = bool(np.all(x < self.xmax))
		tmin = bool(np.all(x > self.xmin))
		return tmax and tmin

######################################################################33
if __name__ == "__main__":
	
	writer = {}
	elements = {}
	
	df  = pd.read_csv("calibration_data_median.csv")
	
	#########################
	## create csv for calibrated parameters
	header = ["year", "doses(millions)", "vacEFficacy", "incidence(millions)", "perc_H1", "perc_H3", "perc_B", "hospitalizations(thousands)", "mortality(thousands)", "beta_H1", "beta_H3", "beta_B", "vac_eff_hospitalization", "vac_eff_mortality"]
	writer = csv.writer(open('results_calibrated_parameters_median.csv','wb'))
	writer.writerow(header)
	
	
	BL0_incidence = [0.02245943,  0.05793087,  0.05943181]
	BL0_burden = [0, 0]
	
	year = df.at[0,'year']
	total_doses = df.at[0,'doses']
	efficacy = df.at[0,'vacEfficacy']/100.
	data_infections = df.at[0,'incidence']/1e6
	data_hospitalization = df.at[0,'hospitalizations']/1e3
	data_mortality = df.at[0,'mortality']/1e3
	data_perc_H1 = df.at[0,'perc_H1']
	data_perc_H3 = df.at[0,'perc_H3']
	data_perc_B = df.at[0,'perc_B']
	
	data_incidence = [data_infections, data_perc_H1, data_perc_H3, data_perc_B]
	
	vacDoses = [total_doses, 0]
	vacEfficacy = [efficacy, 0]
	print ("doses, efficacy"), vacDoses, vacEfficacy
	
	cons = (lambda x:  x[0], lambda x:  x[1], lambda x:  x[2])		
	bounds_incidence = [(0.00001, 1.), (0.00001,1.), (0.00001,1.)]

	minimizer_kwargs = {"method": "TNC", 'bounds': bounds_incidence, "args": (vacEfficacy, vacDoses, data_incidence,)}
	mybounds = MyBounds_incidence()
	beta_list_opt = basinhopping(evaluateObjective_incidence, BL0_incidence, minimizer_kwargs=minimizer_kwargs, niter=10, stepsize=0.05, disp = True,niter_success=3)
	print ("final beta ================"), beta_list_opt["x"]

	
	bounds_burden = [(0,1), (0,1)]
	data_burden = list(beta_list_opt["x"]) + [data_hospitalization, data_mortality]
	print data_burden
	mybounds = MyBounds_burden()
	minimizer_kwargs = {"method": "TNC", 'bounds': bounds_burden , "args": (vacEfficacy, vacDoses, data_burden,)}
	eff_list_opt = basinhopping(evaluateObjective_burden, BL0_burden, minimizer_kwargs=minimizer_kwargs, niter=2, stepsize=0.05,accept_test=mybounds, disp = True,niter_success=3)
	print ("final efficacy params ================"), eff_list_opt["x"]
	
	elements = [year, total_doses, efficacy, data_infections, data_perc_H1, data_perc_H3, data_perc_B, data_hospitalization, data_mortality] + list(beta_list_opt['x']) + list(eff_list_opt['x'])
	writer.writerow(elements)
			
		
		
	