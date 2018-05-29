#!/usr/bin/python
import sys
sys.path.insert(0, r'./Influenza')
import numpy as np
import pandas as pd
import Simulation
import csv
import multiprocessing as mp
from scipy.optimize import fmin_cobyla
from scipy.optimize import basinhopping
import time
###########################################################################

def run_efficacy_simulation_incidence(paramList, vacEfficacy, vacDoses, vaccination_type, sub_iter):
	s = Simulation.run_Simulation(paramValues = {"vacEfficacy":vacEfficacy, "betaList": paramList[0:3]}, index = sub_iter, calibration = True)
	vaccineCoverage = s.compute_typical_vaccination(vacDoses)
	vacsUsedTypical, vacsUsedUniversal = s.simulateWithVaccine(vaccineCoverage, vacEfficacy, vacDoses)
	incidenceL, incidenceH, infections_H1, infections_H3, infections_B, perc_H1, perc_H3, perc_B,hospitalizationsL, hospitalizationsH, deathsL, deathsH = s.calibration_output()
	return sum(incidenceL)/1e6 + sum(incidenceH)/1e6, perc_H1, perc_H3, perc_B

###########################################################################
def evaluateObjective_incidence(paramList, vacEfficacy, vacDoses, data):
	
	vax_incidence = []
	vax_perc_H1 = []
	vax_perc_H3 = []
	vax_perc_B = []
	time1 = time.time()
	for sub_index in xrange(1000):
		incidence, perc_H1, perc_H3, perc_B = run_efficacy_simulation_incidence(paramList, vacEfficacy, vacDoses, "typical", sub_index)
		vax_incidence.append(incidence)
		vax_perc_H1.append(perc_H1)
		vax_perc_H3.append(perc_H3)
		vax_perc_B.append(perc_B)
		
	print ("run time"), time.time() - time1	
	[data_incidence, data_perc_H1, data_perc_H3, data_perc_B] = data
	## only incidence, H1 and H3 perc and not B because perc B = 100 - (perc_H1+H3)
	deviation = (np.mean(vax_incidence) - data_incidence)**2 + (np.mean(vax_perc_H1) - data_perc_H1)**2 + (np.mean(vax_perc_H3) - data_perc_H3)**2 
	
	#if (np.median(vax_incidence)- data_incidence >5): deviation = 10**10
	
	print ("betalist"), paramList, deviation
	print ("param fit --"), np.mean(vax_incidence), np.mean(vax_perc_H1), np.mean(vax_perc_H3), np.mean(vax_perc_B)
	print ("data=="), data
	return deviation

###########################################################################

def run_efficacy_simulation_burden(paramList, betalist, vacEfficacy, vacDoses, vaccination_type, sub_iter):
	s = Simulation.run_Simulation(paramValues = {"vacEfficacy":vacEfficacy, "betaList": betalist, "vac_eff_hospitalization": paramList[0], "vac_eff_mortality": paramList[1]}, index = sub_iter, calibration = True)
	vaccineCoverage = s.compute_typical_vaccination(vacDoses)
	vacsUsedTypical, vacsUsedUniversal = s.simulateWithVaccine(vaccineCoverage, vacEfficacy, vacDoses)
	incidenceL, incidenceH, infections_H1, infections_H3, infections_B, perc_H1, perc_H3, perc_B,hospitalizationsL, hospitalizationsH, deathsL, deathsH = s.calibration_output()
	return sum(hospitalizationsL)/1e3 + sum(hospitalizationsH)/1e3, sum(deathsL)/1e3 + sum(deathsH)/1e3

###########################################################################
def evaluateObjective_burden(paramList, vacEfficacy, vacDoses, data):
	
	vax_hosp = []
	vax_mortality =[]
	betalist = data[0:3]
	for sub_index in xrange(1000):
		hospitalizations, deaths = run_efficacy_simulation_burden(paramList, betalist, vacEfficacy, vacDoses, "typical", sub_index)
		vax_hosp.append(hospitalizations)
		vax_mortality.append(deaths)
		
	[data_hopitalization, data_mortality] = data[3:]
	deviation = (np.median(vax_hosp) - data_hospitalization)**2 + (np.median(vax_mortality) - data_mortality)**2
		
	#print ("paramlist"), paramList, deviation
	#print ("fit --"), np.median(vax_hosp), np.median(vax_mortality)
	#print ("data=="), data[3:]
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
	
	#########################
	## create csv for calibrated parameters
	header = ["year", "doses(millions)", "vacEFficacy", "incidence(millions)", "perc_H1", "perc_H3", "perc_B", "hospitalizations(thousands)", "mortality(thousands)", "beta_H1", "beta_H3", "beta_B", "vac_eff_hospitalization", "vac_eff_mortality"]
	writer = csv.writer(open('calibrated_parameters.csv','wb'))
	writer.writerow(header)
	
	df  = pd.read_csv("calibration_data_median.csv")
	BL0_incidence = [0.05, 0.05, 0.04]
	BL0_burden = [0, 0]
	row = df.iloc[[0]]
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
	#beta_list_opt = fmin_cobyla(evaluateObjective, BL0,  cons, args=  (vacEfficacy, vacDoses, data),  consargs=(),maxfun = 1e10, rhobeg = 0.01,rhoend=1e-7,catol= 0, disp = 0)
	#print ("final beta ================"), beta_list_opt
	# first calibration results [0.0336752  0.04400087 0.04449421]
	
	bounds_incidence = [(0.0001, 1.), (0.0001,1.), (0.0001,1.)]
	

	minimizer_kwargs = {"method": "TNC", 'bounds': bounds_incidence, "args": (vacEfficacy, vacDoses, data_incidence,)}
	mybounds = MyBounds_incidence()
	beta_list_opt = basinhopping(evaluateObjective_incidence, BL0_incidence, minimizer_kwargs=minimizer_kwargs, niter=10, stepsize=0.05,accept_test=mybounds, disp = True)
	print ("final beta ================"), beta_list_opt["x"]

	
	bounds_burden = [(0,1), (0,1)]
	data_burden = list(beta_list_opt["x"]) + [data_hospitalization, data_mortality]
	print data_burden
	mybounds = MyBounds_burden()
	minimizer_kwargs = {"method": "TNC", 'bounds': bounds_burden , "args": (vacEfficacy, vacDoses, data_burden,)}
	eff_list_opt = basinhopping(evaluateObjective_burden, BL0_burden, minimizer_kwargs=minimizer_kwargs, niter=2, stepsize=0.05,accept_test=mybounds, disp = True)
	print ("final efficacy params ================"), eff_list_opt["x"]
	
	elements1 = [year, total_doses, efficacy, data_infections, data_perc_H1, data_perc_H3, data_perc_B, data_hospitalization, data_mortality] + list(beta_list_opt['x']) + list(eff_list_opt['x'])
	writer.writerow(elements1)
			
		
		
	