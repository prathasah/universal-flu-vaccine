#!/usr/bin/python
import sys
sys.path.insert(0, r'./Influenza')
import numpy as np
import pandas as pd
import Simulation 
from scipy.optimize import fmin_cobyla
from scipy.optimize import basinhopping
###########################################################################

def run_efficacy_simulation(betaList, vacEfficacy, vacDoses, vaccination_type, sub_iter):
	s = Simulation.run_Simulation(paramValues = {"vacEfficacy":vacEfficacy, "betaList": betaList}, index = sub_iter, calibration = True)
	if vaccination_type == "random": vaccineCoverage = s.comupte_random_vaccination(vacDoses)
	elif vaccination_type =="typical": vaccineCoverage = s.compute_typical_vaccination(vacDoses)
	vacsUsedTypical, vacsUsedUniversal = s.simulateWithVaccine(vaccineCoverage, vacEfficacy, vacDoses)
	incidenceL, incidenceH, infections_H1, infections_H3, infections_B, perc_H1, perc_H3, perc_B = s.calibration_output()
	#list_prop_vax_low, list_prop_vax_high, list_total_vax_low, list_total_vax_high = s.vaccinated_output()
	return incidenceL, incidenceH, infections_H1, infections_H3, infections_B, perc_H1, perc_H3, perc_B


###########################################################################
def evaluateObjective(betaList, vacEfficacy, vacDoses, data):
	
	vax_incidence = []
	vax_perc_H1 = []
	vax_perc_H3 = []
	vax_perc_B = []
	for sub_index in xrange(10):
		incidenceL, incidenceH, infections_H1, infections_H3, infections_B, perc_H1, perc_H3, perc_B = run_efficacy_simulation(betaList, vacEfficacy, vacDoses, "typical", sub_index)
		vax_incidence.append(sum(incidenceL)/1e6)
		vax_incidence.append(sum(incidenceH)/1e6)
		vax_perc_H1.append(perc_H1)
		vax_perc_H3.append(perc_H3)
		vax_perc_B.append(perc_B)
		
		
	[data_incidence, data_perc_H1, data_perc_H3, data_perc_B] = data
	## only incidence, H1 and H3 perc and not B because perc B = 100 - (perc_H1+H3)
	deviation = (np.median(vax_incidence) - data_incidence)**2 + (np.median(vax_perc_H1) - data_perc_H1)**2 + (np.median(vax_perc_H3) - data_perc_H3)**2
	
	print ("betalist"), betaList, deviation, np.median(vax_incidence), np.median(vax_perc_H1), np.median(vax_perc_H3), np.median(vax_perc_B), data
	return deviation


######################################################################33
class MyBounds(object):
	def __init__(self, xmax=[1  , 1, 1], xmin=[0.001,0.001,0.001] ):
		self.xmax = np.array(xmax)
		self.xmin = np.array(xmin)
        def __call__(self, **kwargs):
		x = kwargs["x_new"]
		tmax = bool(np.all(x < self.xmax))
		tmin = bool(np.all(x > self.xmin))
		return tmax and tmin


######################################################################33
if __name__ == "__main__":
	
	
	df  = pd.read_csv("calibration_data_median.csv")
	BL0 = [0.04, 0.05, 0.05]
	for index, row in df.iterrows():
		year = row['year']
		total_doses = row['doses']
		efficacy = row['vacEfficacy']/100.
		data_incidence = row['incidence']/1e6
		data_perc_H1 = row['perc_H1']
		data_perc_H3 = row['perc_H3']
		data_perc_B = row['perc_B']
		
		data = [data_incidence, data_perc_H1, data_perc_H3, data_perc_B]
		
		vacDoses = [total_doses, 0]
		vacEfficacy = [efficacy, 0]
		print ("doses, efficacy"), vacDoses, vacEfficacy
		
		cons = (lambda x:  x[0], lambda x:  x[1], lambda x:  x[2])
		#beta_list_opt = fmin_cobyla(evaluateObjective, BL0,  cons, args=  (vacEfficacy, vacDoses, data),  consargs=(),maxfun = 1e10, rhobeg = 0.01,rhoend=1e-7,catol= 0, disp = 0)
		#print ("final beta ================"), beta_list_opt
		# first calibration results [0.0336752  0.04400087 0.04449421]
		
		bounds = [(0.001, 1.), (0.001,1.), (0.001,1.)]
		
			
		minimizer_kwargs = {"method": "TNC", 'bounds': bounds, "args": (vacEfficacy, vacDoses, data,)}
		mybounds = MyBounds()
		beta_list_opt = basinhopping(evaluateObjective, BL0, minimizer_kwargs=minimizer_kwargs, niter=10, stepsize=0.05,accept_test=mybounds, disp = True)
		print ("final beta ================"), beta_list_opt
		#first caibration results [0.05292967 0.05628356 0.06006498]
		incidence, infections_H1, infections_H3, infections_B, perc_H1, perc_H3, perc_B = run_efficacy_simulation(beta_list_opt, vacEfficacy, vacDoses, "typical", 0)
		print ("final values================="), incidence, perc_H1, perc_H3, perc_B
			
		
		
	