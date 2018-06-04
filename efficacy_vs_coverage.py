#!/usr/bin/python
import sys
sys.path.insert(0, r'./Influenza')
import numpy as np
import Simulation 

def run_efficacy_simulation(vacEfficacy, vacDoses, vaccination_type, sub_iter):
	s = Simulation.run_Simulation(paramValues = {"vacEfficacy":vacEfficacy}, index = sub_iter)
	if vaccination_type == "random": vaccineCoverage = s.comupte_random_vaccination(vacDoses)
	elif vaccination_type =="typical": vaccineCoverage = s.compute_typical_vaccination(vacDoses)
	vacsUsedTypical, vacsUsedUniversal = s.simulateWithVaccine(vaccineCoverage, vacEfficacy, vacDoses)
	incidenceL, incidenceH, infections_H1, infections_H3, infections_B, perc_H1, perc_H3, perc_B,hospitalizationsL, hospitalizationsH, deathsL, deathsH = s.calibration_output()
	return incidenceL, incidenceH, infections_H1, infections_H3, infections_B, perc_H1, perc_H3, perc_B,hospitalizationsL, hospitalizationsH, deathsL, deathsH
		

######################################################################33
if __name__ == "__main__":
	
	
	incidence = []
	hosp = []
	total_doses = 141.35e6
	mort = []
	perc_h1 =[]
	perc_h3 = []
	perc_b = []
	for sub_index in xrange(1000):
		print sub_index
		scenario = 0
		dose2 = scenario * total_doses
		dose1 = total_doses - dose2
		doses = [dose1, dose2]
		#print ("doses"), doses, total_doses
		incidenceL, incidenceH, infections_H1, infections_H3, infections_B, perc_H1, perc_H3, perc_B,hospitalizationsL, hospitalizationsH, deathsL, deathsH = run_efficacy_simulation([0.48,0], doses, "typical", sub_index)
		incidence.append(sum(incidenceL) + sum(incidenceH))
		hosp.append(sum(hospitalizationsL) + sum(hospitalizationsH))
		mort.append(sum(deathsL) + sum(deathsH))
		perc_h1.append(perc_H1)
		perc_h3.append(perc_H3)
		perc_b.append(perc_B)
		
		#print ("----->"), sum(incidenceL)/1e6 + sum(incidenceH)/1e6, sum(hospitalizationsL)/1e3 + sum(hospitalizationsH)/1e3
		
	print ("infections"), np.mean(incidence)/1e6, np.median(incidence)/1e3
	print ("proportion of H1, h3, B"), np.mean(perc_h1), np.mean(perc_h3), np.mean(perc_b)
	print ("hospitalizations"), np.mean(hosp)/1e3, np.median(hosp)/1e3
	print ("mortality"), np.mean(mort)/1e3, np.median(mort)/1e3
	
	
