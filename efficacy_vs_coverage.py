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
	#print ("check seasonal, universal doses "), vacsUsedTypical, vacsUsedUniversal
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
	incidence_H1_0 = []
	incidence_H1_5 = []
	incidence_H1_25 = []
	incidence_H1_65 = []
	incidence_H3_0 = []
	incidence_H3_5 = []
	incidence_H3_25 = []
	incidence_H3_65 = []
	incidence_B_0 = []
	incidence_B_5 = []
	incidence_B_25 = []
	incidence_B_65 = []
	
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
		incidence_H1_0.append((100*sum(list(infections_H1)[0:2]))/(1.*sum(list(infections_H1))))
		incidence_H1_5.append((100*sum(list(infections_H1)[2:6]))/(1.*sum(list(infections_H1))))
		incidence_H1_25.append((100*sum(list(infections_H1)[6:14]))/(1.*sum(list(infections_H1))))
		incidence_H1_65.append((100*sum(list(infections_H1)[14:]))/(1.*sum(list(infections_H1))))
		
		incidence_H3_0.append((100*sum(list(infections_H3)[0:2]))/(1.*sum(list(infections_H3))))
		incidence_H3_5.append((100*sum(list(infections_H3)[2:6]))/(1.*sum(list(infections_H3))))
		incidence_H3_25.append((100*sum(list(infections_H3)[6:14]))/(1.*sum(list(infections_H3))))
		incidence_H3_65.append((100*sum(list(infections_H3)[14:]))/(1.*sum(list(infections_H3))))
		
		incidence_B_0.append((100*sum(list(infections_B)[0:2]))/(1.*sum(list(infections_B))))
		incidence_B_5.append((100*sum(list(infections_B)[2:6]))/(1.*sum(list(infections_B))))
		incidence_B_25.append((100*sum(list(infections_B)[6:14]))/(1.*sum(list(infections_B))))
		incidence_B_65.append((100*sum(list(infections_B)[14:]))/(1.*sum(list(infections_B))))
		#print ("----->"), sum(incidenceL)/1e6 + sum(incidenceH)/1e6, sum(hospitalizationsL)/1e3 + sum(hospitalizationsH)/1e3
		
	print ("infections"), np.mean(incidence)/1e6
	print ("proportion of H1, h3, B"), np.mean(perc_h1), np.mean(perc_h3), np.mean(perc_b)
	print ("hospitalizations"), np.mean(hosp)/1e3
	print ("mortality"), np.mean(mort)/1e3
	print ("age specific H1 infections"), np.mean(incidence_H1_0), np.mean(incidence_H1_5), np.mean(incidence_H1_25),np.mean(incidence_H1_65)
	print ("age specific H3 infections"), np.mean(incidence_H3_0), np.mean(incidence_H3_5), np.mean(incidence_H3_25),np.mean(incidence_H3_65)
	print ("age specific B infections"), np.mean(incidence_B_0), np.mean(incidence_B_5), np.mean(incidence_B_25),np.mean(incidence_B_65)
	
	
