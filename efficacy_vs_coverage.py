#!/usr/bin/python
import sys
sys.path.insert(0, r'./Influenza')
import numpy as np
import Simulation 

def run_efficacy_simulation(vacEfficacy, vacDoses, sub_iter):
	s = Simulation.run_Simulation(paramValues = {"vacEfficacy":vacEfficacy}, index = sub_iter)
	vaccineCoverage = s.compute_typical_vaccination(vacDoses)
	vacsUsedTypical, vacsUsedUniversal = s.simulateWithVaccine(vaccineCoverage, vacEfficacy, vacDoses)
	#print ("check seasonal, universal doses "), vacsUsedTypical, vacsUsedUniversal
	incidenceL, incidenceH, infections_H1, infections_H3, infections_B, perc_H1, perc_H3, perc_B,hospitalizationsL, hospitalizationsH, deathsL, deathsH = s.calibration_output()
	return incidenceL, incidenceH, infections_H1, infections_H3, infections_B, hospitalizationsL, hospitalizationsH, deathsL, deathsH
		

######################################################################33
if __name__ == "__main__":
	
	
	incidence = []
	hosp = []
	total_doses = 141950000.0
	mort = []
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
	
	for sub_index in xrange(5000):
		
		scenario = 0
		dose2 = scenario * total_doses
		dose1 = total_doses - dose2
		doses = [dose1, dose2]
		#print ("doses"), doses, total_doses
		incidenceL, incidenceH, infections_H1, infections_H3, infections_B, hospitalizationsL, hospitalizationsH, deathsL, deathsH = run_efficacy_simulation([0.4583,0], doses, sub_index)
		incidence.append(sum(incidenceL) + sum(incidenceH))
		
		hosp.append(sum(hospitalizationsL) + sum(hospitalizationsH))
		mort.append(sum(deathsL) + sum(deathsH))
		print sub_index
		incidence_raw = sum(incidenceL) + sum(incidenceH)
		incidence_H1_0.append((100*sum(list(infections_H1)[0:2]))/(1.*incidence_raw))
                incidence_H1_5.append((100*sum(list(infections_H1)[2:6]))/(1.*incidence_raw))
                incidence_H1_25.append((100*sum(list(infections_H1)[6:14]))/(1.*incidence_raw))
                incidence_H1_65.append((100*sum(list(infections_H1)[14:]))/(1.*incidence_raw))
                        
                incidence_H3_0.append((100*sum(list(infections_H3)[0:2]))/(1.*incidence_raw))
                incidence_H3_5.append((100*sum(list(infections_H3)[2:6]))/(1.*incidence_raw))
                incidence_H3_25.append((100*sum(list(infections_H3)[6:14]))/(1.*incidence_raw))
                incidence_H3_65.append((100*sum(list(infections_H3)[14:]))/(1.*incidence_raw))
                        
                incidence_B_0.append((100*sum(list(infections_B)[0:2]))/(1.*incidence_raw))
                incidence_B_5.append((100*sum(list(infections_B)[2:6]))/(1.*incidence_raw))
                incidence_B_25.append((100*sum(list(infections_B)[6:14]))/(1.*incidence_raw))
                incidence_B_65.append((100*sum(list(infections_B)[14:]))/(1.*incidence_raw))
        #print ("----->"), incidence_H3_5
		
	print ("infections"), np.mean(incidence)/1e6, np.std(incidence)/1e6
	print ("hospitalizations"), np.mean(hosp)/1e3
	print ("mortality"), np.mean(mort)/1e3
	print ("age specific H1 infections"), np.mean(incidence_H1_0), np.mean(incidence_H1_5), np.mean(incidence_H1_25),np.mean(incidence_H1_65)
	print ("age specific H3 infections"), np.mean(incidence_H3_0), np.mean(incidence_H3_5), np.mean(incidence_H3_25),np.mean(incidence_H3_65)
	print ("age specific B infections"), np.mean(incidence_B_0), np.mean(incidence_B_5), np.mean(incidence_B_25),np.mean(incidence_B_65)
	print ("std H1 infections"), np.std(incidence_H1_0), np.std(incidence_H1_5), np.std(incidence_H1_25),np.std(incidence_H1_65)
	print ("std H3 infections"), np.std(incidence_H3_0), np.std(incidence_H3_5), np.std(incidence_H3_25),np.std(incidence_H3_65)
	print ("std B infections"), np.std(incidence_B_0), np.std(incidence_B_5), np.std(incidence_B_25),np.std(incidence_B_65)
	
	
