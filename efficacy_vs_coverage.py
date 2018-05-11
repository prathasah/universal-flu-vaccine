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
	infectionsL, infectionsH, hospitalizationsL, hospitalizationsH,mortalityL, mortalityH = s.short_output()
	prop_vax_TL, prop_vax_TH, prop_vax_NL, prop_vax_NH, doses_TL, doses_TH, doses_NL, doses_NH = s.vaccinated_output()
	return infectionsL, infectionsH, hospitalizationsL, hospitalizationsH,mortalityL, mortalityH, prop_vax_TL, prop_vax_TH, prop_vax_NL, prop_vax_NL, prop_vax_NH, doses_TL, doses_TH, doses_NL, doses_NH
		

######################################################################33
if __name__ == "__main__":
	
	
	incidence = []
	hosp = []
	total_doses = 141.35e6
	mort = []
	for sub_index in xrange(100):
		print sub_index
		scenario = 0
		dose2 = scenario * total_doses
		dose1 = total_doses - dose2
		doses = [dose1, dose2]
		#print ("doses"), doses, total_doses
		infectionsL, infectionsH, hospitalizationsL, hospitalizationsH,mortalityL, mortalityH, prop_vax_TL, prop_vax_TH, prop_vax_NL, prop_vax_NL, prop_vax_NH, doses_TL, doses_TH, doses_NL, doses_NH = run_efficacy_simulation([0.48,0], doses, "typical", sub_index)
		incidence.append(sum(infectionsL))
		incidence.append(sum(infectionsH))
		hosp.append(sum(hospitalizationsL))
		hosp.append(sum(hospitalizationsH))
		mort.append(sum(mortalityL))
		mort.append(sum(mortalityH))
		#print [round(a/1e3,2) + round(b/1e3,2) for (a,b) in zip(infectionsL, infectionsH)]
		#print ("----->"), sum(infectionsL)/1e6 + sum(infectionsH)/1e6
		
	print ("infections"), np.median(incidence)/1e6
	print ("hospitalizations"), np.median(hosp)/1e3
	print ("mortality"), np.median(mort)/1e3
	
