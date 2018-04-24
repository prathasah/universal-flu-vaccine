#!/usr/bin/python
import sys
sys.path.insert(0, r'./Influenza')
import numpy as np
import Simulation 
## typical vaccination proportions: [0.7012, 0.5514, 0.3264, 0.4528, 0.6532]
##for age groups: 6m-4y, 5y - 17y,18y-49y,50y - 64y,65y+
##typical_coverage = [0.7012, 0.5514, 0.3264, 0.4528, 0.6532]


def run_efficacy_simulation(vacEfficacy, vacDoses, vaccination_type, sub_iter):
	s = Simulation.run_Simulation(paramValues = {"vacEfficacy":vacEfficacy}, index = sub_iter)
	if vaccination_type == "random": vaccineCoverage = s.comupte_random_vaccination(doses)
	elif vaccination_type =="typical": vaccineCoverage = s.compute_typical_vaccination(doses)
	vacsUsed = s.simulateWithVaccine(vaccineCoverage, vacEfficacy, vacDoses)
	infections, hospitalizations, mortality,DALY = s.short_output()
	#check = s.debug_info()
	list_prop_vax_low, list_prop_vax_high, list_total_vax_low, list_total_vax_high = s.vaccinated_output()
	return list_prop_vax_low, list_prop_vax_high, list_total_vax_low, list_total_vax_high, infections, hospitalizations, mortality, DALY
		

######################################################################33
if __name__ == "__main__":
	
	
	vax_60_incidence = []
	total_doses = 150e6
	for sub_index in xrange(1):
		scenario = 0
		dose2 = scenario * total_doses
		dose1 = total_doses - dose2
		doses = [dose1, dose2]
		print ("doses"), doses, total_doses
		list_prop_vax_low, list_prop_vax_high, list_total_vax_low, list_total_vax_high, infections, hospitalizations, mortality, DALY = run_efficacy_simulation([0,0], doses, "typical", sub_index)
		vax_60_incidence.append(sum(infections)/1e6)
		
	print ("typical vaccination incidence"), np.mean(vax_60_incidence), np.std(vax_60_incidence)
	