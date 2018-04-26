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
	vacsUsedTypical, vacsUsedUniversal = s.simulateWithVaccine(vaccineCoverage, vacEfficacy, vacDoses)
	#infections, hospitalizations, mortality,DALY = s.short_output()
	incidence, infections_H1, infections_H3, infections_B, perc_H1, perc_H3, perc_B = s.calibration_output()
	list_prop_vax_low, list_prop_vax_high, list_total_vax_low, list_total_vax_high = s.vaccinated_output()
	return vacsUsedTypical, vacsUsedUniversal, incidence, infections_H1, infections_H3, infections_B, perc_H1, perc_H3, perc_B
		

######################################################################33
if __name__ == "__main__":
	
	
	vax_incidence = []
	vax_perc_H1 = []
	vax_perc_H3 = []
	vax_perc_B = []
	incidence_H1 =[]
	incidence_H3 = []
	incidence_B = []
	total_doses = 141.35e6
	for sub_index in xrange(1000):
		print sub_index
		scenario = 0
		dose2 = scenario * total_doses
		dose1 = total_doses - dose2
		doses = [dose1, dose2]
		#print ("doses"), doses, total_doses
		vacsUsedTypical, vacsUsedUniversal, incidence, infections_H1, infections_H3, infections_B, perc_H1, perc_H3, perc_B = run_efficacy_simulation([0.48,0], doses, "typical", sub_index)
		vax_incidence.append(incidence)
		vax_perc_H1.append(perc_H1)
		vax_perc_H3.append(perc_H3)
		vax_perc_B.append(perc_B)
		incidence_H1.append(infections_H1)
		incidence_H3.append(infections_H3)
		incidence_B.append(infections_B)
	
	print ("total vaccines used"), vacsUsedTypical/1e6, vacsUsedUniversal/1e6	
	print ("typical vaccination incidence"), np.mean(vax_incidence)/1e6, np.std(vax_incidence)
	print ("percentage H1 infections"), np.mean(vax_perc_H1), np.mean(incidence_H1)
	print ("percentage H3 infections"), np.mean(vax_perc_H3), np.mean(incidence_H3)
	print ("percentage B infections"), np.mean(vax_perc_B), np.mean(incidence_B)
	