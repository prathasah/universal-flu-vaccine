import csv
import numpy


######################################################################################		
if __name__ == "__main__":
    
    
    header = ["iter", "recovery_rate", "proportionHighRisk_0", "proportionHighRisk_2","proportionHighRisk_5","proportionHighRisk_19", "proportionHighRisk_25", "proportionHighRisk_50","proportionHighRisk_65", "susceptibility_0", "susceptibility_5", "susceptibility_25", "susceptibility_50", "susceptibility_65", "relative_vaccineEfficacyVsInfection_0", "relative_vaccineEfficacyVsInfection_0.5","relative_vaccineEfficacyVsInfection_16","relative_vaccineEfficacyVsInfection_65","vaccineEfficacyVsDeath_0", "vaccineEfficacyVsDeath_20","vaccineEfficacyVsDeath_65", "vaccineEfficacyVsHospitalization_0", "vaccineEfficacyVsHospitalization_5", "vaccineEfficacyVsHospitalization_16", "vaccineEfficacyVsHospitalization_65","caseMortality_0","caseMortality_5","caseMortality_18","caseMortality_50","caseMortality_65",  "caseHospitalization_0","caseHospitalization_5","caseHospitalization_18","caseHospitalization_50","caseHospitalization_65", "highRiskRelativeCaseMortality_0" ,"highRiskRelativeCaseMortality_20" ,"highRiskRelativeCaseMortality_65" ,"R0"]
    writer = csv.writer(open('sampled_parameter_set.csv','wb'))
    writer.writerow(header)
    
    for num in xrange(1000):
        rec_rate = 1. / numpy.random.triangular(0.9, 2.5, 4.6)
        prop_high_risk_0 = numpy.random.normal(0.0415, 0.0044)
        prop_high_risk_2 = numpy.random.normal(0.0883, 0.0051)
        prop_high_risk_5 = numpy.random.normal(0.1168, 0.0030)
        prop_high_risk_19 = numpy.random.normal(0.1235, 0.0055)
        prop_high_risk_25 = numpy.random.normal(0.1570, 0.0027)
        prop_high_risk_50 = numpy.random.normal(0.3056, 0.0044)
        prop_high_risk_65 = numpy.random.normal(0.4701, 0.0050)
        susc_0 = numpy.random.uniform(0.9, 1.)
        susc_5 =      numpy.random.uniform(0.7, 1.)
        susc_25 =     numpy.random.uniform(0.6, 1.)
        susc_50 =     numpy.random.uniform(0.5, 1.)
        susc_65 =      numpy.random.uniform(0.3, 1)
        vac_eff_inf_0 = 0
        vac_eff_inf_0_5 = numpy.random.triangular(0.4, 0.7, 1.)
        vac_eff_inf_16 = numpy.random.triangular(0.5, 0.7, 0.9)
        vac_eff_inf_65 = numpy.random.triangular(0.4, 0.5, 0.6)
        vac_eff_death_0 = numpy.random.uniform(0.4, 0.75)
        vac_eff_death_20 = numpy.random.uniform(0.4, 0.7)
        vac_eff_death_65 = numpy.random.uniform(0.3, 0.7)
        vac_eff_hosp_0 = numpy.random.triangular(0.44, 0.6, 0.72)
        vac_eff_hosp_5 = numpy.random.triangular(0.35, 0.53,0.66)
        vac_eff_hosp_16 = numpy.random.triangular(0.24, 0.43, 0.57)
        vac_eff_hosp_65 = numpy.random.triangular(0.37, 0.50, 0.60)
        case_mort_0 = numpy.random.triangular(0.00002, 0.00004, 0.00026)
        case_mort_5 = numpy.random.triangular(0.00001, 0.00002, 0.00010)
        case_mort_18 = numpy.random.triangular(0.00009, 0.00010, 0.00159)
        case_mort_50 = numpy.random.triangular(0.00010, 0.00134, 0.00159)
        case_mort_65 = numpy.random.triangular(0.00010, 0.00090, 0.01170)
        case_hosp_0 = numpy.random.triangular(0.0033, 0.0141, 0.0245)
        case_hosp_5 = numpy.random.triangular(0.0006, 0.0011, 0.0061)
        case_hosp_18 = numpy.random.triangular(0.0015, 0.0042, 0.0300)
        case_hosp_50 = numpy.random.triangular(0.0015, 0.0193, 0.0300)
        case_hosp_65 = numpy.random.triangular(0.0016, 0.0184, 0.0421)
        ##check the figures below
        high_risk_mortality_0 = numpy.random.triangular(0.4, 0.6, 21.9) / numpy.random.triangular(0.041, 0.07, 0.30)
        high_risk_mortality_20 =  numpy.random.uniform(0.8, 24.9) / numpy.random.triangular(0.21, 0.31, 0.41),
        high_risk_mortality_65 = numpy.random.uniform(23, 29.6) / numpy.random.triangular(2.3, 3.51, 4.52)
        
        R0 = 1.2375
        
        elements = [num, rec_rate,  susc_0, susc_5, susc_25, susc_50, susc_65, vac_eff_inf_0, vac_eff_inf_0_5, vac_eff_inf_16, vac_eff_inf_65, vac_eff_death_0, vac_eff_death_20, vac_eff_death_65, vac_eff_hosp_0, vac_eff_hosp_5, vac_eff_hosp_16, vac_eff_hosp_65, case_mort_0, case_mort_5, case_mort_18, case_mort_50, case_hosp_65, case_hosp_0, case_hosp_5, case_hosp_18, case_hosp_50, case_hosp_65, high_risk_mortality_0, high_risk_mortality_20, high_risk_mortality_65, R0]
        writer.writerow(elements)
        
        


        
        
        
