import csv
import numpy
import scipy.stats as stats
######################################################################################		
if __name__ == "__main__":
    
    
    header = ["iter", "recovery_rate", "proportionHighRisk_0", "proportionHighRisk_2","proportionHighRisk_5","proportionHighRisk_19", "proportionHighRisk_25", "proportionHighRisk_50","proportionHighRisk_65", "susceptibility_H1_0", "susceptibility_H1_4", "susceptibility_H1_18", "susceptibility_H1_65", "susceptibility_H3_0", "susceptibility_H3_4", "susceptibility_H3_18", "susceptibility_H3_65","susceptibility_B_0", "susceptibility_B_4", "susceptibility_B_18", "susceptibility_B_65",
              "relative_vaccineEfficacyVsInfection_H1_0", "relative_vaccineEfficacyVsInfection_H1_0.5",  "relative_vaccineEfficacyVsInfection_H1_5","relative_vaccineEfficacyVsInfection_H1_18","relative_vaccineEfficacyVsInfection_H1_50",
              "relative_vaccineEfficacyVsInfection_H3_0", "relative_vaccineEfficacyVsInfection_H3_0.5",  "relative_vaccineEfficacyVsInfection_H3_5","relative_vaccineEfficacyVsInfection_H3_18","relative_vaccineEfficacyVsInfection_H3_50",
              "relative_vaccineEfficacyVsInfection_B_0", "relative_vaccineEfficacyVsInfection_B_0.5",  "relative_vaccineEfficacyVsInfection_B_5","relative_vaccineEfficacyVsInfection_B_18","relative_vaccineEfficacyVsInfection_B_50",
               "vaccineEfficacyVsHospitalization_H1_0", "vaccineEfficacyVsHospitalization_H1_0.5", "vaccineEfficacyVsHospitalization_H1_16", "vaccineEfficacyVsHospitalization_H1_65",
              "vaccineEfficacyVsHospitalization_H3_0", "vaccineEfficacyVsHospitalization_H3_0.5", "vaccineEfficacyVsHospitalization_H3_16", "vaccineEfficacyVsHospitalization_H3_65",
              "vaccineEfficacyVsHospitalization_B_0", "vaccineEfficacyVsHospitalization_B_0.5", "vaccineEfficacyVsHospitalization_B_16", "vaccineEfficacyVsHospitalization_B_65",
              "vaccineEfficacyVsDeath_H1_0", "vaccineEfficacyVsDeath_H1_0.5", "vaccineEfficacyVsDeath_H1_18", "vaccineEfficacyVsDeath_H1_65",
              "vaccineEfficacyVsDeath_H3_0", "vaccineEfficacyVsDeath_H3_0.5", "vaccineEfficacyVsDeath_H3_18", "vaccineEfficacyVsDeath_H3_65",
              "vaccineEfficacyVsDeath_B_0", "vaccineEfficacyVsDeath_B_0.5", "vaccineEfficacyVsDeath_B_18", "vaccineEfficacyVsDeath_B_65",
              "highRiskvaccineEfficacyVsDeath_H1_0", "highRiskvaccineEfficacyVsDeath_H1_0.5", "highRiskvaccineEfficacyVsDeath_H1_18", "highRiskvaccineEfficacyVsDeath_H1_65",
              "highRiskvaccineEfficacyVsDeath_H3_0", "highRiskvaccineEfficacyVsDeath_H3_0.5", "highRiskvaccineEfficacyVsDeath_H3_18", "highRiskvaccineEfficacyVsDeath_H3_65",
              "highRiskvaccineEfficacyVsDeath_B_0", "highRiskvaccineEfficacyVsDeath_B_0.5", "highRiskvaccineEfficacyVsDeath_B_18", "highRiskvaccineEfficacyVsDeath_B_65",

              "lowRiskcaseMortality_H1_0","lowRiskcaseMortality_H1_5","lowRiskcaseMortality_H1_18","lowRiskcaseMortality_H1_50","lowRiskcaseMortality_H1_65", "lowRiskcaseMortality_H1_75",
              "lowRiskcaseMortality_H3_0","lowRiskcaseMortality_H3_5","lowRiskcaseMortality_H3_18","lowRiskcaseMortality_H3_50","lowRiskcaseMortality_H3_65", "lowRiskcaseMortality_H3_75",
              "lowRiskcaseMortality_B_0","lowRiskcaseMortality_B_5","lowRiskcaseMortality_B_18","lowRiskcaseMortality_B_50","lowRiskcaseMortality_B_65", "lowRiskcaseMortality_B_75",
              "highRiskcaseMortality_H1_0","highRiskcaseMortality_H1_5","highRiskcaseMortality_H1_18","highRiskcaseMortality_H1_50","highRiskcaseMortality_H1_65", "highRiskcaseMortality_H1_75",
              "highRiskcaseMortality_H3_0","highRiskcaseMortality_H3_5","highRiskcaseMortality_H3_18","highRiskcaseMortality_H3_50","highRiskcaseMortality_H3_65", "highRiskcaseMortality_H3_75",
              "highRiskcaseMortality_B_0","highRiskcaseMortality_B_5","highRiskcaseMortality_B_18","highRiskcaseMortality_B_50","highRiskcaseMortality_B_65", "highRiskcaseMortality_B_75",
              
              "lowRiskcaseHospitalization_H1_0","lowRiskcaseHospitalization_H1_5","lowRiskcaseHospitalization_H1_18","lowRiskcaseHospitalization_H1_50","lowRiskcaseHospitalization_H1_65", "lowRiskcaseHospitalization_H1_75",
              "lowRiskcaseHospitalization_H3_0","lowRiskcaseHospitalization_H3_5","lowRiskcaseHospitalization_H3_18","lowRiskcaseHospitalization_H3_50","lowRiskcaseHospitalization_H3_65", "lowRiskcaseHospitalization_H3_75",
              "lowRiskcaseHospitalization_B_0","lowRiskcaseHospitalization_B_5","lowRiskcaseHospitalization_B_18","lowRiskcaseHospitalization_B_50","lowRiskcaseHospitalization_B_65", "lowRiskcaseHospitalization_B_75",
              "highRiskcaseHospitalization_H1_0","highRiskcaseHospitalization_H1_5","highRiskcaseHospitalization_H1_18","highRiskcaseHospitalization_H1_50","highRiskcaseHospitalization_H1_65", "highRiskcaseHospitalization_H1_75",
              "highRiskcaseHospitalization_H3_0","highRiskcaseHospitalization_H3_5","highRiskcaseHospitalization_H3_18","highRiskcaseHospitalization_H3_50","highRiskcaseHospitalization_H3_65", "highRiskcaseHospitalization_H3_75",
              "highRiskcaseHospitalization_B_0","highRiskcaseHospitalization_B_5","highRiskcaseHospitalization_B_18","highRiskcaseHospitalization_B_50","highRiskcaseHospitalization_B_65", "highRiskcaseHospitalization_B_75", "crossImmunity", "R0"]
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
        
        #####
        ##susceptibility
        ##ref: Datasheet downloaded from Age Dependence and Isotype Specificity of Influenza Virus Hemagglutinin Stalk-Reactive Antibodies in Humans
        ## http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.798.7813&rep=rep1&type=pdf
        susc_H1_0 = stats.truncnorm((0 - 0.961538) / 0.05, (1 - 0.961538) / 0.05, loc= 0.961538, scale=0.05).rvs(1)[0]
        susc_H1_4 = stats.truncnorm((0 - 0.714286) / 0.05, (1 - 0.714286) / 0.05, loc= 0.714286, scale=0.05).rvs(1)[0]
        susc_H1_18 =  stats.truncnorm((0 - 0.42623) / 0.05, (1 - 0.42623) / 0.05, loc= 0.42623, scale=0.05).rvs(1)[0]
        susc_H1_65 =   stats.truncnorm((0 - 0.176471) / 0.05, (1 - 0.176471) / 0.05, loc= 0.176471, scale=0.05).rvs(1)[0]
        susc_H3_0 = stats.truncnorm((0 - 0.730769) / 0.05, (1 - 0.730769) / 0.05, loc= 0.730769, scale=0.05).rvs(1)[0]
        susc_H3_4 =    stats.truncnorm((0 - 0.214286) / 0.05, (1 - 0.214286) / 0.05, loc= 0.214286, scale=0.05).rvs(1)[0]
        susc_H3_18 =  stats.truncnorm((0 - 0.704918) / 0.05, (1 - 0.704918) / 0.05, loc= 0.704918, scale=0.05).rvs(1)[0]
        susc_H3_65 =   stats.truncnorm((0 - 0.235294) / 0.05, (1 -  0.235294) / 0.05, loc=  0.235294, scale=0.05).rvs(1)[0]
        susc_B_0 = stats.truncnorm((0 - 0.884615) / 0.05, (1 - 0.884615) / 0.05, loc= 0.884615, scale=0.05).rvs(1)[0]
        susc_B_4 =     stats.truncnorm((0 - 0.714286) / 0.05, (1 - 0.714286) / 0.05, loc= 0.714286, scale=0.05).rvs(1)[0]
        susc_B_18 =   stats.truncnorm((0 - 0.295082) / 0.05, (1 - 0.295082) / 0.05, loc= 0.295082, scale=0.05).rvs(1)[0]
        susc_B_65 =     stats.truncnorm((0 - 0.235294) / 0.05, (1 - 0.235294) / 0.05, loc= 0.235294, scale=0.05).rvs(1)[0]
        ########
        
        ##################3
        ##vaccine efficacy against infection.
        ##ref Table 2 of Assessment of influenza vaccine effectiveness in a sentinel surveillance network 2010-13, United States
        ##Benjamin J. Cowlinga, Shuo Fenga, Lyn Finelli, Andrea Steffens, Ashley Fowlkes
        
        vac_eff_inf_H1_0 = 0
        vac_eff_inf_H1_6mo = numpy.random.triangular(0.48, 0.7, 0.83)
        vac_eff_inf_H1_5 = numpy.random.triangular(0.31, 0.56, 0.72)
        vac_eff_inf_H1_18 = numpy.random.triangular(0.39, 0.57, 0.70)
        vac_eff_inf_H1_50 = numpy.random.triangular(0.1, 0.3, 0.63)
        
        vac_eff_inf_H3_0 = 0
        vac_eff_inf_H3_6mo = numpy.random.triangular(0.47, 0.60, 0.69)
        vac_eff_inf_H3_5 = numpy.random.triangular(0.16, 0.33, 0.47)
        vac_eff_inf_H3_18 = numpy.random.triangular(0.08, 0.26, 0.41)
        vac_eff_inf_H3_50 = numpy.random.triangular(0.09, 0.34, 0.52)
        
        vac_eff_inf_B_0 = 0
        vac_eff_inf_B_6mo = numpy.random.triangular(0.35, 0.51, 0.64)
        vac_eff_inf_B_5 = numpy.random.triangular(0.4, 0.52, 0.62)
        vac_eff_inf_B_18 = numpy.random.triangular(0.09, 0.33, 0.50)
        vac_eff_inf_B_50 = numpy.random.triangular(0.10, 0.37, 0.65)
        ##################
        
        #############################
        ##vaccine efficacy against hospitalization
        ## for 6mo - 5years: Table 3 of Vaccine effectiveness against laboratoryconfirmed influenza
        ## hospitalizations among young children during the 2010-11 to 2013-14 influenza seasons in Ontario, Canada
        
        ## for 16-64 and 65+ : Table 2 of
        ## Effectiveness of influenza vaccines in preventing severe influenza illness among adults: A systematic review and
        ##meta-analysis of test-negative design case-control studies
        
        vac_eff_hosp_H1_0 = 0
        vac_eff_hosp_H1_6mo = numpy.random.triangular(0.273, 0.821,0.956)
        vac_eff_hosp_H1_16 = numpy.random.triangular(0.34, 0.55, 0.76)
        vac_eff_hosp_H1_65 = numpy.random.triangular(0.26, 0.54, 0.82)
        
        vac_eff_hosp_H3_0 = 0
        vac_eff_hosp_H3_6mo = numpy.random.triangular(0.035, 0.533,0.774)
        vac_eff_hosp_H3_16 = numpy.random.triangular(0.38, 0.50, 0.62)
        vac_eff_hosp_H3_65 = numpy.random.triangular(0.21, 0.33, 0.45)
        
        
        vac_eff_hosp_B_0 = 0
        vac_eff_hosp_B_6mo = numpy.random.triangular(0.283, 0.580,0.754)
        vac_eff_hosp_B_16 = numpy.random.triangular(0.08, 0.45, 0.81)
        vac_eff_hosp_B_65 = numpy.random.triangular(0.11, 0.31, 0.51)
        ###########################
        
        #################################
        ## vaccine efficacy against death
        ##for 6mo-17 years: Table 3 of nfluenza Vaccine Effectiveness Against Pediatric Deaths: 2010-2014; Flanerry et al.
        # for elderly: Table 4 of A Cohort Study of the Effectiveness of Influenza Vaccine in Older People, Performed Using the United Kingdom General Practice Research Database
        # for healthy adults: text in Prioritization of Influenza Pandemic Vaccination to Minimize Years of Life Lost
        # for high risk adults: Table 3 of Clinical Effectiveness of Influenza Vaccination in Persons Younger Than 65 Years With High-Risk Medical Conditions
        
        vac_eff_death_H1_0 = 0
        vac_eff_death_H1_6mo = numpy.random.triangular(0.31, 0.59,0.77)
        vac_eff_death_H1_18 = numpy.random.uniform(0.7, 0.9)
        vac_eff_death_H1_65 = numpy.random.triangular(0.11, 0.21, 0.29)
        
        vac_eff_death_H3_0 = 0
        vac_eff_death_H3_6mo = numpy.random.triangular(0.31, 0.59,0.77)
        vac_eff_death_H3_18 = numpy.random.uniform(0.7, 0.9)
        vac_eff_death_H3_65 = numpy.random.triangular(0.11, 0.21, 0.29)
        
        vac_eff_death_B_0 = 0
        vac_eff_death_B_6mo = numpy.random.triangular(0.43, 0.71,0.87)
        vac_eff_death_B_18 = numpy.random.uniform(0.7, 0.9)
        vac_eff_death_B_65 = numpy.random.triangular(0.11, 0.21, 0.29)
        
        high_risk_vac_eff_death_H1_0 = 0
        high_risk_vac_eff_death_H1_6mo = numpy.random.triangular(0.35, 0.59,0.74)
        high_risk_vac_eff_death_H1_18 = numpy.random.triangular(0.39, 0.78,0.92)
        high_risk_vac_eff_death_H1_65 = numpy.random.triangular(0.22, 0.29, 0.34)
        
        high_risk_vac_eff_death_H3_0 = 0
        high_risk_vac_eff_death_H3_6mo = numpy.random.triangular(0.35, 0.59,0.74)
        high_risk_vac_eff_death_H3_18 = numpy.random.triangular(0.39, 0.78,0.92)
        high_risk_vac_eff_death_H3_65 = numpy.random.triangular(0.22, 0.29, 0.34)
        
        high_risk_vac_eff_death_B_0 = 0
        high_risk_vac_eff_death_B_6mo = numpy.random.triangular(0.13, 0.35,0.63)
        high_risk_vac_eff_death_B_18 = numpy.random.triangular(0.39, 0.78,0.92)
        high_risk_vac_eff_death_B_65 = numpy.random.triangular(0.22, 0.29, 0.34)
        ######################
        
        
        #################
        ##case mortality
        ## rate calculated from Table 3 and 4 of
        ##Estimates of mortality attributable to influenza and RSV in the United States during 1997-2009 by influenza type or subtype, age, cause of death, and risk status
        ## Goncalo Matias,Robert Taylor,Francois Haguinet, Cynthia Schuck-Paim, Roger Lustig, Vivek Shinde
        ##see mortality rate data sheet for calculations
        low_risk_case_mort_H1_0 = numpy.random.uniform(0.00000052, 0.00000056)
        low_risk_case_mort_H1_5 = numpy.random.uniform(0.00000020, 0.00000024)
        low_risk_case_mort_H1_18 = numpy.random.uniform(0.00000034, 0.00000038)
        low_risk_case_mort_H1_50 = numpy.random.uniform(0.00000000, 0.00000002)
        low_risk_case_mort_H1_65 =numpy.random.uniform(0.00000000, 0.00000002)
        low_risk_case_mort_H1_75 = numpy.random.uniform(0.00000000, 0.00000002)
        
        low_risk_case_mort_H3_0 = numpy.random.uniform(0.00000153, 0.00000157)
        low_risk_case_mort_H3_5 = numpy.random.uniform(0.00000041, 0.00000045)
        low_risk_case_mort_H3_18 = numpy.random.uniform(0.00000240, 0.00000244)
        low_risk_case_mort_H3_50 = numpy.random.uniform(0.00001220, 0.00001224)
        low_risk_case_mort_H3_65 =numpy.random.uniform(0.00004622, 0.00004626)
        low_risk_case_mort_H3_75 = numpy.random.uniform(0.00035179, 0.00035183)
        
        low_risk_case_mort_B_0 = numpy.random.uniform(0.00000089, 0.00000093)
        low_risk_case_mort_B_5 = numpy.random.uniform(0.00000033, 0.00000037)
        low_risk_case_mort_B_18 = numpy.random.uniform(0.00000121, 0.00000123)
        low_risk_case_mort_B_50 = numpy.random.uniform(0.00000476, 0.00000480)
        low_risk_case_mort_B_65 =numpy.random.uniform(0.00001274, 0.00001278)
        low_risk_case_mort_B_75 = numpy.random.uniform(0.00012317, 0.00012321)
        
        
        high_risk_mort_H1_0 = numpy.random.uniform(0.00000018, 0.00000022)
        high_risk_mort_H1_5 = numpy.random.uniform(0.00000022, 0.00000026)
        high_risk_mort_H1_18 = numpy.random.uniform(0.00000050, 0.00000054)
        high_risk_mort_H1_50 = numpy.random.uniform(0.00000005, 0.00000009)
        high_risk_mort_H1_65 =numpy.random.uniform(0.00000010, 0.00000014)
        high_risk_mort_H1_75 = numpy.random.uniform(0.00000004, 0.00000008)
        
        high_risk_mort_H3_0 = numpy.random.uniform(0.00000057, 0.00000061)
        high_risk_mort_H3_5 = numpy.random.uniform(0.00000044, 0.00000048)
        high_risk_mort_H3_18 = numpy.random.uniform(0.00000349, 0.00000353)
        high_risk_mort_H3_50 = numpy.random.uniform(0.00004544, 0.00004548)
        high_risk_mort_H3_65 =numpy.random.uniform(0.00027000, 0.00027004)
        high_risk_mort_H3_75 = numpy.random.uniform(0.00096745, 0.00096749)
        
        high_risk_mort_B_0 = numpy.random.uniform(0.00000033, 0.00000035)
        high_risk_mort_B_5 = numpy.random.uniform(0.00000036, 0.00000040)
        high_risk_mort_B_18 = numpy.random.uniform(0.00000176, 0.00000180)
        high_risk_mort_B_50 = numpy.random.uniform(0.00001776, 0.00001780)
        high_risk_mort_B_65 =numpy.random.uniform(0.00007452, 0.00007456)
        high_risk_mort_B_75 = numpy.random.uniform(0.00033876, 0.00033880)
        ###################
        
        ##################
        #case hospitalization
        ##ref Table 4 of Estimates of hospitalization attributable to influenza and RSV in the US during 1997-2009, by age and risk status
        # Goncalo Matias, Robert Taylor, Francois Haguinet, Cynthia Schuck-Paim, Roger Lustig and Vivek Shinde
        low_risk_case_hosp_H1_0 = numpy.random.triangular(0.00000, 0.00008, 0.00023)
        low_risk_case_hosp_H1_5 = numpy.random.triangular(0.00000, 0.00003, 0.00009)
        low_risk_case_hosp_H1_18 = numpy.random.triangular(0.00000, 0.00001, 0.00003)
        low_risk_case_hosp_H1_50 = numpy.random.triangular(0.00000, 0.00001, 0.00001)
        low_risk_case_hosp_H1_65 = 0
        low_risk_case_hosp_H1_75 = 0
        
        low_risk_case_hosp_H3_0 = numpy.random.triangular(0.00001, 0.00061, 0.00124)
        low_risk_case_hosp_H3_5 = numpy.random.triangular(0.00000, 0.00008, 0.00017)
        low_risk_case_hosp_H3_18 = numpy.random.triangular(0.00000, 0.00010, 0.00020)
        low_risk_case_hosp_H3_50 = numpy.random.triangular(0.00000, 0.00017, 0.00035)
        low_risk_case_hosp_H3_65 = numpy.random.triangular(0.00001, 0.00039, 0.00080)
        low_risk_case_hosp_H3_75 = numpy.random.triangular(0.00003, 0.00117, 0.00235)
        
        low_risk_case_hosp_B_0 = numpy.random.triangular(0.00002, 0.00041, 0.00086)
        low_risk_case_hosp_B_5 = numpy.random.triangular(0.00000, 0.00008, 0.00015)
        low_risk_case_hosp_B_18 = numpy.random.triangular(0.00000, 0.00006, 0.00011)
        low_risk_case_hosp_B_50 = numpy.random.triangular(0.00000, 0.00007, 0.00015)
        low_risk_case_hosp_B_65 = numpy.random.triangular(0.00000, 0.00009, 0.00018)
        low_risk_case_hosp_B_75 = numpy.random.triangular(0.00001, 0.00040, 0.00080)
        
        high_risk_hosp_H1_0 = 0
        high_risk_hosp_H1_5 = numpy.random.triangular(0.00000, 0.00003, 0.00009)
        high_risk_hosp_H1_18 = numpy.random.triangular(0.00000, 0.00009, 0.00030)
        high_risk_hosp_H1_50 = numpy.random.triangular(0.00000, 0.00010, 0.00033)
        high_risk_hosp_H1_65 = numpy.random.triangular(0.00000, 0.00020, 0.00081)
        high_risk_hosp_H1_75 = numpy.random.triangular(0.00000, 0.00021, 0.00086)
        
        high_risk_hosp_H3_0 = numpy.random.triangular(0.00000, 0.00011, 0.00026)
        high_risk_hosp_H3_5 = numpy.random.triangular(0.00000, 0.00004, 0.00008)
        high_risk_hosp_H3_18 = numpy.random.triangular(0.00001, 0.00052, 0.00110)
        high_risk_hosp_H3_50 = numpy.random.triangular(0.00003, 0.00149, 0.00313)
        high_risk_hosp_H3_65 = numpy.random.triangular(0.00006, 0.00286, 0.00591)
        high_risk_hosp_H3_75 = numpy.random.triangular(0.00012, 0.00587, 0.01198)
        
        high_risk_hosp_B_0 = numpy.random.triangular(0.00000, 0.00010, 0.00024)
        high_risk_hosp_B_5 = numpy.random.triangular(0.00000, 0.00003, 0.00007)
        high_risk_hosp_B_18 = numpy.random.triangular(0.00001, 0.00030, 0.00064)
        high_risk_hosp_B_50 = numpy.random.triangular(0.00002, 0.00060, 0.00134)
        high_risk_hosp_B_65 = numpy.random.triangular(0.00001, 0.00044, 0.00105)
        high_risk_hosp_B_75 = numpy.random.triangular(0.00004, 0.00161, 0.00380)
        ##################
        ## from http://rsif.royalsocietypublishing.org/content/early/2011/06/22/rsif.2011.0309#F7
        cross_immunity = numpy.random.uniform(0.2, 0.6)
        R0 = 1.18
        
        elements = [num, rec_rate,  prop_high_risk_0, prop_high_risk_2,prop_high_risk_5,prop_high_risk_19, prop_high_risk_25, prop_high_risk_50, prop_high_risk_65, susc_H1_0, susc_H1_4, susc_H1_18,  susc_H1_65, susc_H3_0, susc_H3_4, susc_H3_18,  susc_H3_65,susc_B_0, susc_B_4, susc_B_18,  susc_B_65,
                    vac_eff_inf_H1_0, vac_eff_inf_H1_6mo, vac_eff_inf_H1_5, vac_eff_inf_H1_18, vac_eff_inf_H1_50,
                    vac_eff_inf_H3_0, vac_eff_inf_H3_6mo, vac_eff_inf_H3_5, vac_eff_inf_H3_18, vac_eff_inf_H3_50,
                    vac_eff_inf_B_0, vac_eff_inf_B_6mo, vac_eff_inf_B_5, vac_eff_inf_B_18, vac_eff_inf_B_50,
                
                    vac_eff_hosp_H1_0, vac_eff_hosp_H1_6mo, vac_eff_hosp_H1_16, vac_eff_hosp_H1_65,
                    vac_eff_hosp_H3_0, vac_eff_hosp_H3_6mo, vac_eff_hosp_H3_16, vac_eff_hosp_H3_65,
                    vac_eff_hosp_B_0, vac_eff_hosp_B_6mo, vac_eff_hosp_B_16, vac_eff_hosp_B_65,

                    high_risk_vac_eff_death_H1_0, high_risk_vac_eff_death_H1_6mo, high_risk_vac_eff_death_H1_18, high_risk_vac_eff_death_H1_65,
                    high_risk_vac_eff_death_H3_0, high_risk_vac_eff_death_H3_6mo, high_risk_vac_eff_death_H3_18, high_risk_vac_eff_death_H3_65,
                    high_risk_vac_eff_death_B_0, high_risk_vac_eff_death_B_6mo, high_risk_vac_eff_death_B_18, high_risk_vac_eff_death_B_65,
                    
                    vac_eff_death_H1_0, vac_eff_death_H1_6mo, vac_eff_death_H1_18, vac_eff_death_H1_65,
                    vac_eff_death_H3_0, vac_eff_death_H3_6mo, vac_eff_death_H3_18, vac_eff_death_H3_65,
                    vac_eff_death_B_0, vac_eff_death_B_6mo, vac_eff_death_B_18, vac_eff_death_B_65,

                    
                    low_risk_case_mort_H1_0, low_risk_case_mort_H1_5, low_risk_case_mort_H1_18, low_risk_case_mort_H1_50, low_risk_case_mort_H1_65, low_risk_case_mort_H1_75,
                    low_risk_case_mort_H3_0, low_risk_case_mort_H3_5, low_risk_case_mort_H3_18, low_risk_case_mort_H3_50, low_risk_case_mort_H3_65, low_risk_case_mort_H3_75,
                    low_risk_case_mort_B_0, low_risk_case_mort_B_5, low_risk_case_mort_B_18, low_risk_case_mort_B_50, low_risk_case_mort_B_65, low_risk_case_mort_B_75,
                    high_risk_mort_H1_0, high_risk_mort_H1_5, high_risk_mort_H1_18, high_risk_mort_H1_50, high_risk_mort_H1_65, high_risk_mort_H1_75,
                    high_risk_mort_H3_0, high_risk_mort_H3_5, high_risk_mort_H3_18, high_risk_mort_H3_50, high_risk_mort_H3_65, high_risk_mort_H3_75,
                    high_risk_mort_B_0, high_risk_mort_B_5, high_risk_mort_B_18, high_risk_mort_B_50, high_risk_mort_B_65, high_risk_mort_B_75,
                    
                    low_risk_case_hosp_H1_0, low_risk_case_hosp_H1_5, low_risk_case_hosp_H1_18, low_risk_case_hosp_H1_50, low_risk_case_hosp_H1_65, low_risk_case_hosp_H1_75,
                    low_risk_case_hosp_H3_0, low_risk_case_hosp_H3_5, low_risk_case_hosp_H3_18, low_risk_case_hosp_H3_50, low_risk_case_hosp_H3_65, low_risk_case_hosp_H3_75,
                    low_risk_case_hosp_B_0, low_risk_case_hosp_B_5, low_risk_case_hosp_B_18, low_risk_case_hosp_B_50, low_risk_case_hosp_B_65, low_risk_case_hosp_B_75,
                    high_risk_hosp_H1_0, high_risk_hosp_H1_5, high_risk_hosp_H1_18, high_risk_hosp_H1_50, high_risk_hosp_H1_65, high_risk_hosp_H1_75,
                    high_risk_hosp_H3_0, high_risk_hosp_H3_5, high_risk_hosp_H3_18, high_risk_hosp_H3_50, high_risk_hosp_H3_65, high_risk_hosp_H3_75,
                    high_risk_hosp_B_0, high_risk_hosp_B_5, high_risk_hosp_B_18, high_risk_hosp_B_50, high_risk_hosp_B_65, high_risk_hosp_B_75, cross_immunity,               
                     R0]
        writer.writerow(elements)
        
        


