import csv
import numpy
import scipy.stats as stats
import random
#####################################
def gamma_random_sample(mean, variance):
    """Yields a list of random numbers following a gamma distribution defined by mean and variance"""
    g_alpha = mean*mean/variance
    g_beta = mean/variance
    return random.gammavariate(g_alpha,1/g_beta)
######################################################################################		
if __name__ == "__main__":
    
    
    header = ["iter", "infectious_period_0", "infectious_period_15", "proportionHighRisk_0", "proportionHighRisk_2","proportionHighRisk_5","proportionHighRisk_19", "proportionHighRisk_25", "proportionHighRisk_50","proportionHighRisk_65", "susceptibility_H1_0", "susceptibility_H1_4", "susceptibility_H1_18", "susceptibility_H1_65", "susceptibility_H3_0", "susceptibility_H3_4", "susceptibility_H3_18", "susceptibility_H3_65","susceptibility_B_0", "susceptibility_B_4", "susceptibility_B_18", "susceptibility_B_65",
              "relative_vaccineEfficacyVsInfection_H1_0", "relative_vaccineEfficacyVsInfection_H1_0.5",  "relative_vaccineEfficacyVsInfection_H1_5","relative_vaccineEfficacyVsInfection_H1_18","relative_vaccineEfficacyVsInfection_H1_50",
              "relative_vaccineEfficacyVsInfection_H3_0", "relative_vaccineEfficacyVsInfection_H3_0.5",  "relative_vaccineEfficacyVsInfection_H3_5","relative_vaccineEfficacyVsInfection_H3_18","relative_vaccineEfficacyVsInfection_H3_50",
              "relative_vaccineEfficacyVsInfection_B_0", "relative_vaccineEfficacyVsInfection_B_0.5",  "relative_vaccineEfficacyVsInfection_B_5","relative_vaccineEfficacyVsInfection_B_18","relative_vaccineEfficacyVsInfection_B_50",
               "relative_vaccineEfficacyVsHospitalization_H1_0", "relative_vaccineEfficacyVsHospitalization_H1_0.5", "relative_vaccineEfficacyVsHospitalization_H1_16", "relative_vaccineEfficacyVsHospitalization_H1_65",
              "relative_vaccineEfficacyVsHospitalization_H3_0", "relative_vaccineEfficacyVsHospitalization_H3_0.5", "relative_vaccineEfficacyVsHospitalization_H3_16", "relative_vaccineEfficacyVsHospitalization_H3_65",
              "relative_vaccineEfficacyVsHospitalization_B_0", "relative_vaccineEfficacyVsHospitalization_B_0.5", "relative_vaccineEfficacyVsHospitalization_B_16", "relative_vaccineEfficacyVsHospitalization_B_65",
              "relative_vaccineEfficacyVsDeath_H1_0", "relative_vaccineEfficacyVsDeath_H1_0.5", "relative_vaccineEfficacyVsDeath_H1_18", "relative_vaccineEfficacyVsDeath_H1_65",
              "relative_vaccineEfficacyVsDeath_H3_0", "relative_vaccineEfficacyVsDeath_H3_0.5", "relative_vaccineEfficacyVsDeath_H3_18", "relative_vaccineEfficacyVsDeath_H3_65",
              "relative_vaccineEfficacyVsDeath_B_0", "relative_vaccineEfficacyVsDeath_B_0.5", "relative_vaccineEfficacyVsDeath_B_18", "relative_vaccineEfficacyVsDeath_B_65",
              "relative_highRiskvaccineEfficacyVsDeath_H1_0", "relative_highRiskvaccineEfficacyVsDeath_H1_0.5", "relative_highRiskvaccineEfficacyVsDeath_H1_18", "relative_highRiskvaccineEfficacyVsDeath_H1_65",
              "relative_highRiskvaccineEfficacyVsDeath_H3_0", "relative_highRiskvaccineEfficacyVsDeath_H3_0.5", "relative_highRiskvaccineEfficacyVsDeath_H3_18", "relative_highRiskvaccineEfficacyVsDeath_H3_65",
              "relative_highRiskvaccineEfficacyVsDeath_B_0", "relative_highRiskvaccineEfficacyVsDeath_B_0.5", "relative_highRiskvaccineEfficacyVsDeath_B_18", "relative_highRiskvaccineEfficacyVsDeath_B_65",
                "prob_death_0", "prob_death_5", "prob_death_18", "prob_death_50", "prob_death_65",
                "ratio_death_strain_H1_0", "ratio_death_strain_H1_5", "ratio_death_strain_H1_18", "ratio_death_strain_H1_50", "ratio_death_strain_H1_65", "ratio_death_strain_H1_75",
                "ratio_death_strain_H3_0", "ratio_death_strain_H3_5", "ratio_death_strain_H3_18", "ratio_death_strain_H3_50", "ratio_death_strain_H3_65", "ratio_death_strain_H3_75",
                "ratio_death_highrisk_0", "ratio_death_highrisk_5", "ratio_death_highrisk_18", "ratio_death_highrisk_50", "ratio_death_highrisk_65", "ratio_death_highrisk_75",
              "prob_hosp_0", "prob_hosp_5", "prob_hosp_18", "prob_hosp_50", "prob_hosp_65", 
              "ratio_hosp_highrisk_0", "ratio_hosp_highrisk_5", "ratio_hosp_highrisk_18", "ratio_hosp_highrisk_50", "ratio_hosp_highrisk_65","ratio_hosp_highrisk_75",
              "lowRiskhospitalizationRate_H1_0","lowRiskhospitalizationRate_H1_5","lowRiskhospitalizationRate_H1_18","lowRiskhospitalizationRate_H1_50","lowRiskhospitalizationRate_H1_65", "lowRiskhospitalizationRate_H1_75",
              "lowRiskhospitalizationRate_H3_0","lowRiskhospitalizationRate_H3_5","lowRiskhospitalizationRate_H3_18","lowRiskhospitalizationRate_H3_50","lowRiskhospitalizationRate_H3_65", "lowRiskhospitalizationRate_H3_75",
              "lowRiskhospitalizationRate_B_0","lowRiskhospitalizationRate_B_5","lowRiskhospitalizationRate_B_18","lowRiskhospitalizationRate_B_50","lowRiskhospitalizationRate_B_65", "lowRiskhospitalizationRate_B_75",
              "highRiskhospitalizationRate_H1_0","highRiskhospitalizationRate_H1_5","highRiskhospitalizationRate_H1_18","highRiskhospitalizationRate_H1_50","highRiskhospitalizationRate_H1_65", "highRiskhospitalizationRate_H1_75",
              "highRiskhospitalizationRate_H3_0","highRiskhospitalizationRate_H3_5","highRiskhospitalizationRate_H3_18","highRiskhospitalizationRate_H3_50","highRiskhospitalizationRate_H3_65", "highRiskhospitalizationRate_H3_75",
              "highRiskhospitalizationRate_B_0","highRiskhospitalizationRate_B_5","highRiskhospitalizationRate_B_18","highRiskhospitalizationRate_B_50","highRiskhospitalizationRate_B_65", "highRiskhospitalizationRate_B_75", "beta_H1", "beta_H3", "beta_B", "vac_eff_hospitalization", "vac_eff_mortality"]
    writer = csv.writer(open('sampled_parameter_set.csv','wb'))
    writer.writerow(header)
    
    for num in xrange(1000):
        # infectious period from https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.1912
        infectious_period_0 = numpy.random.triangular(2.3,3.6,5.2)
        infectious_period_15 = numpy.random.triangular(3.2,3.9,4.9)
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
        ## for 6mo - 15years: Table 3 of Vaccine effectiveness against laboratoryconfirmed influenza
        ## hospitalizations among young children during the 2010-11 to 2013-14 influenza seasons in Ontario, Canada
        
        ## for 16-64 and 65+ : Table 2 of
        ## Effectiveness of influenza vaccines in preventing severe influenza illness among adults: A systematic review and
        ##meta-analysis of test-negative design case-control studies
        
        relative_vac_eff_hosp_H1_0 = 0
        relative_vac_eff_hosp_H1_6mo = 1
        relative_vac_eff_hosp_H1_16 = numpy.random.triangular(0.34, 0.55, 0.76)/numpy.random.triangular(0.273, 0.821,0.956)
        relative_vac_eff_hosp_H1_65 = numpy.random.triangular(0.26, 0.54, 0.82)/numpy.random.triangular(0.273, 0.821,0.956)
        
        relative_vac_eff_hosp_H3_0 = 0
        relative_vac_eff_hosp_H3_6mo = 1
        relative_vac_eff_hosp_H3_16 = numpy.random.triangular(0.38, 0.50, 0.62)/numpy.random.triangular(0.273, 0.821,0.956)
        relative_vac_eff_hosp_H3_65 = numpy.random.triangular(0.21, 0.33, 0.45)/numpy.random.triangular(0.273, 0.821,0.956)
        
        
        relative_vac_eff_hosp_B_0 = 0
        relative_vac_eff_hosp_B_6mo = 1
        relative_vac_eff_hosp_B_16 = numpy.random.triangular(0.08, 0.45, 0.81)/numpy.random.triangular(0.273, 0.821,0.956)
        relative_vac_eff_hosp_B_65 = numpy.random.triangular(0.11, 0.31, 0.51)/numpy.random.triangular(0.273, 0.821,0.956)
        ###########################
        
        #################################
        ## vaccine efficacy against death
        ##for 6mo-17 years: Table 3 of nfluenza Vaccine Effectiveness Against Pediatric Deaths: 2010-2014; Flanerry et al.
        # for elderly: Table 4 of A Cohort Study of the Effectiveness of Influenza Vaccine in Older People, Performed Using the United Kingdom General Practice Research Database
        # for healthy adults: text in Prioritization of Influenza Pandemic Vaccination to Minimize Years of Life Lost
        # for high risk adults: Table 3 of Clinical Effectiveness of Influenza Vaccination in Persons Younger Than 65 Years With High-Risk Medical Conditions
        
        relative_vac_eff_death_H1_0 = 0
        relative_vac_eff_death_H1_6mo = 1
        relative_vac_eff_death_H1_18 = numpy.random.uniform(0.7, 0.9)/ numpy.random.triangular(0.31, 0.59,0.77)
        relative_vac_eff_death_H1_65 = numpy.random.triangular(0.11, 0.21, 0.29)/ numpy.random.triangular(0.31, 0.59,0.77)
        
        relative_vac_eff_death_H3_0 = 0
        relative_vac_eff_death_H3_6mo = 1
        relative_vac_eff_death_H3_18 = numpy.random.uniform(0.7, 0.9)/ numpy.random.triangular(0.31, 0.59,0.77)
        relative_vac_eff_death_H3_65 = numpy.random.triangular(0.11, 0.21, 0.29)/ numpy.random.triangular(0.31, 0.59,0.77)
        
        relative_vac_eff_death_B_0 = 0
        relative_vac_eff_death_B_6mo = 1
        relative_vac_eff_death_B_18 = numpy.random.uniform(0.7, 0.9)/numpy.random.triangular(0.43, 0.71,0.87)
        relative_vac_eff_death_B_65 = numpy.random.triangular(0.11, 0.21, 0.29)/numpy.random.triangular(0.43, 0.71,0.87)
        
        relative_high_risk_vac_eff_death_H1_0 = 0
        relative_high_risk_vac_eff_death_H1_6mo = 1
        relative_high_risk_vac_eff_death_H1_18 = numpy.random.triangular(0.39, 0.78,0.92)/numpy.random.triangular(0.35, 0.59,0.74)
        relative_high_risk_vac_eff_death_H1_65 = numpy.random.triangular(0.22, 0.29, 0.34)/numpy.random.triangular(0.35, 0.59,0.74)
        
        relative_high_risk_vac_eff_death_H3_0 = 0
        relative_high_risk_vac_eff_death_H3_6mo = 1
        relative_high_risk_vac_eff_death_H3_18 = numpy.random.triangular(0.39, 0.78,0.92)/numpy.random.triangular(0.35, 0.59,0.74)
        relative_high_risk_vac_eff_death_H3_65 = numpy.random.triangular(0.22, 0.29, 0.34)/numpy.random.triangular(0.35, 0.59,0.74)
        
        relative_high_risk_vac_eff_death_B_0 = 0
        relative_high_risk_vac_eff_death_B_6mo = 1
        relative_high_risk_vac_eff_death_B_18 = numpy.random.triangular(0.39, 0.78,0.92)/numpy.random.triangular(0.13, 0.35,0.63)
        relative_high_risk_vac_eff_death_B_65 = numpy.random.triangular(0.22, 0.29, 0.34)/numpy.random.triangular(0.13, 0.35,0.63)
        ######################
        # probability of death from https://www.sciencedirect.com/science/article/pii/S0264410X07003854?via%3Dihub
        prob_death_0 = min(max(numpy.random.normal(loc=0.00004, scale =0.00001),0),1)
        prob_death_5 =  min(max(numpy.random.normal(loc=0.00001, scale = 0),0),1)
        prob_death_18 = min(max(numpy.random.normal(loc=0.00009, scale = 0.00003),0),1)
        prob_death_50 = min(max(numpy.random.normal(loc=0.00134, scale = 0.00045),0),1)
        prob_death_65 = min(max(numpy.random.normal(loc=0.01170, scale = 0.00390),0),1)
        ###########################################################
         #################
        ##case mortality
        ## rate calculated from Table 5 of
        ##Estimates of mortality attributable to influenza and RSV in the United States during 1997-2009 by influenza type or subtype, age, cause of death, and risk status
        ## Goncalo Matias,Robert Taylor,Francois Haguinet, Cynthia Schuck-Paim, Roger Lustig, Vivek Shinde
        ##see mortality rate data sheet for calculations
        
        ratio_death_strain_H1_0 = 0.444
        ratio_death_strain_H1_5 = 0.5882
        ratio_death_strain_H1_18 = 0.1245
        ratio_death_strain_H1_50 = 0.0277
        ratio_death_strain_H1_65 = 0.0016
        ratio_death_strain_H1_75 = 0
        
        ratio_death_strain_H3_0 = 1.333
        ratio_death_strain_H3_5 = 0.8823
        ratio_death_strain_H3_18 = 1.5897
        ratio_death_strain_H3_50 = 2.7944
        ratio_death_strain_H3_65 = 3.1693
        ratio_death_strain_H3_75 = 1.4822
        

        ratio_death_highrisk_0 = 0.38
        ratio_death_highrisk_5 = 1.09
        ratio_death_highrisk_18 = 1.45
        ratio_death_highrisk_50 = 3.72
        ratio_death_highrisk_65 = 5.84
        ratio_death_highrisk_75 = 2.75
        
        ##################
        # probability of hospitalizatio
        
        prob_hosp_0 = numpy.random.triangular(0.0049, 0.0141,0.0233)
        prob_hosp_5 =  numpy.random.triangular(0.0002, 0.0006,0.001)
        prob_hosp_18 = numpy.random.triangular( 0.0015, 0.0042,0.0069)
        prob_hosp_50 = numpy.random.triangular( 0.00676, 0.0193,0.0318)
        prob_hosp_65 = numpy.random.triangular(0.0147, 0.0421,0.0695)
        ###########################################################
        
        ratio_hosp_highrisk_0 = numpy.random.triangular(5,21,42)/numpy.random.triangular(39,110,147)
        ratio_hosp_highrisk_5 = numpy.random.triangular(6,10,16)/numpy.random.triangular(11,19,29)
        ratio_hosp_highrisk_18 = numpy.random.triangular(41,91,145)/numpy.random.triangular(6,17,25)
        ratio_hosp_highrisk_50 = numpy.random.triangular(76,220,358)/numpy.random.triangular(7,25,37)
        ratio_hosp_highrisk_65 = numpy.random.triangular(82,350,594)/numpy.random.triangular(8,49,80)
        ratio_hosp_highrisk_75 = numpy.random.triangular(205,769,1209)/numpy.random.triangular(29,157,238)

        ############################################################
        #case hospitalization
        ##ref Table 4 of Estimates of hospitalization attributable to influenza and RSV in the US during 1997-2009, by age and risk status
        # Goncalo Matias, Robert Taylor, Francois Haguinet, Cynthia Schuck-Paim, Roger Lustig and Vivek Shinde
        #######################################
        
        
        low_risk_hosp_rate_H1_0 = numpy.random.triangular(0,8,23)
        low_risk_hosp_rate_H1_5 = numpy.random.triangular(0,3,9)
        low_risk_hosp_rate_H1_18 = numpy.random.triangular(0,1,3)
        low_risk_hosp_rate_H1_50 =numpy.random.triangular(0,1,1)
        low_risk_hosp_rate_H1_65 = 0
        low_risk_hosp_rate_H1_75 = 0
        
        low_risk_hosp_rate_H3_0 = numpy.random.triangular(1,61,124)
        low_risk_hosp_rate_H3_5 = numpy.random.triangular(0,8,17)
        low_risk_hosp_rate_H3_18 =numpy.random.triangular(0,10,20)
        low_risk_hosp_rate_H3_50 = numpy.random.triangular(0,17,35)
        low_risk_hosp_rate_H3_65 = numpy.random.triangular(1,39,80)
        low_risk_hosp_rate_H3_75 =numpy.random.triangular(3,117,235)
        
        low_risk_hosp_rate_B_0 = numpy.random.triangular(2,41,86)
        low_risk_hosp_rate_B_5 = numpy.random.triangular(0,8,15)
        low_risk_hosp_rate_B_18 = numpy.random.triangular(0,6,11)
        low_risk_hosp_rate_B_50 = numpy.random.triangular(0,7,15)
        low_risk_hosp_rate_B_65 = numpy.random.triangular(0,9,18)
        low_risk_hosp_rate_B_75 =numpy.random.triangular(1,40,80)
        
        high_risk_hosp_rate_H1_0 = 0
        high_risk_hosp_rate_H1_5 = numpy.random.triangular(0,3,9)
        high_risk_hosp_rate_H1_18 = numpy.random.triangular(0,9,30)
        high_risk_hosp_rate_H1_50 = numpy.random.triangular(0,10,33)
        high_risk_hosp_rate_H1_65 = numpy.random.triangular(0,20,81)
        high_risk_hosp_rate_H1_75 = numpy.random.triangular(0,21,86)
        ## the parameter realistically cannot be zero, so compute value based on highrisk/lowrisk ratio
        low_risk_hosp_rate_H1_65 = high_risk_hosp_rate_H1_65/7.9
        low_risk_hosp_rate_H1_75 = high_risk_hosp_rate_H1_65/4.9
        
        high_risk_hosp_rate_H3_0 = numpy.random.triangular(0,11,26)
        high_risk_hosp_rate_H3_5 = numpy.random.triangular(0,4,8)
        high_risk_hosp_rate_H3_18 = numpy.random.triangular(1,52,110)
        high_risk_hosp_rate_H3_50 =numpy.random.triangular(3,149,313)
        high_risk_hosp_rate_H3_65 = numpy.random.triangular(6,286,591)
        high_risk_hosp_rate_H3_75 = numpy.random.triangular(12,587,1198)
        
        high_risk_hosp_rate_B_0 =numpy.random.triangular(0,10,24)
        high_risk_hosp_rate_B_5 = numpy.random.triangular(0,3,7)
        high_risk_hosp_rate_B_18 =numpy.random.triangular(1,30,64)
        high_risk_hosp_rate_B_50 = numpy.random.triangular(2,60,134)
        high_risk_hosp_rate_B_65 =numpy.random.triangular(1,44,105)
        high_risk_hosp_rate_B_75 = numpy.random.triangular(4,161,380)
        ##################
        ## from http://rsif.royalsocietypublishing.org/content/early/2011/06/22/rsif.2011.0309#F7
        #cross_immunity = numpy.random.uniform(0.2, 0.6)
        beta_H1 = 0.05602748
        beta_H3 = 0.05777868
        beta_B = 0.05965933
        vac_eff_hospitalization = 0.1
        vac_eff_mortality = 0.2
        
        elements = [num, infectious_period_0,  infectious_period_15,  prop_high_risk_0, prop_high_risk_2,prop_high_risk_5,prop_high_risk_19, prop_high_risk_25, prop_high_risk_50, prop_high_risk_65, susc_H1_0, susc_H1_4, susc_H1_18,  susc_H1_65, susc_H3_0, susc_H3_4, susc_H3_18,  susc_H3_65,susc_B_0, susc_B_4, susc_B_18,  susc_B_65,
                    vac_eff_inf_H1_0, vac_eff_inf_H1_6mo, vac_eff_inf_H1_5, vac_eff_inf_H1_18, vac_eff_inf_H1_50,
                    vac_eff_inf_H3_0, vac_eff_inf_H3_6mo, vac_eff_inf_H3_5, vac_eff_inf_H3_18, vac_eff_inf_H3_50,
                    vac_eff_inf_B_0, vac_eff_inf_B_6mo, vac_eff_inf_B_5, vac_eff_inf_B_18, vac_eff_inf_B_50,
                
                    relative_vac_eff_hosp_H1_0, relative_vac_eff_hosp_H1_6mo, relative_vac_eff_hosp_H1_16, relative_vac_eff_hosp_H1_65,
                    relative_vac_eff_hosp_H3_0, relative_vac_eff_hosp_H3_6mo, relative_vac_eff_hosp_H3_16, relative_vac_eff_hosp_H3_65,
                    relative_vac_eff_hosp_B_0, relative_vac_eff_hosp_B_6mo, relative_vac_eff_hosp_B_16, relative_vac_eff_hosp_B_65,

                    
                    relative_vac_eff_death_H1_0, relative_vac_eff_death_H1_6mo, relative_vac_eff_death_H1_18, relative_vac_eff_death_H1_65,
                    relative_vac_eff_death_H3_0, relative_vac_eff_death_H3_6mo, relative_vac_eff_death_H3_18, relative_vac_eff_death_H3_65,
                    relative_vac_eff_death_B_0, relative_vac_eff_death_B_6mo, relative_vac_eff_death_B_18, relative_vac_eff_death_B_65,

                    relative_high_risk_vac_eff_death_H1_0, relative_high_risk_vac_eff_death_H1_6mo, relative_high_risk_vac_eff_death_H1_18, relative_high_risk_vac_eff_death_H1_65,
                    relative_high_risk_vac_eff_death_H3_0, relative_high_risk_vac_eff_death_H3_6mo, relative_high_risk_vac_eff_death_H3_18, relative_high_risk_vac_eff_death_H3_65,
                    relative_high_risk_vac_eff_death_B_0, relative_high_risk_vac_eff_death_B_6mo, relative_high_risk_vac_eff_death_B_18, relative_high_risk_vac_eff_death_B_65,
                    
                    prob_death_0, prob_death_5, prob_death_18, prob_death_50, prob_death_65,
                ratio_death_strain_H1_0, ratio_death_strain_H1_5, ratio_death_strain_H1_18, ratio_death_strain_H1_50, ratio_death_strain_H1_65, ratio_death_strain_H1_75,
                ratio_death_strain_H3_0, ratio_death_strain_H3_5, ratio_death_strain_H3_18, ratio_death_strain_H3_50, ratio_death_strain_H3_65, ratio_death_strain_H3_75,
                ratio_death_highrisk_0, ratio_death_highrisk_5, ratio_death_highrisk_18, ratio_death_highrisk_50, ratio_death_highrisk_65, ratio_death_highrisk_75,
                    
                    prob_hosp_0, prob_hosp_5, prob_hosp_18, prob_hosp_50, prob_hosp_65, 
                    ratio_hosp_highrisk_0, ratio_hosp_highrisk_5, ratio_hosp_highrisk_18, ratio_hosp_highrisk_50, ratio_hosp_highrisk_65,ratio_hosp_highrisk_75,
                    low_risk_hosp_rate_H1_0, low_risk_hosp_rate_H1_5, low_risk_hosp_rate_H1_18, low_risk_hosp_rate_H1_50, low_risk_hosp_rate_H1_65, low_risk_hosp_rate_H1_75,
                    low_risk_hosp_rate_H3_0, low_risk_hosp_rate_H3_5, low_risk_hosp_rate_H3_18, low_risk_hosp_rate_H3_50, low_risk_hosp_rate_H3_65, low_risk_hosp_rate_H3_75,
                    low_risk_hosp_rate_B_0, low_risk_hosp_rate_B_5, low_risk_hosp_rate_B_18, low_risk_hosp_rate_B_50, low_risk_hosp_rate_B_65, low_risk_hosp_rate_B_75,
                    high_risk_hosp_rate_H1_0, high_risk_hosp_rate_H1_5, high_risk_hosp_rate_H1_18, high_risk_hosp_rate_H1_50, high_risk_hosp_rate_H1_65, high_risk_hosp_rate_H1_75,
                    high_risk_hosp_rate_H3_0, high_risk_hosp_rate_H3_5, high_risk_hosp_rate_H3_18, high_risk_hosp_rate_H3_50, high_risk_hosp_rate_H3_65, high_risk_hosp_rate_H3_75,
                    high_risk_hosp_rate_B_0, high_risk_hosp_rate_B_5, high_risk_hosp_rate_B_18, high_risk_hosp_rate_B_50, high_risk_hosp_rate_B_65, high_risk_hosp_rate_B_75,                
                    beta_H1, beta_H3, beta_B, vac_eff_hospitalization, vac_eff_mortality]
        writer.writerow(elements)
        
        


