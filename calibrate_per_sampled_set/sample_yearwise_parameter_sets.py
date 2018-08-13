import csv
import numpy
import scipy.stats as stats
import random
import pandas as pd
#####################################
def gamma_random_sample(mean, variance):
    """Yields a list of random numbers following a gamma distribution defined by mean and variance"""
    g_alpha = mean*mean/variance
    g_beta = mean/variance
    return random.gammavariate(g_alpha,1/g_beta)
######################################################################################		
if __name__ == "__main__":
    
    
    header = ["iter", "year",  "data_incidence", "data_hospitalizations", "data_mortality", "data_vacEfficacy", "data_vacDoses", "data_H1_0", "data_H1_5", "data_H1_25", "data_H1_65", "data_H3_0", "data_H3_5", "data_H3_25", "data_H3_65", "data_B_0", "data_B_5", "data_B_25", "data_B_65", "infectious_period_0", "infectious_period_15", "proportionHighRisk_0", "proportionHighRisk_2","proportionHighRisk_5","proportionHighRisk_19", "proportionHighRisk_25", "proportionHighRisk_50","proportionHighRisk_65", 
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
              "highRiskhospitalizationRate_B_0","highRiskhospitalizationRate_B_5","highRiskhospitalizationRate_B_18","highRiskhospitalizationRate_B_50","highRiskhospitalizationRate_B_65", "highRiskhospitalizationRate_B_75", "vac_eff_hospitalization", "vac_eff_mortality"]
    writer = csv.writer(open('sampled_yearwise_parameter_set_with_data.csv','wb'))
    writer.writerow(header)
    
    #https://www.cdc.gov/flu/about/disease/2015-16.htm
    #list is min, mode, max for the random triangular distributions
    incidence_data = {}
    incidence_data[2010] = [17582319, 21096749 , 27698870]
    incidence_data[2011] = [7281179, 9231004, 13835345]
    incidence_data[2012] = [30113616, 35590424, 44250092]
    incidence_data[2013] = [24968054, 28445377, 33040119]
    incidence_data[2014] = [30332937, 34292299, 40051029]
    incidence_data[2015] = [21504826, 24577163,28626313]
    
    
    hospitalizations_data = {}
    hospitalizations_data[2010] = [239013, 281589, 373931]
    hospitalizations_data[2011] = [115865, 139497, 206066]
    hospitalizations_data[2012] = [509813, 592688,733307]
    hospitalizations_data[2013] = [283230, 322123, 376646]
    hospitalizations_data[2014] = [624149, 707155, 838516]
    hospitalizations_data[2015] = [271143, 308232, 362029]
        
        
    mortality_data = {}
    mortality_data[2010] = [12111,13541,15372]
    mortality_data[2011] = [3691,4154,4747]
    mortality_data[2012] = [18006,19962, 22434]
    mortality_data[2013] = [12252,13590,15307]
    mortality_data[2014] = [17718,19490,21740]
    mortality_data[2015] = [10634,11995,13914]
    
    #https://www.cdc.gov/flu/professionals/vaccination/effectiveness-studies.htm
    vacEfficacy_data = {}
    vacEfficacy_data[2010] = [53,60,66]
    vacEfficacy_data[2011] = [36,47,56]
    vacEfficacy_data[2012] = [43,49,55]
    vacEfficacy_data[2013] = [44,52,59]
    vacEfficacy_data[2014] = [10,19,27]
    vacEfficacy_data[2015] = [41,48,55]
    
    vacDoses_data = {}
    vacDoses_data[2010] = 155.1e6
    vacDoses_data[2011] = 132e6
    vacDoses_data[2012] = 134.9e6
    vacDoses_data[2013] = 135.5e6
    vacDoses_data[2014] = 147.8e6
    vacDoses_data[2015] = 146.4e6
    
    yearlist = [2010, 2011, 2012, 2013, 2014, 2015]
    
    H1 = {}
    H3 = {}
    B = {}
    
    #from datasheet
    df = pd.read_csv("data_age_specific_incidence.csv")
    
    for year in yearlist:
        H1[year] = {}
        H3[year] = {}
        B[year] = {}
        for age_groups in [0,5,25,65]:
            H1[year][age_groups] = df.loc[(df["year"] == year) & (df['age_group'] ==age_groups), 'H1'].item()
            H3[year][age_groups] = df.loc[(df["year"] == year) & (df['age_group'] ==age_groups), 'H3'].item()
            B[year][age_groups] = df.loc[(df["year"] == year) & (df['age_group'] ==age_groups), 'B'].item()
    
    
    
    
    #for year in yearlist:
    
    for year in yearlist:
        #incidence, hospitalizations, mortality, vacEfficacy, vacDoses, H1_0, H1_5, H1_25, H1_65, H3_0, H3_5, H3_25, H3_65, B_0, B_5, B_25, B_65,
        for num in xrange(1000):
            
            incidence = numpy.random.triangular(incidence_data[year][0], incidence_data[year][1], incidence_data[year][2])
            hospitalizations = numpy.random.triangular(hospitalizations_data[year][0],hospitalizations_data[year][1],hospitalizations_data[year][2])
            mortality = numpy.random.triangular(mortality_data[year][0], mortality_data[year][1], mortality_data[year][2])
            vacEfficacy = numpy.random.triangular(vacEfficacy_data[year][0], vacEfficacy_data[year][1],vacEfficacy_data[year][2])
            
            vacDoses = vacDoses_data[year]
            
            H1_0 = H1[year][0]
            H1_5 = H1[year][5]
            H1_25 = H1[year][25]
            H1_65 = H1[year][65]
            
            H3_0 = H3[year][0]
            H3_5 = H3[year][5]
            H3_25 = H3[year][25]
            H3_65 = H3[year][65]
            
            B_0 = B[year][0]
            B_5 = B[year][5]
            B_25 = B[year][25]
            B_65 = B[year][65]
            ######################
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
            # probability of hospitalization https://bmcinfectdis.biomedcentral.com/articles/10.1186/s12879-015-1193-4 and
            # https://www.sciencedirect.com/science/article/pii/S0264410X07003854
            
            prob_hosp_0 = numpy.random.triangular(0.00708, 0.01410 ,0.02027)
            prob_hosp_5 =  numpy.random.triangular(0.0003,0.0006,0.0008)
            prob_hosp_18 =  numpy.random.triangular(0.002,0.00420,0.005)
            prob_hosp_50 =  numpy.random.triangular(0.00973,0.01930,0.02050)
            prob_hosp_65 = numpy.random.triangular(0.02118,0.04210,0.05037)
            
            
      
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
            
            vac_eff_hospitalization = 0
            vac_eff_mortality = 0
            ##################

        
            elements = [num, year, incidence, hospitalizations, mortality, vacEfficacy, vacDoses, H1_0, H1_5, H1_25, H1_65, H3_0, H3_5, H3_25, H3_65, B_0, B_5, B_25, B_65, infectious_period_0, infectious_period_15,  prop_high_risk_0, prop_high_risk_2,prop_high_risk_5,prop_high_risk_19, prop_high_risk_25, prop_high_risk_50, prop_high_risk_65, 
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
                        
                        prob_hosp_0, prob_hosp_5,  prob_hosp_18,  prob_hosp_50,  prob_hosp_65, 
                        ratio_hosp_highrisk_0, ratio_hosp_highrisk_5, ratio_hosp_highrisk_18, ratio_hosp_highrisk_50, ratio_hosp_highrisk_65,ratio_hosp_highrisk_75,
                        low_risk_hosp_rate_H1_0, low_risk_hosp_rate_H1_5, low_risk_hosp_rate_H1_18, low_risk_hosp_rate_H1_50, low_risk_hosp_rate_H1_65, low_risk_hosp_rate_H1_75,
                        low_risk_hosp_rate_H3_0, low_risk_hosp_rate_H3_5, low_risk_hosp_rate_H3_18, low_risk_hosp_rate_H3_50, low_risk_hosp_rate_H3_65, low_risk_hosp_rate_H3_75,
                        low_risk_hosp_rate_B_0, low_risk_hosp_rate_B_5, low_risk_hosp_rate_B_18, low_risk_hosp_rate_B_50, low_risk_hosp_rate_B_65, low_risk_hosp_rate_B_75,
                        high_risk_hosp_rate_H1_0, high_risk_hosp_rate_H1_5, high_risk_hosp_rate_H1_18, high_risk_hosp_rate_H1_50, high_risk_hosp_rate_H1_65, high_risk_hosp_rate_H1_75,
                        high_risk_hosp_rate_H3_0, high_risk_hosp_rate_H3_5, high_risk_hosp_rate_H3_18, high_risk_hosp_rate_H3_50, high_risk_hosp_rate_H3_65, high_risk_hosp_rate_H3_75,
                        high_risk_hosp_rate_B_0, high_risk_hosp_rate_B_5, high_risk_hosp_rate_B_18, high_risk_hosp_rate_B_50, high_risk_hosp_rate_B_65, high_risk_hosp_rate_B_75, vac_eff_hospitalization, vac_eff_mortality]
            writer.writerow(elements)
        
        


