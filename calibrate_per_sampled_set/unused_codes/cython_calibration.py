import time
#import pyximport; pyximport.install()
import calibrate_model_cluster_per_set as cb
import pandas as pd
import csv
from multiprocessing import Pool


###################################################

def cons1(x):
    return 1 - max(x)-0.001

def cons2(x):
    return min(x)-0.001

######################################################################33
 
def main():
    start = time.time()
    writer = {}
    elements = {}
    
    df  = pd.read_csv("sampled_parameter_10000_set_with_data.csv")
    
    #########################
    ## create csv for calibrated parameters
    header = ["iter", "incidence(millions)",  "hospitalizations(thousands)", "mortality(thousands)", "beta_H1", "beta_H3", "beta_B", "susceptibility_H1_0", "susceptibility_H1_5", "susceptibility_H1_25", "susceptibility_H1_65", "susceptibility_H3_0", "susceptibility_H3_5", "susceptibility_H3_25", "susceptibility_H3_65","susceptibility_B_0", "susceptibility_B_5", "susceptibility_B_25", "susceptibility_B_65", "vac_eff_hospitalization", "vac_eff_mortality"]
    writer = csv.writer(open('results_calibrated_parameters.csv','wb'))
    writer.writerow(header)

    BL0_burden = [0, 0]
    
    total_doses = df['data_vacDoses'].tolist()
    efficacy = df['data_vacEfficacy'].tolist()
    efficacy = [num/100. for num in efficacy]
    data_infections = df['data_incidence'].tolist()
    data_infections = [num/1e6 for num in data_infections]
    data_hospitalization = df['data_hospitalizations'].tolist()
    data_hospitalization = [num/1e3 for num in data_hospitalization]
    data_mortality = df['data_mortality'].tolist()
    data_mortality = [num/1e3 for num in data_mortality]
    data_perc_H1 = df['data_perc_H1'].tolist()
    data_perc_H3 = df['data_perc_H3'].tolist()
    data_perc_B = df['data_perc_B'].tolist()
    data_incidence_H1_0 = df['data_H1_0'].tolist()
    data_incidence_H1_5 = df['data_H1_5'].tolist()
    data_incidence_H1_25 = df['data_H1_25'].tolist()
    data_incidence_H1_65 = df['data_H1_65'].tolist()
    data_incidence_H3_0 = df['data_H3_0'].tolist()
    data_incidence_H3_5 = df['data_H3_5'].tolist()
    data_incidence_H3_25 = df['data_H3_25'].tolist()
    data_incidence_H3_65 = df['data_H3_65'].tolist()
    data_incidence_B_0 = df['data_B_0'].tolist()
    data_incidence_B_5 = df['data_B_5'].tolist()
    data_incidence_B_25 = df['data_B_25'].tolist()
    data_incidence_B_65 = df['data_B_65'].tolist()
    
    
    
    data_incidence = [total_doses, efficacy, data_infections, data_perc_H1, data_perc_H3, data_perc_B, data_incidence_H1_0, data_incidence_H1_5, data_incidence_H1_25, data_incidence_H1_65, data_incidence_H3_0, data_incidence_H3_5, data_incidence_H3_25, data_incidence_H3_65, data_incidence_B_0, data_incidence_B_5, data_incidence_B_25, data_incidence_B_65]
    
    
    beta_H1, beta_H3, beta_B, H1_0, H1_5, H1_25, H1_65, H3_0, H3_5, H3_25, H3_65, B_0, B_5, B_25, B_65 = cb.optimize(cons1, cons2, data_incidence, beta_dict_opt, sub_indexes[0])
       

    elem1 = [num, data_infections, data_hospitalization, data_mortality, beta_H1, beta_H3, beta_B, H1_0, H1_5, H1_25, H1_65, H3_0, H3_5, H3_25, H3_65,B_0, B_5, B_25, B_65]
    writer.writerow(elem1)

    duration = time.time() - start
    print(result, duration)
 
 
if __name__ == '__main__':
    main()
