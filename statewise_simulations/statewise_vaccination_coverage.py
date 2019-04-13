import pandas as pd
################################

def return_state_vax_coverage(state):
    df = pd.read_csv("./statewise_data/Clean_data_statewise_avg_vax_coverage_2010-11_2017-18_March20_2019_PS_edited.csv")
    
    
    ## vax coverage for 0.5-4, 5-12, 13-17y, 18-49y, 50-64y and 65+y
    vax_coverage_list_raw = (df.loc[df['State'] == state]).values.tolist()
    vax_coverage_list_raw = vax_coverage_list_raw[0][2:]
    
    ##converting to vax coverage for [0.5, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75]
    vax_coverage_less_5 =  (vax_coverage_list_raw[0])/100
    vax_coverage_5_10 = (vax_coverage_list_raw[1])/100
    vax_coverage_10_15 = (0.5*vax_coverage_list_raw[1]+ 0.5*vax_coverage_list_raw[2])/100
    vax_coverage_15_20 = (0.5*vax_coverage_list_raw[2]+ 0.5*vax_coverage_list_raw[3])/100
    vax_coverage_20_25 = (vax_coverage_list_raw[3])/100
    vax_coverage_25_30 = (vax_coverage_list_raw[3])/100
    vax_coverage_30_35 = (vax_coverage_list_raw[3])/100
    vax_coverage_35_40 = (vax_coverage_list_raw[3])/100
    vax_coverage_40_45 = (vax_coverage_list_raw[3])/100
    vax_coverage_45_50 = (vax_coverage_list_raw[3])/100
    vax_coverage_50_55 = (vax_coverage_list_raw[4])/100
    vax_coverage_55_60 = (vax_coverage_list_raw[4])/100
    vax_coverage_60_65 = (vax_coverage_list_raw[4])/100
    vax_coverage_65_70 = (vax_coverage_list_raw[5])/100
    vax_coverage_70_75 = (vax_coverage_list_raw[5])/100
    vax_coverage_75 = (vax_coverage_list_raw[5])/100
    
    vax_coverage_list = [vax_coverage_less_5, vax_coverage_5_10, vax_coverage_10_15, vax_coverage_15_20, vax_coverage_20_25, vax_coverage_25_30, vax_coverage_30_35, vax_coverage_35_40, vax_coverage_40_45, vax_coverage_45_50, vax_coverage_50_55, vax_coverage_55_60, vax_coverage_60_65, vax_coverage_65_70, vax_coverage_70_75, vax_coverage_75]
    return vax_coverage_list

####################################################
if __name__ == "__main__":

    print return_state_vax_coverage("New York")