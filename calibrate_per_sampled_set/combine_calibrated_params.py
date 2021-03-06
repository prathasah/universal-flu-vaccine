import pandas as pd

################################
df1 = pd.read_csv("./calibration_results_March20_2019/results_calibrated_parameters_COMBINED_March20_2019.csv")
print (list(df1))
df1.rename(columns={'incidence(millions)':'data_incidence'}, inplace=True)
df1['data_incidence'] = df1['data_incidence'] *1e6
df2 = pd.read_csv("sampled_parameter_5000_set_with_data_20March2019.csv")
df2.drop(['vac_eff_hospitalization', 'vac_eff_mortality'], axis=1, inplace=True)
###
#combine datasets

merged_df = df2.merge(df1, how = 'inner', on = ['iter'])
merged_df.drop(['Unnamed: 0'], axis=1, inplace=True)
################################
#to drop rows
print merged_df.shape
merged_df.dropna(inplace=True)
#to drop rows
print merged_df.shape
merged_df = merged_df[:5000]
print merged_df.shape
#################################
merged_df.to_csv("sampled_parameters_with_calibration_20March_2019.csv")



	

