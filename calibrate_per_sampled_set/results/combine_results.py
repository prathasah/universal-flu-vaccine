import pandas as pd

################################
df = {}
df1={}
for num in xrange(50):
	filename = "results_calibrated_parameters_num_"+str(num)+".csv"
	df[num] = pd.read_csv(filename)
	filename2 = "./results_efficacy/efficacy_parameters_num_"+str(num)+".csv"
	df1[num] = pd.read_csv(filename2)
	#remove square brackets
	df[num].drop(['vac_eff_hospitalization', 'vac_eff_mortality'], axis=1, inplace=True)
	df[num]['beta_H1'] = pd.DataFrame([str(line).strip('[').strip(']') for line in df[num]['beta_H1']])
	df[num]['beta_H3'] = pd.DataFrame([str(line).strip('[').strip(']') for line in df[num]['beta_H3']])
	df[num]['beta_B'] = pd.DataFrame([str(line).strip('[').strip(']') for line in df[num]['beta_B']])
	df[num]['susceptibility_H1_0'] = pd.DataFrame([str(line).strip('[').strip(']') for line in df[num]['susceptibility_H1_0']])
	df[num]['susceptibility_H1_5'] = pd.DataFrame([str(line).strip('[').strip(']') for line in df[num]['susceptibility_H1_5']])
	df[num]['susceptibility_H1_25'] = pd.DataFrame([str(line).strip('[').strip(']') for line in df[num]['susceptibility_H1_25']])
	df[num]['susceptibility_H1_65'] = pd.DataFrame([str(line).strip('[').strip(']') for line in df[num]['susceptibility_H1_65']])
	df[num]['susceptibility_H3_0'] = pd.DataFrame([str(line).strip('[').strip(']') for line in df[num]['susceptibility_H3_0']])
	df[num]['susceptibility_H3_5'] = pd.DataFrame([str(line).strip('[').strip(']') for line in df[num]['susceptibility_H3_5']])
	df[num]['susceptibility_H3_25'] = pd.DataFrame([str(line).strip('[').strip(']') for line in df[num]['susceptibility_H3_25']])
	df[num]['susceptibility_H3_65'] = pd.DataFrame([str(line).strip('[').strip(']') for line in df[num]['susceptibility_H3_65']])
	df[num]['susceptibility_B_0'] = pd.DataFrame([str(line).strip('[').strip(']') for line in df[num]['susceptibility_B_0']])
	df[num]['susceptibility_B_5'] = pd.DataFrame([str(line).strip('[').strip(']') for line in df[num]['susceptibility_B_5']])
	df[num]['susceptibility_B_25'] = pd.DataFrame([str(line).strip('[').strip(']') for line in df[num]['susceptibility_B_25']])
	df[num]['susceptibility_B_65'] = pd.DataFrame([str(line).strip('[').strip(']') for line in df[num]['susceptibility_B_65']])
	df1[num]['vac_eff_hospitalization'] = pd.DataFrame([str(line).strip('[').strip(']') for line in df1[num]['vac_eff_hospitalization']])
	df1[num]['vac_eff_mortality'] = pd.DataFrame([str(line).strip('[').strip(']') for line in df1[num]['vac_eff_mortality']])
	df1[num]['prob_hosp_scaling'] = pd.DataFrame([str(line).strip('[').strip(']') for line in df1[num]['prob_hosp_scaling']])
	df1[num]['prob_death_scaling'] = pd.DataFrame([str(line).strip('[').strip(']') for line in df1[num]['prob_death_scaling']])
	
	

#combine data-frames
df2 = df[0].copy()
for num in xrange(1,50):
	df2 = df2.append(df[num])


df11 = df1[0].copy()
for num in xrange(1,50):
	df11 = df11.append(df1[num])
	
	
merged_df = df2.merge(df11, how = 'inner', on = ['iter'])
print len(df2['iter'].tolist()), len(df11['iter'].tolist()), len(merged_df['iter'].tolist())

merged_df.drop(['incidence(millions)_y'], axis=1, inplace=True)
merged_df.to_csv("results_calibrated_parameters_COMBINED.csv")