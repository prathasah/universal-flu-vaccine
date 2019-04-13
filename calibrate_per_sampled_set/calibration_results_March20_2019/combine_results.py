import pandas as pd

################################
df = {}
for num in xrange(500):
	filename = "results_calibrated_parameters_num_"+str(num)+".csv"
	df[num] = pd.read_csv(filename)
	

#combine data-frames
df2 = df[0].copy()
for num in xrange(1,500):
	df2 = df2.append(df[num])

	

#merged_df.drop(['incidence(millions)_y'], axis=1, inplace=True)
df2.to_csv("results_calibrated_parameters_COMBINED_March20_2019.csv")