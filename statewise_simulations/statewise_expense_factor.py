import sys
sys.path.insert(0, r'../Influenza/Parameters')
from PiecewiseAgeParameter import PiecewiseAgeNumber
import pandas as pd

################################

def return_state_expense(state):
    df = pd.read_csv("./statewise_data/statewise_medical_cost_PS_edited.csv")


    extracted_row = (df.loc[df['state'] == state]).values.tolist()
    multiplication_factor = extracted_row[0][2]
    
    
    return multiplication_factor
    
    
####################################################
if __name__ == "__main__":

    print return_state_expense("New York")
