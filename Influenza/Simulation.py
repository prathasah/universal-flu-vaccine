import fileIO
import random
import numpy
import sys
sys.path.insert(0, r'./Influenza/Parameters')
import Parameters
        
class run_Simulation:
    def __init__(self, options = None, tMin=0 , tMax = 365, paramValues = {}, index=None):
        self.tMax = tMax
	self.tMin = tMin

        if options != None:
            self.options = options
	   
        else:
            from getOptions import getOptions
	    

        # Must wait until after options, where RNG seed is set
       # Parameters = __import__('Parameters', globals())
	
        self.parameters = Parameters.Parameters(index, **paramValues)

        # Initial condition
        self.Y0 = numpy.zeros(60 * self.parameters.ages.size)

        self.hasSolution = False
	
    def comupte_random_vaccination(self, doses):
	
	    self.vacNumbers = doses
	    typical_vax_doses_condition = False
	    universal_vax_doses_condition = False
	    ## chose a random between between 0 and min of (population size of subgroup and doses remaining). Reqd for when max doses = 0
	    typical_vax_remaining = self.vacNumbers[0]
	    vacsUsedTypical = {}
	    ##randomize age groups that go first to avoid universal values chosen for initial age groups
	    random_age_group_list = [num for num in range(1,17)]
	    random.shuffle(random_age_group_list)
	    for num in random_age_group_list:
		vacsUsedTypical[num] = random.randint(0, min(self.parameters.population[num], typical_vax_remaining))
		##find the last value by substracting the sum of vacsUsedTypical
		typical_vax_remaining = self.vacNumbers[0] - sum(vacsUsedTypical.values())
	    population_unvaccinated = [self.parameters.population[0]] + [(a-b) for (a,b) in zip(self.parameters.population[1:], [vacsUsedTypical[age] for age in range(1,17)])]
	   
	    
	    universal_vax_remaining = self.vacNumbers[1]
	    vacsUsedUniversal = {}
	    random.shuffle(random_age_group_list)
	    for num in random_age_group_list:
		vacsUsedUniversal[num] = random.randint(0, min(population_unvaccinated[num],universal_vax_remaining))
		##find the last value by substracting the sum of vacsUsedTypical
		universal_vax_remaining = self.vacNumbers[1] - sum(vacsUsedUniversal.values())
		
	    
	    return [a/(1.*b) for (a,b) in zip([vacsUsedTypical[age] for age in range(1,17)], self.parameters.population[1:])] +  [a/(1.*b) for (a,b) in zip([vacsUsedUniversal[age] for age in range(1,17)], self.parameters.population[1:])]
	    
    def compute_typical_vaccination(self, doses):
	
    
	    
	    # source :  https://www.cdc.gov/flu/fluvaxview/index.htm	
	    #------------------------------------------------------------
	    #Year	|  6m-4y    |5y - 17y  | 18y-49y |50y - 64y|	65y+|
	    #-----------------------------------------------------------
	    #2012-13	|	69.8|	53.1   | 31.1	 |45.1	   |	66.2
	    #2013-14	|	70.4|	55.3   |32.3	 |45.3     |	65
	    #2014-15	|	70.4|	55.8   |33.5	 |47	   |	66.7
	    #2015-16	|	70  |	55.9   |32.7	 |43.6	   | 	63.4
	    #2016-17	|	70  | 	55.6   |33.6	 |45.4	   |	65.3
	    #----------------------------------------------------------
	    
	    ## vaccination uptake in 2016 + adjustment to ensure total vaccines = 150million
	    ##assuming vaccination coverage of high risk age group = 70% ------- TO CHANGE
	    
	
	    vaccination_coverage = [ 0.70,  0.56,  0.56, 0.5152, 0.34, 0.34,0.34,0.34,0.34,0.34,
				   0.46,  0.46, 0.46, 0.653, 0.653, 0.653] +[0.6]*16

	    return vaccination_coverage + vaccination_coverage
	    

    def computeR0(self):
        return self.parameters.computeR0()

    def getLastValues(self):
        return (self.SUL[-1, :], self.IUL_H1[-1, :], self.IUL_H3[-1, :], self.IUL_B[-1, :], self.RUL_H1[-1, :], self.RUL_H3[-1, :], self.RUL_B[-1, :], self.DUL_H1[-1, :], self.DUL_H3[-1, :], self.DUL_B[-1, :], 
		self.SUH[-1, :], self.IUH_H1[-1, :], self.IUH_H3[-1, :], self.IUH_B[-1, :], self.RUH_H1[-1, :], self.RUH_H3[-1, :], self.RUH_B[-1, :], self.DUH_H1[-1, :], self.DUH_H3[-1, :], self.DUH_B[-1, :],
		self.STL[-1, :], self.ITL_H1[-1, :], self.ITL_H3[-1, :], self.ITL_B[-1, :], self.RTL_H1[-1, :], self.RTL_H3[-1, :], self.RTL_B[-1, :], self.DTL_H1[-1, :], self.DTL_H3[-1, :], self.DTL_B[-1, :],
		self.STH[-1, :], self.ITH_H1[-1, :], self.ITH_H3[-1, :], self.ITH_B[-1, :], self.RTH_H1[-1, :], self.RTH_H3[-1, :], self.RTH_B[-1, :],self.DTH_H1[-1, :], self.DTH_H3[-1, :], self.DTH_B[-1, :],    
		self.SNL[-1, :], self.INL_H1[-1, :], self.INL_H3[-1, :], self.INL_B[-1, :], self.RNL_H1[-1, :], self.RNL_H3[-1, :], self.RNL_B[-1, :], self.DNL_H1[-1, :], self.DNL_H3[-1, :], self.DNL_B[-1, :],
		self.SNH[-1, :], self.INH_H1[-1, :], self.INH_H3[-1, :], self.INH_B[-1, :], self.RNH_H1[-1, :], self.RNH_H3[-1, :], self.RNH_B[-1, :], self.DNH_H1[-1, :], self.DNH_H3[-1, :], self.DNH_B[-1, :])
    
    def updateIC(self):
        if not self.hasSolution:
            # S
	    ## SUL
            self.Y0[ 0: : 60] =  (1 - self.parameters.proportionVaccinatedTL -  self.parameters.proportionVaccinatedNL) * self.parameters.population_lowrisk 

	    ## SUH
            self.Y0[ 10: : 60] =  (1 - self.parameters.proportionVaccinatedTH -  self.parameters.proportionVaccinatedNH) * self.parameters.population_highrisk 
  
	    ## STL
            self.Y0[ 20: : 60] = self.parameters.proportionVaccinatedTL * self.parameters.population_lowrisk 
	    
	    ## STH
            self.Y0[ 30: : 60] =  self.parameters.proportionVaccinatedTH * self.parameters.population_highrisk

	    ## SNL
            self.Y0[ 40: : 60] = self.parameters.proportionVaccinatedNL * self.parameters.population_lowrisk 
	    ## SNH 
            self.Y0[ 50: : 60] = self.parameters.proportionVaccinatedNH  * self.parameters.population_highrisk
	
            # I: Add a single infectious person in each age
	    # IUL_H1, IUL_H3, IUL_B
	    self.Y0[ 1: : 60] = numpy.full(self.parameters.ages.size,1000)
	    self.Y0[ 2: : 60] = numpy.full(self.parameters.ages.size,1000)
	    self.Y0[ 3: : 60] = numpy.full(self.parameters.ages.size,1000)
	    
	    # IUH_H1, IUH_H3, IUH_B
            self.Y0[ 11: : 60] = numpy.full(self.parameters.ages.size, 1000)
	    self.Y0[ 12: : 60] = numpy.full(self.parameters.ages.size, 1000)
	    self.Y0[ 13: : 60] = numpy.full(self.parameters.ages.size, 1000)
	    
	    # ITL_H1, ITL_H3, ITL_B
	    self.Y0[ 21: : 60] = numpy.full(self.parameters.ages.size, 1000)
	    self.Y0[ 22: : 60] = numpy.full(self.parameters.ages.size, 1000)
	    self.Y0[ 23: : 60] = numpy.full(self.parameters.ages.size, 1000)
	    
	    # ITH_H1, ITH_H3, ITH_B
	    self.Y0[ 31: : 60] = numpy.full(self.parameters.ages.size, 1000)
	    self.Y0[ 32: : 60] = numpy.full(self.parameters.ages.size, 1000)
	    self.Y0[ 33: : 60] = numpy.full(self.parameters.ages.size, 1000)
	    
	    # INL_H1, INL_H3, INL_B
	    self.Y0[ 41: : 60] = numpy.full(self.parameters.ages.size, 1000)
	    self.Y0[ 42: : 60] = numpy.full(self.parameters.ages.size, 1000)
	    self.Y0[ 43: : 60] = numpy.full(self.parameters.ages.size, 1000)
	    
	    # INH_H1, INH_H3, INH_B
	    self.Y0[ 51: : 60] = numpy.full(self.parameters.ages.size, 1000)
	    self.Y0[ 52: : 60] = numpy.full(self.parameters.ages.size, 1000)
	    self.Y0[ 53: : 60] = numpy.full(self.parameters.ages.size, 1000)

            # S: Remove those new infectious people from the susceptibles
            self.Y0[ 0: : 60] -= (self.Y0[ 1: : 60] + self.Y0[ 2: : 60] + self.Y0[ 3: : 60])
	    self.Y0[ 10: : 60] -= (self.Y0[ 11: : 60] + self.Y0[ 12: : 60] + self.Y0[ 13: : 60])
	    self.Y0[ 20: : 60] -= (self.Y0[ 21: : 60] + self.Y0[ 22: : 60] + self.Y0[ 23: : 60])
	    self.Y0[ 30: : 60] -= (self.Y0[ 31: : 60] + self.Y0[ 32: : 60] + self.Y0[ 33: : 60])
	    self.Y0[ 40: : 60] -= (self.Y0[ 41: : 60] + self.Y0[ 42: : 60] + self.Y0[ 43: : 60])
	    self.Y0[ 50: : 60] -= (self.Y0[ 51: : 60] + self.Y0[ 52: : 60] + self.Y0[ 53: : 60])
 
            # R
	    self.Y0[ 4:  : 60] = 0.
	    self.Y0[ 5:  : 60] = 0.
	    self.Y0[ 6:  : 60] = 0.
	    self.Y0[ 14:  : 60] = 0.
	    self.Y0[ 15:  : 60] = 0.
	    self.Y0[ 16:  : 60] = 0.
	    self.Y0[ 24:  : 60] = 0.
	    self.Y0[ 25:  : 60] = 0.
	    self.Y0[ 26:  : 60] = 0.
	    self.Y0[ 34:  : 60] = 0.
	    self.Y0[ 35:  : 60] = 0.
	    self.Y0[ 36:  : 60] = 0.
	    self.Y0[ 44:  : 60] = 0.
	    self.Y0[ 45:  : 60] = 0.
	    self.Y0[ 46:  : 60] = 0.
	    self.Y0[ 54:  : 60] = 0.
	    self.Y0[ 55:  : 60] = 0.
	    self.Y0[ 56:  : 60] = 0.
	    
	    #Dead class
	    self.Y0[ 7:  : 60] = 0.
	    self.Y0[ 8:  : 60] = 0.
	    self.Y0[ 9:  : 60] = 0.
	    self.Y0[ 17:  : 60] = 0.
	    self.Y0[ 18:  : 60] = 0.
	    self.Y0[ 19:  : 60] = 0.
	    self.Y0[ 27:  : 60] = 0.
	    self.Y0[ 28:  : 60] = 0.
	    self.Y0[ 29:  : 60] = 0.
	    self.Y0[ 37:  : 60] = 0.
	    self.Y0[ 38:  : 60] = 0.
	    self.Y0[ 39:  : 60] = 0.
	    self.Y0[ 47:  : 60] = 0.
	    self.Y0[ 48:  : 60] = 0.
	    self.Y0[ 49:  : 60] = 0.
	    self.Y0[ 57:  : 60] = 0.
	    self.Y0[ 58:  : 60] = 0.
	    self.Y0[ 59:  : 60] = 0.
	    
        else:
	    SUL, IUL_H1, IUL_H3, IUL_B, RUL_H1, RUL_H3, RUL_B, DUL_H1, DUL_H3, DUL_B, 
	    SUH, IUH_H1, IUH_H3, IUH_B, RUH_H1, RUH_H3, RUH_B, DUH_H1, DUH_H3, DUH_B, 
	    STL, ITL_H1, ITL_H3, ITL_B, RTL_H1, RTL_H3, RTL_B, DTL_H1, DTL_H3, DTL_B, 
	    STH, ITH_H1, ITH_H3, ITH_B, RTH_H1, RTH_H3, RTH_B, DTH_H1, DTH_H3, DTH_B, 
	    SNL, INL_H1, INL_H3, INL_B, RNL_H1, RNL_H3, RNL_B, DNL_H1, DNL_H3, DNL_B, 
	    SNH, INH_H1, INH_H3, INH_B, RNH_H1, RNH_H3, RNH_B, DNH_H1, DNH_H3, DNH_B = self.getLastValues()
	    
            self.Y0[ 0 : : 60] = (1 - self.parameters.proportionVaccinatedL) * SUL
	    self.Y0[ 10 : : 60] = (1 - self.parameters.proportionVaccinatedH) * SUH
	    self.Y0[ 20 : : 60] = STL + (1 - self.parameters.proportionVaccinatedTL) * SUL
	    self.Y0[ 30 : : 60] = STH + (1 - self.parameters.proportionVaccinatedTH) * SUH
	    self.Y0[ 40 : : 60] = SNL + (1 - self.parameters.proportionVaccinatedNL) * SUL
	    self.Y0[ 50 : : 60] = SNH + (1 - self.parameters.proportionVaccinatedNH) * SUH
	    
	    #I
	    self.Y0[ 1 : : 60] = IUL_H1
	    self.Y0[ 2 : : 60] = IUL_H3
	    self.Y0[ 3 : : 60] = IUL_B
	    
	    self.Y0[ 11 : : 60] = IUH_H1
	    self.Y0[ 12 : : 60] = IUH_H3
	    self.Y0[ 13 : : 60] = IUH_B
	    
	    self.Y0[ 21 : : 60] = ITL_H1
	    self.Y0[ 22 : : 60] = ITL_H3
	    self.Y0[ 23 : : 60] = ITL_B
	    
	    self.Y0[ 31 : : 60] = ITH_H1
	    self.Y0[ 32 : : 60] = ITH_H3
	    self.Y0[ 33 : : 60] = ITH_B
	    
	    self.Y0[ 41 : : 60] = INL_H1
	    self.Y0[ 42 : : 60] = INL_H3
	    self.Y0[ 43 : : 60] = INL_B
	    
	    self.Y0[ 51 : : 60] = INH_H1
	    self.Y0[ 52 : : 60] = INH_H3
	    self.Y0[ 53 : : 60] = INH_B
	    
	    #R class
	    self.Y0[ 4 : : 60] = RUL_H1
	    self.Y0[ 5 : : 60] = RUL_H3
	    self.Y0[ 6 : : 60] = RUL_B
	    
	    self.Y0[ 14 : : 60] = RUH_H1
	    self.Y0[ 15 : : 60] = RUH_H3
	    self.Y0[ 16 : : 60] = RUH_B
	    
	    self.Y0[ 24 : : 60] = RTL_H1
	    self.Y0[ 25 : : 60] = RTL_H3
	    self.Y0[ 26 : : 60] = RTL_B
	    
	    self.Y0[ 34 : : 60] = RTH_H1
	    self.Y0[ 35 : : 60] = RTH_H3
	    self.Y0[ 36 : : 60] = RTH_B
	    
	    self.Y0[ 44 : : 60] = RNL_H1
	    self.Y0[ 45 : : 60] = RNL_H3
	    self.Y0[ 46 : : 60] = RNL_B
	    
	    self.Y0[ 54 : : 60] = RNH_H1
	    self.Y0[ 55 : : 60] = RNH_H3
	    self.Y0[ 56 : : 60] = RNH_B
	    
	    #D class
	    self.Y0[ 7 : : 60] = DUL_H1
	    self.Y0[ 8 : : 60] = DUL_H3
	    self.Y0[ 9 : : 60] = DUL_B
	    
	    self.Y0[ 17 : : 60] = DUH_H1
	    self.Y0[ 18 : : 60] = DUH_H3
	    self.Y0[ 19 : : 60] = DUH_B
	    
	    self.Y0[ 27 : : 60] = DTL_H1
	    self.Y0[ 28 : : 60] = DTL_H3
	    self.Y0[ 29 : : 60] = DTL_B
	    
	    self.Y0[ 37 : : 60] = DTH_H1
	    self.Y0[ 38 : : 60] = DTH_H3
	    self.Y0[ 39 : : 60] = DTH_B
	    
	    self.Y0[ 47 : : 60] = DNL_H1
	    self.Y0[ 48 : : 60] = DNL_H3
	    self.Y0[ 49 : : 60] = DNL_B
	    
	    self.Y0[ 57 : : 60] = DNH_H1
	    self.Y0[ 58 : : 60] = DNH_H3
	    self.Y0[ 59 : : 60] = DNH_B
	    

    def RHS(self, Y, t):
        '''
        SIR model with multiple host types.
        
        This function gives the right-hand sides of the ODEs.
        '''
        
        # Convert vector to meaningful component vectors

	SUL    = Y[ 0 : : 60]
        IUL_H1 = Y[ 1 : : 60]
	IUL_H3 = Y[ 2 : : 60]
	IUL_B  = Y[ 3 : : 60]
	RUL_H1 = Y[ 4 : : 60]
	RUL_H3 = Y[ 5 : : 60]
	RUL_B =  Y[ 6 : : 60]
	DUL_H1 = Y[ 7 : : 60]
	DUL_H3 = Y[ 8 : : 60]
	DUL_B =  Y[ 9 : : 60]
	
	SUH    = Y[ 10 : : 60]
        IUH_H1 = Y[ 11 : : 60]
	IUH_H3 = Y[ 12 : : 60]
	IUH_B  = Y[ 13 : : 60]
	RUH_H1 = Y[ 14 : : 60]
	RUH_H3 = Y[ 15 : : 60]
	RUH_B =  Y[ 16 : : 60]
	DUH_H1 = Y[ 17 : : 60]
	DUH_H3 = Y[ 18 : : 60]
	DUH_B =  Y[ 19 : : 60]
	
	STL    = Y[ 20 : : 60]
        ITL_H1 = Y[ 21 : : 60]
	ITL_H3 = Y[ 22 : : 60]
	ITL_B  = Y[ 23 : : 60]
	RTL_H1 = Y[ 24 : : 60]
	RTL_H3 = Y[ 25 : : 60]
	RTL_B =  Y[ 26 : : 60]
	DTL_H1 = Y[ 27 : : 60]
	DTL_H3 = Y[ 28 : : 60]
	DTL_B =  Y[ 29 : : 60]
	
	STH    = Y[ 30 : : 60]
        ITH_H1 = Y[ 31 : : 60]
	ITH_H3 = Y[ 32 : : 60]
	ITH_B  = Y[ 33 : : 60]
	RTH_H1 = Y[ 34 : : 60]
	RTH_H3 = Y[ 35 : : 60]
	RTH_B =  Y[ 36 : : 60]
	DTH_H1 = Y[ 37 : : 60]
	DTH_H3 = Y[ 38 : : 60]
	DTH_B =  Y[ 39 : : 60]
	
	SNL    = Y[ 40 : : 60]
        INL_H1 = Y[ 41 : : 60]
	INL_H3 = Y[ 42 : : 60]
	INL_B  = Y[ 43 : : 60]
	RNL_H1 = Y[ 44 : : 60]
	RNL_H3 = Y[ 45 : : 60]
	RNL_B =  Y[ 46 : : 60]
	DNL_H1 = Y[ 47 : : 60]
	DNL_H3 = Y[ 48 : : 60]
	DNL_B =  Y[ 49 : : 60]
	
	SNH    = Y[ 50 : : 60]
        INH_H1 = Y[ 51 : : 60]
	INH_H3 = Y[ 52 : : 60]
	INH_B  = Y[ 53 : : 60]
	RNH_H1 = Y[ 54 : : 60]
	RNH_H3 = Y[ 55 : : 60]
	RNH_B =  Y[ 56 : : 60]
	DNH_H1 = Y[ 57 : : 60]
	DNH_H3 = Y[ 58 : : 60]
	DNH_B =  Y[ 59 : : 60]
	
	
        N =  sum(SUL+ IUL_H1+ IUL_H3+ IUL_B+ RUL_H1 + RUL_H3 + RUL_B + DUL_H1 + DUL_H3 + DUL_B + 
	    SUH+ IUH_H1+ IUH_H3+ IUH_B+ RUH_H1 + RUH_H3 + RUH_B + DUH_H1 + DUH_H3 + DUH_B +
	    STL+ ITL_H1+ ITL_H3+ ITL_B+ RTL_H1 + RTL_H3 + RTL_B + DTL_H1 + DTL_H3 + DTL_B +
	    STH+ ITH_H1+ ITH_H3+ ITH_B+ RTH_H1 + RTH_H3 + RTH_B + DTH_H1 + DTH_H3 + DTH_B + 
	    SNL+ INL_H1+ INL_H3+ INL_B+ RNL_H1 + RNL_H3 + RNL_B + DNL_H1 + DNL_H3 + DNL_B +  
	    SNH+ INH_H1+ INH_H3+ INH_B+ RNH_H1 + RNH_H3 + RNH_B + DNH_H1 + DNH_H3 + DNH_B) 
      
        # The force of infection
        Lambda_H1 = self.parameters.transmissionScaling_H1 * self.parameters.susceptibility_H1\
		    * numpy.dot(self.parameters.contactMatrix, self.parameters.transmissibility * (IUL_H1 + IUH_H1 + ITL_H1 + ITH_H1+ INL_H1+ INH_H1)) / N
	
	Lambda_H3 = self.parameters.transmissionScaling_H3 * self.parameters.susceptibility_H3 \
                 * numpy.dot(self.parameters.contactMatrix, 
                             self.parameters.transmissibility * (IUL_H3 + IUH_H3 + ITL_H3 + ITH_H3 + INL_H3+ INH_H3)) / N
		
	Lambda_B = self.parameters.transmissionScaling_B * self.parameters.susceptibility_B \
                 * numpy.dot(self.parameters.contactMatrix,
                             self.parameters.transmissibility * (IUL_B + IUH_B + ITL_B + ITH_B+ INL_B+ INH_B)) / N
	
	
        
        # The right-hand sides
	
	#UL
        dSUL    = - (Lambda_H1 + Lambda_H3 + Lambda_B) * SUL 
        dIUL_H1 = (Lambda_H1 * SUL) - (self.parameters.recoveryRate + self.parameters.deathRateUL_H1) * IUL_H1
	dIUL_H3 = (Lambda_H3 * SUL) - (self.parameters.recoveryRate + self.parameters.deathRateUL_H3) * IUL_H3
	dIUL_B  = (Lambda_B * SUL) - (self.parameters.recoveryRate + self.parameters.deathRateUL_B) * IUL_B
	dRUL_H1    = self.parameters.recoveryRate * IUL_H1
	dRUL_H3    = self.parameters.recoveryRate * IUL_H3
	dRUL_B    = self.parameters.recoveryRate * IUL_B
	dDUL_H1  = self.parameters.deathRateUL_H1 * IUL_H1
	dDUL_H3  = self.parameters.deathRateUL_H3 * IUL_H3
	dDUL_B  = self.parameters.deathRateUL_B * IUL_B

	
	#UH
        dSUH    = - (Lambda_H1 + Lambda_H3 + Lambda_B) * SUH 
        dIUH_H1 = (Lambda_H1 * SUH) - (self.parameters.recoveryRate + self.parameters.deathRateUH_H1) * IUH_H1
	dIUH_H3 = (Lambda_H3 * SUH) - (self.parameters.recoveryRate + self.parameters.deathRateUH_H3) * IUH_H3
	dIUH_B  = (Lambda_B * SUH) - (self.parameters.recoveryRate + self.parameters.deathRateUH_B) * IUH_B
	dRUH_H1    = self.parameters.recoveryRate * IUH_H1
	dRUH_H3    = self.parameters.recoveryRate * IUH_H3
	dRUH_B    = self.parameters.recoveryRate * IUH_B
	dDUH_H1  = self.parameters.deathRateUH_H1 * IUH_H1
	dDUH_H3  = self.parameters.deathRateUH_H3 * IUH_H3
	dDUH_B  = self.parameters.deathRateUH_B * IUH_B
	
	#TL
	dSTL = - ((1 - self.parameters.vaccineEfficacyVsInfectionTypical_H1) * Lambda_H1 + (1 - self.parameters.vaccineEfficacyVsInfectionTypical_H3) * Lambda_H3+ (1 - self.parameters.vaccineEfficacyVsInfectionTypical_B) * Lambda_B)  *STL

        dITL_H1 = ((1 - self.parameters.vaccineEfficacyVsInfectionTypical_H1) * Lambda_H1 * STL) - (self.parameters.recoveryRate + self.parameters.deathRateVL_H1) * ITL_H1
	dITL_H3 = ((1 - self.parameters.vaccineEfficacyVsInfectionTypical_H3) * Lambda_H3 * STL) - (self.parameters.recoveryRate + self.parameters.deathRateVL_H3) * ITL_H3
	dITL_B  = ((1 - self.parameters.vaccineEfficacyVsInfectionTypical_B) * Lambda_B * STL) - (self.parameters.recoveryRate + self.parameters.deathRateVL_B) * ITL_B
	dRTL_H1    = self.parameters.recoveryRate * ITL_H1
	dRTL_H3    = self.parameters.recoveryRate * ITL_H3
	dRTL_B    = self.parameters.recoveryRate * ITL_B
	dDTL_H1  = self.parameters.deathRateVL_H1 * ITL_H1
	dDTL_H3  =  self.parameters.deathRateVL_H3 * ITL_H3
	dDTL_B   = self.parameters.deathRateVL_B * ITL_B
	
	#TH
	dSTH =- ((1 - self.parameters.vaccineEfficacyVsInfectionTypical_H1) * Lambda_H1 + (1 - self.parameters.vaccineEfficacyVsInfectionTypical_H3) * Lambda_H3+ (1 - self.parameters.vaccineEfficacyVsInfectionTypical_B) * Lambda_B)  *STH
	
        dITH_H1 = ((1 - self.parameters.vaccineEfficacyVsInfectionTypical_H1) *Lambda_H1 * STH) - (self.parameters.recoveryRate + self.parameters.deathRateVH_H1) * ITH_H1
	dITH_H3 = ((1 - self.parameters.vaccineEfficacyVsInfectionTypical_H3) *Lambda_H3 * STH) - (self.parameters.recoveryRate + self.parameters.deathRateVH_H3) * ITH_H3
	dITH_B = ((1 - self.parameters.vaccineEfficacyVsInfectionTypical_B) * Lambda_B * STH) - (self.parameters.recoveryRate +  self.parameters.deathRateVH_B) * ITH_B
	dRTH_H1    = self.parameters.recoveryRate * ITH_H1
	dRTH_H3    = self.parameters.recoveryRate * ITH_H3
	dRTH_B    = self.parameters.recoveryRate * ITH_B
	dDTH_H1   =  self.parameters.deathRateVH_H1 * ITH_H1
	dDTH_H3   = self.parameters.deathRateVH_H3 * ITH_H3
	dDTH_B    = self.parameters.deathRateVH_B * ITH_B
	
	#NL
	dSNL = - ((1 - self.parameters.vaccineEfficacyVsInfectionUniversal_H1) * Lambda_H1 + (1 - self.parameters.vaccineEfficacyVsInfectionUniversal_H3) * Lambda_H3+ (1 - self.parameters.vaccineEfficacyVsInfectionUniversal_B) * Lambda_B)  *SNL
	
        dINL_H1 = ((1 - self.parameters.vaccineEfficacyVsInfectionUniversal_H1) *Lambda_H1 * SNL) - (self.parameters.recoveryRate + self.parameters.deathRateVL_H1) * INL_H1
	dINL_H3 = ((1 - self.parameters.vaccineEfficacyVsInfectionUniversal_H3) *Lambda_H3 * SNL) - (self.parameters.recoveryRate + self.parameters.deathRateVL_H3) * INL_H3
	dINL_B = ((1 - self.parameters.vaccineEfficacyVsInfectionUniversal_B) * Lambda_B * SNL) - (self.parameters.recoveryRate + self.parameters.deathRateVL_B) * INL_B
	dRNL_H1    = self.parameters.recoveryRate * INL_H1
	dRNL_H3    = self.parameters.recoveryRate * INL_H3
	dRNL_B    = self.parameters.recoveryRate * INL_B
	dDNL_H1  =  self.parameters.deathRateVL_H1 * INL_H1
	dDNL_H3  =  self.parameters.deathRateVL_H3 * INL_H3
	dDNL_B  =   self.parameters.deathRateVL_B * INL_B
	
	#NH
	dSNH =  - ((1 - self.parameters.vaccineEfficacyVsInfectionUniversal_H1) * Lambda_H1 + (1 - self.parameters.vaccineEfficacyVsInfectionUniversal_H3) * Lambda_H3+ (1 - self.parameters.vaccineEfficacyVsInfectionUniversal_B) * Lambda_B)  *SNH
        dINH_H1 = ((1 - self.parameters.vaccineEfficacyVsInfectionUniversal_H1) *Lambda_H1 * SNH) - (self.parameters.recoveryRate + self.parameters.deathRateVH_H1) * INH_H1
	dINH_H3 = ((1 - self.parameters.vaccineEfficacyVsInfectionUniversal_H3) *Lambda_H3 * SNH) - (self.parameters.recoveryRate + self.parameters.deathRateVH_H3) * INH_H3
	dINH_B = ((1 - self.parameters.vaccineEfficacyVsInfectionUniversal_B) *Lambda_B * SNH) - (self.parameters.recoveryRate +  self.parameters.deathRateVH_B) * INH_B
	dRNH_H1    = self.parameters.recoveryRate * INH_H1
	dRNH_H3    = self.parameters.recoveryRate * INH_H3
	dRNH_B    = self.parameters.recoveryRate * INH_B
	dDNH_H1  = self.parameters.deathRateVH_H1 * INH_H1
	dDNH_H3  = self.parameters.deathRateVH_H3 * INH_H3
	dDNH_B  = self.parameters.deathRateVH_B * INH_B
	
	
        # Convert meaningful component vectors into a single vector
        dY = numpy.empty(Y.size, dtype = float)
        dY[ 0 : : 60] = dSUL
        dY[ 1 : : 60] = dIUL_H1
        dY[ 2 : : 60] = dIUL_H3
	dY[ 3 : : 60] = dIUL_B
	dY[ 4 : : 60] = dRUL_H1
	dY[ 5 : : 60] = dRUL_H3
	dY[ 6 : : 60] = dRUL_B
	dY[ 7 : : 60] = dDUL_H1
	dY[ 8 : : 60] = dDUL_H3
	dY[ 9 : : 60] = dDUL_B
	
	dY[ 10 : : 60] = dSUH
        dY[ 11 : : 60] = dIUH_H1
        dY[ 12 : : 60] = dIUH_H3
	dY[ 13 : : 60] = dIUH_B
	dY[ 14 : : 60] = dRUH_H1
	dY[ 15 : : 60] = dRUH_H3
	dY[ 16 : : 60] = dRUH_B
	dY[ 17 : : 60] = dDUH_H1
	dY[ 18 : : 60] = dDUH_H3
	dY[ 19 : : 60] = dDUH_B
	
	dY[ 20 : : 60] = dSTL
        dY[ 21 : : 60] = dITL_H1
        dY[ 22 : : 60] = dITL_H3
	dY[ 23 : : 60] = dITL_B
	dY[ 24 : : 60] = dRTL_H1
	dY[ 25 : : 60] = dRTL_H3
	dY[ 26 : : 60] = dRTL_B
	dY[ 27 : : 60] = dDTL_H1
	dY[ 28 : : 60] = dDTL_H3
	dY[ 29 : : 60] = dDTL_B
	
	dY[ 30 : : 60] = dSTH
        dY[ 31 : : 60] = dITH_H1
        dY[ 32 : : 60] = dITH_H3
	dY[ 33 : : 60] = dITH_B
	dY[ 34 : : 60] = dRTH_H1
	dY[ 35 : : 60] = dRTH_H3
	dY[ 36 : : 60] = dRTH_B
	dY[ 37 : : 60] = dDTH_H1
	dY[ 38 : : 60] = dDTH_H3
	dY[ 39 : : 60] = dDTH_B
	
	dY[ 40 : : 60] = dSNL
        dY[ 41 : : 60] = dINL_H1
        dY[ 42 : : 60] = dINL_H3
	dY[ 43 : : 60] = dINL_B
	dY[ 44 : : 60] = dRNL_H1
	dY[ 45 : : 60] = dRNL_H3
	dY[ 46 : : 60] = dRNL_B
	dY[ 47 : : 60] = dDNL_H1
	dY[ 48 : : 60] = dDNL_H3
	dY[ 49 : : 60] = dDNL_B
	
	dY[ 50 : : 60] = dSNH
        dY[ 51 : : 60] = dINH_H1
        dY[ 52 : : 60] = dINH_H3
	dY[ 53 : : 60] = dINH_B
	dY[ 54 : : 60] = dRNH_H1
	dY[ 55 : : 60] = dRNH_H3
	dY[ 56 : : 60] = dRNH_B
	dY[ 57 : : 60] = dDNH_H1
	dY[ 58 : : 60] = dDNH_H3
	dY[ 59 : : 60] = dDNH_B
	
        return dY
    
    def resetSolution(self):
        self.hasSolution = False

    def solve(self, tStart = 0., tEnd = None, tStep = 1.):
        if tEnd == None:
            tEnd = self.tMax

        if self.hasSolution:
            TOld  = self.T.copy()
	    
            SUL_Old = self.SUL.copy()
            IUL_H1_Old = self.IUL_H1.copy()
	    IUL_H3_Old = self.IUL_H3.copy()
	    IUL_B_Old = self.IUL_B.copy()
	    RUL_H1_Old = self.RUL_H1.copy()
	    RUL_H3_Old = self.RUL_H3.copy()
	    RUL_B_Old = self.RUL_B.copy()
	    DUL_H1_Old = self.DUL_H1.copy()
	    DUL_H3_Old = self.DUL_H3.copy()
	    DUL_B_Old = self.DUL_B.copy()
	    
	    SUH_Old = self.SUH.copy()
            IUH_H1_Old = self.IUH_H1.copy()
	    IUH_H3_Old = self.IUH_H3.copy()
	    IUH_B_Old = self.IUH_B.copy()
	    RUH_H1_Old = self.RUH_H1.copy()
	    RUH_H3_Old = self.RUH_H3.copy()
	    RUH_B_Old = self.RUH_B.copy()
	    DUH_H1_Old = self.DUH_H1.copy()
	    DUH_H3_Old = self.DUH_H3.copy()
	    DUH_B_Old = self.DUH_B.copy()
	    
	    STL_Old = self.STL.copy()
            ITL_H1_Old = self.ITL_H1.copy()
	    ITL_H3_Old = self.ITL_H3.copy()
	    ITL_B_Old = self.ITL_B.copy()
	    RTL_H1_Old = self.RTL_H1.copy()
	    RTL_H3_Old = self.RTL_H3.copy()
	    RTL_B_Old = self.RTL_B.copy()
	    DTL_H1_Old = self.DTL_H1.copy()
	    DTL_H3_Old = self.DTL_H3.copy()
	    DTL_B_Old = self.DTL_B.copy()
	    
	    STH_Old = self.STH.copy()
            ITH_H1_Old = self.ITH_H1.copy()
	    ITH_H3_Old = self.ITH_H3.copy()
	    ITH_B_Old = self.ITH_B.copy()
	    RTH_H1_Old = self.RTH_H1.copy()
	    RTH_H3_Old = self.RTH_H3.copy()
	    RTH_B_Old = self.RTH_B.copy()
	    DTH_H1_Old = self.DTH_H1.copy()
	    DTH_H3_Old = self.DTH_H3.copy()
	    DTH_B_Old = self.DTH_B.copy()
	    
	    SNL_Old = self.SNL.copy()
            INL_H1_Old = self.INL_H1.copy()
	    INL_H3_Old = self.INL_H3.copy()
	    INL_B_Old = self.INL_B.copy()
	    RNL_H1_Old = self.RNL_H1.copy()
	    RNL_H3_Old = self.RNL_H3.copy()
	    RNL_B_Old = self.RNL_B.copy()
	    DNL_H1_Old = self.DNL_H1.copy()
	    DNL_H3_Old = self.DNL_H3.copy()
	    DNL_B_Old = self.DNL_B.copy()
	    
	    SNH_Old = self.SNH.copy()
            INH_H1_Old = self.INH_H1.copy()
	    INH_H3_Old = self.INH_H3.copy()
	    INH_B_Old = self.INH_B.copy()
	    RNH_H1_Old = self.RNH_H1.copy()
	    RNH_H3_Old = self.RNH_H3.copy()
	    RNH_B_Old = self.RNH_B.copy()
	    DNH_H1_Old = self.DNH_H1.copy()
	    DNH_H3_Old = self.DNH_H3.copy()
	    DNH_B_Old = self.DNH_B.copy()
	
        # Time vector for solution
        self.T = numpy.hstack((numpy.arange(tStart, tEnd, tStep), tEnd))
        
        # Integrate the ODE
        from scipy.integrate import odeint
        self.Y = odeint(self.RHS,
                        self.Y0.copy(),
                        self.T,
                        mxstep = 1000)
        Z = self.Y.copy()
	
	self.SUL    = Z[:, 0 : : 60]
        self.IUL_H1 = Z[:, 1 : : 60]
	self.IUL_H3 = Z[:, 2 : : 60]
	self.IUL_B  = Z[:, 3 : : 60]
	self.RUL_H1 = Z[:, 4 : : 60]
	self.RUL_H3 = Z[:, 5 : : 60]
	self.RUL_B  = Z[:, 6 : : 60]
	self.DUL_H1 = Z[:, 7 : : 60]
	self.DUL_H3 = Z[:, 8 : : 60]
	self.DUL_B  = Z[:, 9 : : 60]
	
	self.SUH    = Z[:, 10 : : 60]
        self.IUH_H1 = Z[:, 11 : : 60]
	self.IUH_H3 = Z[:, 12 : : 60]
	self.IUH_B  = Z[:, 13 : : 60]
	self.RUH_H1 = Z[:, 14 : : 60]
	self.RUH_H3 = Z[:, 15 : : 60]
	self.RUH_B  = Z[:, 16 : : 60]
	self.DUH_H1 = Z[:, 17 : : 60]
	self.DUH_H3 = Z[:, 18 : : 60]
	self.DUH_B  = Z[:, 19 : : 60]
	
	self.STL    = Z[:, 20 : : 60]
        self.ITL_H1 = Z[:, 21 : : 60]
	self.ITL_H3 = Z[:, 22 : : 60]
	self.ITL_B  = Z[:, 23 : : 60]
	self.RTL_H1 = Z[:, 24 : : 60]
	self.RTL_H3 = Z[:, 25 : : 60]
	self.RTL_B  = Z[:, 26 : : 60]
	self.DTL_H1 = Z[:, 27 : : 60]
	self.DTL_H3 = Z[:, 28 : : 60]
	self.DTL_B  = Z[:, 29 : : 60]
	
	self.STH    = Z[:, 30 : : 60]
        self.ITH_H1 = Z[:, 31 : : 60]
	self.ITH_H3 = Z[:, 32 : : 60]
	self.ITH_B  = Z[:, 33 : : 60]
	self.RTH_H1 = Z[:, 34 : : 60]
	self.RTH_H3 = Z[:, 35 : : 60]
	self.RTH_B  = Z[:, 36 : : 60]
	self.DTH_H1 = Z[:, 37 : : 60]
	self.DTH_H3 = Z[:, 38 : : 60]
	self.DTH_B  = Z[:, 39 : : 60]
	
	self.SNL    = Z[:, 40 : : 60]
        self.INL_H1 = Z[:, 41 : : 60]
	self.INL_H3 = Z[:, 42 : : 60]
	self.INL_B  = Z[:, 43 : : 60]
	self.RNL_H1 = Z[:, 44 : : 60]
	self.RNL_H3 = Z[:, 45 : : 60]
	self.RNL_B  = Z[:, 46 : : 60]
	self.DNL_H1 = Z[:, 47 : : 60]
	self.DNL_H3 = Z[:, 48 : : 60]
	self.DNL_B  = Z[:, 49 : : 60]
	
	self.SNH    = Z[:, 50 : : 60]
        self.INH_H1 = Z[:, 51 : : 60]
	self.INH_H3 = Z[:, 52 : : 60]
	self.INH_B  = Z[:, 53 : : 60]
	self.RNH_H1 = Z[:, 54 : : 60]
	self.RNH_H3 = Z[:, 55 : : 60]
	self.RNH_B  = Z[:, 56 : : 60]
	self.DNH_H1 = Z[:, 57 : : 60]
	self.DNH_H3 = Z[:, 58 : : 60]
	self.DNH_B  = Z[:, 59 : : 60]
	
	
        if self.hasSolution:
	   
            self.T = numpy.hstack((TOld, self.T))
	    #UL
            self.SUL = numpy.vstack((SUL_Old, self.SUL))
            self.IUL_H1 = numpy.vstack((IUL_H1_Old, self.IUL_H1))
	    self.IUL_H3 = numpy.vstack((IUL_H3_Old, self.IUL_H3))
	    self.IUL_B = numpy.vstack((IUL_B_Old, self.IUL_B))
	    self.RUL_H1 =  numpy.vstack((RUL_H1_Old, self.RUL_H1))
	    self.RUL_H3 =  numpy.vstack((RUL_H3_Old, self.RUL_H3))
	    self.RUL_B =  numpy.vstack((RUL_B_Old, self.RUL_B))
	    self.DUL_H1 =  numpy.vstack((DUL_H1_Old, self.DUL_H1))
	    self.DUL_H3 =  numpy.vstack((DUL_H3_Old, self.DUL_H3))
	    self.DUL_B =  numpy.vstack((DUL_B_Old, self.DUL_B))
	    
	    #UH
	    self.SUH = numpy.vstack((SUH_Old, self.SUH))
            self.IUH_H1 = numpy.vstack((IUH_H1_Old, self.IUH_H1))
	    self.IUH_H3 = numpy.vstack((IUH_H3_Old, self.IUH_H3))
	    self.IUH_B = numpy.vstack((IUH_B_Old, self.IUH_B))
	    self.RUH_H1 =  numpy.vstack((RUH_H1_Old, self.RUH_H1))
	    self.RUH_H3 =  numpy.vstack((RUH_H3_Old, self.RUH_H3))
	    self.RUH_B =  numpy.vstack((RUH_B_Old, self.RUH_B))
	    self.DUH_H1 =  numpy.vstack((DUH_H1_Old, self.DUH_H1))
	    self.DUH_H3 =  numpy.vstack((DUH_H3_Old, self.DUH_H3))
	    self.DUH_B =  numpy.vstack((DUH_B_Old, self.DUH_B))
	    
	    #TL
	    self.STL = numpy.vstack((STL_Old, self.STL))
            self.ITL_H1 = numpy.vstack((ITL_H1_Old, self.ITL_H1))
	    self.ITL_H3 = numpy.vstack((ITL_H3_Old, self.ITL_H3))
	    self.ITL_B = numpy.vstack((ITL_B_Old, self.ITL_B))
	    self.RTL_H1 =  numpy.vstack((RTL_H1_Old, self.RTL_H1))
	    self.RTL_H3 =  numpy.vstack((RTL_H3_Old, self.RTL_H3))
	    self.RTL_B =  numpy.vstack((RTL_B_Old, self.RTL_B))
	    self.DTL_H1 =  numpy.vstack((DTL_H1_Old, self.DTL_H1))
	    self.DTL_H3 =  numpy.vstack((DTL_H3_Old, self.DTL_H3))
	    self.DTL_B =  numpy.vstack((DTL_B_Old, self.DTL_B))
	    
	    #TH
	    self.STH = numpy.vstack((STH_Old, self.STH))
            self.ITH_H1 = numpy.vstack((ITH_H1_Old, self.ITH_H1))
	    self.ITH_H3 = numpy.vstack((ITH_H3_Old, self.ITH_H3))
	    self.ITH_B = numpy.vstack((ITH_B_Old, self.ITH_B))
	    self.RTH_H1 =  numpy.vstack((RTH_H1_Old, self.RTH_H1))
	    self.RTH_H3 =  numpy.vstack((RTH_H3_Old, self.RTH_H3))
	    self.RTH_B =  numpy.vstack((RTH_B_Old, self.RTH_B))
	    self.DTH_H1 =  numpy.vstack((DTH_H1_Old, self.DTH_H1))
	    self.DTH_H3 =  numpy.vstack((DTH_H3_Old, self.DTH_H3))
	    self.DTH_B =  numpy.vstack((DTH_B_Old, self.DTH_B))
	    
	    #NL
	    self.SNL = numpy.vstack((SNL_Old, self.SNL))
            self.INL_H1 = numpy.vstack((INL_H1_Old, self.INL_H1))
	    self.INL_H3 = numpy.vstack((INL_H3_Old, self.INL_H3))
	    self.INL_B = numpy.vstack((INL_B_Old, self.INL_B))
	    self.RNL_H1 =  numpy.vstack((RNL_H1_Old, self.RNL_H1))
	    self.RNL_H3 =  numpy.vstack((RNL_H3_Old, self.RNL_H3))
	    self.RNL_B =  numpy.vstack((RNL_B_Old, self.RNL_B))
	    self.DNL_H1 =  numpy.vstack((DNL_H1_Old, self.DNL_H1))
	    self.DNL_H3 =  numpy.vstack((DNL_H3_Old, self.DNL_H3))
	    self.DNL_B =  numpy.vstack((DNL_B_Old, self.DNL_B))
	    
	    #NH
	    self.SNH = numpy.vstack((SNH_Old, self.SNH))
            self.INH_H1 = numpy.vstack((INH_H1_Old, self.INH_H1))
	    self.INH_H3 = numpy.vstack((INH_H3_Old, self.INH_H3))
	    self.INH_B = numpy.vstack((INH_B_Old, self.INH_B))
	    self.RNH_H1 =  numpy.vstack((RNH_H1_Old, self.RNH_H1))
	    self.RNH_H3 =  numpy.vstack((RNH_H3_Old, self.RNH_H3))
	    self.RNH_B =  numpy.vstack((RNH_B_Old, self.RNH_B))
	    self.DNH_H1 =  numpy.vstack((DNH_H1_Old, self.DNH_H1))
	    self.DNH_H3 =  numpy.vstack((DNH_H3_Old, self.DNH_H3))
	    self.DNH_B =  numpy.vstack((DNH_B_Old, self.DNH_B))
	    
	    
        self.hasSolution = True

    def updateStats(self):
	    
	self.NUL = self.SUL +  self.IUL_H1 + self.IUL_H3 + self.IUL_B + self.RUL_H1 +  self.RUL_H3 +  self.RUL_B + self.DUL_H1 +  self.DUL_H3 +  self.DUL_B  
	self.NUH = self.SUH +  self.IUH_H1 + self.IUH_H3 + self.IUH_B + self.RUH_H1 +  self.RUH_H3 +  self.RUH_B + self.DUH_H1 +  self.DUH_H3 +  self.DUH_B 
	self.NTL = self.STL +  self.ITL_H1 + self.ITL_H3 + self.ITL_B + self.RTL_H1 +  self.RTL_H3 +  self.RTL_B + self.DTL_H1 +  self.DTL_H3 +  self.DTL_B 
	self.NTH = self.STH +  self.ITH_H1 + self.ITH_H3 + self.ITH_B + self.RTH_H1 +  self.RTH_H3 +  self.RTH_B + self.DTH_H1 +  self.DTH_H3 +  self.DTH_B 
	self.NNL = self.SNL +  self.INL_H1 + self.INL_H3 + self.INL_B + self.RNL_H1 +  self.RNL_H3 +  self.RNL_B + self.DNL_H1 +  self.DNL_H3 +  self.DNL_B 
	self.NNH = self.SNH +  self.INH_H1 + self.INH_H3 + self.INH_B + self.RNH_H1 +  self.RNH_H3 +  self.RNH_B + self.DNH_H1 +  self.DNH_H3 +  self.DNH_B 
	
	self.NU = self.NUL + self.NUH
	self.NT = self.NTL + self.NTH
	self.NN = self.NNL + self.NNH
	
        self.N  = self.NU + self.NT + self.NN

        self.infectionsUL = self.NUL[0,:] - self.SUL[-1,:]
	self.infectionsUH = self.NUH[0,:] - self.SUH[-1,:]
	self.infectionsTL = self.NTL[0,:] - self.STL[-1,:]
	self.infectionsTH = self.NTH[0,:] - self.STH[-1,:]
	self.infectionsNL = self.NNL[0,:] - self.SNL[-1,:]
	self.infectionsNH = self.NNH[0,:] - self.SNH[-1,:]
	
	self.infectionsUL_H1 =  self.RUL_H1[-1,: ] + self.DUL_H1[-1,: ] + self.IUL_H1[-1,: ]
	self.infectionsUL_H3 =  self.RUL_H3[-1,: ] + self.DUL_H3[-1,: ] + self.IUL_H3[-1,: ]
	self.infectionsUL_B  =  self.RUL_B[-1,: ]  + self.DUL_B[-1,: ] + self.IUL_B[-1,: ]
	
	self.infectionsVL_H1 =  self.RTL_H1[-1,: ] + self.RNL_H1[-1,: ] + self.DTL_H1[-1,: ] + self.DNL_H1[-1,: ] + self.ITL_H1[-1,: ] + self.INL_H1[-1,: ]
	self.infectionsVL_H3 =  self.RTL_H3[-1,: ] + self.RNL_H3[-1,: ] + self.DTL_H3[-1,: ] + self.DNL_H3[-1,: ] + self.ITL_H3[-1,: ] + self.INL_H3[-1,: ]
	self.infectionsVL_B  =  self.RTL_B[-1,: ]  + self.RNL_B[-1,: ]  + self.DTL_B[-1,: ] + self.DNL_B[-1,: ] + self.ITL_B[-1,: ] + self.INL_B[-1,: ]
	
	self.infectionsUH_H1 =  self.RUH_H1[-1,: ] + self.DUH_H1[-1,: ] + self.IUH_H1[-1,: ]
	self.infectionsUH_H3 =  self.RUH_H3[-1,: ] + self.DUH_H3[-1,: ] + self.IUH_H3[-1,: ]
	self.infectionsUH_B  =  self.RUH_B[-1,: ]  + self.DUH_B[-1,: ] + self.IUH_B[-1,: ]
	
	self.infectionsVH_H1 =  self.RTH_H1[-1,: ] + self.RNH_H1[-1,: ] + self.DTH_H1[-1,: ] + self.DNH_H1[-1,: ] + self.ITH_H1[-1,: ] + self.INH_H1[-1,: ]
	self.infectionsVH_H3 =  self.RTH_H3[-1,: ] + self.RNH_H3[-1,: ] + self.DTH_H3[-1,: ] + self.DNH_H3[-1,: ] + self.ITH_H3[-1,: ] + self.INH_H3[-1,: ]
	self.infectionsVH_B  =  self.RTH_B[-1,: ]  + self.RNH_B[-1,: ]  + self.DTH_B[-1,: ] + self.DNH_B[-1,: ]   + self.ITH_B[-1,: ] + self.INH_B[-1,: ]

        # Find duplicate times: these are where vaccination occurs
        #for i in numpy.arange(len(self.T)).compress(numpy.diff(self.T) == 0):
            # Update for vaccinated people
        #    self.infectionsUL += self.SUL[i + 1, :] - self.SUL[i, :]
	#    self.infectionsUH += self.SUH[i + 1, :] - self.SUH[i, :]
        #    self.infectionsTL += self.STL[i + 1, :] - self.STL[i, :]
	#    self.infectionsTH += self.STH[i + 1, :] - self.STH[i, :]
	#    self.infectionsNL += self.SNL[i + 1, :] - self.SNL[i, :]
	#    self.infectionsNH += self.SNH[i + 1, :] - self.SNH[i, :]

	self.infections_H1 = self.infectionsUL_H1 + self.infectionsUH_H1 + self.infectionsVL_H1 + self.infectionsVH_H1 
	self.infections_H3 = self.infectionsUL_H3 + self.infectionsUH_H3 + self.infectionsVL_H3 + self.infectionsVH_H3 
	self.infections_B = self.infectionsUL_B + self.infectionsUH_B + self.infectionsVL_B + self.infectionsVH_B 
	
	self.infectionsL  = self.infectionsUL + self.infectionsTL + self.infectionsNL
        self.infectionsH  = self.infectionsUH + self.infectionsTH + self.infectionsNH
	self.infectionsU =  self.infectionsUL + self.infectionsUH
	self.infectionsV  = self.infectionsTL + self.infectionsTH + self.infectionsNL + self.infectionsNH        
	self.infections  = self.infectionsU + self.infectionsV
        self.totalInfections = self.infections.sum()

	self.hospitalizationsUL_H1 = self.infectionsUL_H1 * self.parameters.lowRiskcaseHospitalization_H1
	self.hospitalizationsUL_H3 = self.infectionsUL_H3 * self.parameters.lowRiskcaseHospitalization_H3
	self.hospitalizationsUL_B = self.infectionsUL_B * self.parameters.lowRiskcaseHospitalization_B
	
	self.hospitalizationsUH_H1 = self.infectionsUH_H1 * self.parameters.highRiskcaseHospitalization_H1
	self.hospitalizationsUH_H3 = self.infectionsUH_H3 * self.parameters.highRiskcaseHospitalization_H3
	self.hospitalizationsUH_B = self.infectionsUH_B * self.parameters.highRiskcaseHospitalization_B
	
	self.hospitalizationsVL_H1 = self.infectionsUL_H1 * (1 - self.parameters.vaccineEfficacyVsHospitalization_H1) * self.parameters.lowRiskcaseHospitalization_H1
	self.hospitalizationsVL_H3 = self.infectionsUL_H3 *  (1 - self.parameters.vaccineEfficacyVsHospitalization_H3) *self.parameters.lowRiskcaseHospitalization_H3
	self.hospitalizationsVL_B = self.infectionsUL_B *  (1 - self.parameters.vaccineEfficacyVsHospitalization_B) *self.parameters.lowRiskcaseHospitalization_B
	
	self.hospitalizationsVH_H1 = self.infectionsUH_H1 * (1 - self.parameters.vaccineEfficacyVsHospitalization_H1) *self.parameters.highRiskcaseHospitalization_H1
	self.hospitalizationsVH_H3 = self.infectionsUH_H3 * (1 - self.parameters.vaccineEfficacyVsHospitalization_H3) *self.parameters.highRiskcaseHospitalization_H3
	self.hospitalizationsVH_B = self.infectionsUH_B * (1 - self.parameters.vaccineEfficacyVsHospitalization_B) *self.parameters.highRiskcaseHospitalization_B
	
	
	
	self.hospitalizationsL  = self.hospitalizationsUL_H1 + self.hospitalizationsUL_H3 + self.hospitalizationsUL_B + self.hospitalizationsVL_H1 + self.hospitalizationsVL_H3 + self.hospitalizationsVL_B
	self.hospitalizationsH  = self.hospitalizationsUH_H1 + self.hospitalizationsUH_H3 + self.hospitalizationsUH_B + self.hospitalizationsVH_H1 + self.hospitalizationsVH_H3 + self.hospitalizationsVH_B
        self.hospitalizations = self.hospitalizationsL + self.hospitalizationsH
	self.totalHospitalizations = self.hospitalizations.sum()
	
	
	self.deathsUL_H1 =   self.DUL_H1[-1,: ]
	self.deathsUL_H3 =   self.DUL_H3[-1,: ]
	self.deathsUL_B  =   self.DUL_B[-1,: ]
	
	self.deathsVL_H1 =   self.DTL_H1[-1,: ] + self.DNL_H1[-1,: ]
	self.deathsVL_H3 =   self.DTL_H3[-1,: ] + self.DNL_H3[-1,: ]
	self.deathsVL_B  =   self.DTL_B[-1,: ] + self.DNL_B[-1,: ]
	
	self.deathsUH_H1 =   self.DUH_H1[-1,: ]
	self.deathsUH_H3 =   self.DUH_H3[-1,: ]
	self.deathsUH_B  =   self.DUH_B[-1,: ]
	
	self.deathsVH_H1 =   self.DTH_H1[-1,: ] + self.DNH_H1[-1,: ]
	self.deathsVH_H3 =   self.DTH_H3[-1,: ] + self.DNH_H3[-1,: ]
	self.deathsVH_B  =   self.DTH_B[-1,: ] + self.DNH_B[-1,: ]
	
	
        # Find duplicate times: these are where vaccination occurs
        #for i in numpy.arange(len(self.T)).compress(numpy.diff(self.T) == 0):
            # Update for vaccinated people
        #    self.deathsUL += self.NUL[i + 1, :] - self.NUL[i, :]
	#    self.deathsUH += self.NUH[i + 1, :] - self.NUH[i, :]
        #    self.deathsTL += self.NTL[i + 1, :] - self.NTL[i, :]
	#    self.deathsTH += self.NTH[i + 1, :] - self.NTH[i, :]
	#    self.deathsNL += self.NNL[i + 1, :] - self.NNL[i, :]
	#    self.deathsNH += self.NNH[i + 1, :] - self.NNH[i, :]


	self.deathsUL = self.deathsUL_H1 + self.deathsUL_H3 + self.deathsUL_B
	self.deathsUH = self.deathsUH_H1 + self.deathsUH_H3 + self.deathsUH_B
	self.deathsVL = self.deathsVL_H1 + self.deathsVL_H3 + self.deathsVL_B 
	self.deathsVH = self.deathsVH_H1 + self.deathsVH_H3 + self.deathsVH_B 
        self.deaths   = self.deathsUL + self.deathsUH + self.deathsVL + self.deathsVH
        self.totalDeaths = self.deaths.sum()

	self.YLL = numpy.multiply(self.parameters.expectationOfLife, self.deaths)

	#Years lived with disability
	################################################################################
	### Complications cases that are hospitalized
	self._ARDS = numpy.multiply(self.infections, self.parameters.caseARDSfraction)
	self.YLD_ARDS = self.parameters.disabilityWeightARDS * (numpy.multiply(self._ARDS, self.parameters.expectationOfLife))
	
	self._pneumonia = numpy.multiply(self.infections, self.parameters.casePneumoniafraction)
	self.YLD_pneumonia = self.parameters.disabilityWeightPneumonia * (numpy.multiply(self._pneumonia, self.parameters.durationPneumonia))

	self._HospitalizedComplicated_cases = self._ARDS + self._pneumonia
	self.YLD_HospitalizedComplicated = self.YLD_ARDS + self.YLD_pneumonia
	################################################################################
	## Complications that are not hospitalized
	
	self._otitis = numpy.multiply(self.infections, self.parameters.caseOtitisfraction)
	self.YLD_otitis = self.parameters.disabilityWeightOtitis * (numpy.multiply(self._otitis, self.parameters.durationOtitis))
	
	
	self._deafness = numpy.multiply(self.infections, self.parameters.caseDeafnessfraction)
	self.YLD_deafness = numpy.multiply(self.parameters.disabilityWeightDeafness, (numpy.multiply(self._deafness, self.parameters.expectationOfLife)))
	
	self._NonhospitalizedComplicated_cases =  self._otitis + self._deafness
	self.YLD_NonhospitalizedComplicated = self.YLD_otitis + self.YLD_deafness
	################################################################################
	## Hospitalizations with no complications
	
	self._HospitalizedUncomplicated_cases = self.hospitalizations - self._HospitalizedComplicated_cases
	self.YLD_HospitalizedUncomplicated =  self.parameters.disabilityWeightHospitalizedUncomplicated * numpy.multiply(self._HospitalizedUncomplicated_cases, (1./(self.parameters.recoveryRate*365)))
	
	################################################################################
	## No hospitalizations with no complications 
	self._NonhospitalizedUncomplicated_cases = self.infections - (self._HospitalizedUncomplicated_cases + self._HospitalizedComplicated_cases + self._NonhospitalizedComplicated_cases)	
	self.YLD_NonhospitalizedUncomplicated = self.parameters.disabilityWeightUncomplicated * numpy.multiply(self._NonhospitalizedUncomplicated_cases, (1./(self.parameters.recoveryRate*365)))
	
	########################################################################
	##total YLD
	self.YLD = self.YLD_NonhospitalizedUncomplicated + self.YLD_HospitalizedUncomplicated + self.YLD_HospitalizedComplicated + self.YLD_NonhospitalizedComplicated
	

	self.DALY = (self.YLL + self.YLD)
	self.totalDALY = self.DALY.sum()
	

        
        
    def simulate(self):
        self.updateIC()
        self.solve()
        self.updateStats()

    def updateProportionVaccinated(self, PVPWVal, nVacTypes, vacDoses):
        # Update propotion vaccinated

	 # Convert flat vector to 2-D array
        if numpy.ndim(PVPWVal) != 2:
	  
            PVPWVal = numpy.asarray(PVPWVal).reshape(
                (nVacTypes,
                 self.parameters.proportionVaccinatedLLength + self.parameters.proportionVaccinatedHLength))

	
	## raw coverage among age groups
	proportionVaccinatedTL = PVPWVal[0][: self.parameters.proportionVaccinatedLLength]
	proportionVaccinatedTH = list(PVPWVal[0][self.parameters.proportionVaccinatedLLength:])*len(proportionVaccinatedTL)
	proportionVaccinatedNL = PVPWVal[1][: self.parameters.proportionVaccinatedLLength]
	proportionVaccinatedNH = list(PVPWVal[1][self.parameters.proportionVaccinatedLLength:])*len(proportionVaccinatedNL)
	
	
	#raw dose uptake among age groups
	dosesVaccinatedTL =  [(a*b) for a, b in zip(proportionVaccinatedTL, self.parameters.population_lowrisk[1:])]
	dosesVaccinatedTH =  [(a*b) for a, b in zip(proportionVaccinatedTH, self.parameters.population_highrisk[1:])]
	dosesVaccinatedNL =  [(a*b) for a, b in zip(proportionVaccinatedNL, self.parameters.population_lowrisk[1:])]
	dosesVaccinatedNH =  [(a*b) for a, b in zip(proportionVaccinatedNH, self.parameters.population_highrisk[1:])]
	
	
	
	##coverage RELATIVE to each other (based is age-group 0 for Low risk age groups)
	relative_coverage_TL = [num/dosesVaccinatedTL[0] for num in dosesVaccinatedTL]
	relative_coverage_TH = [num/dosesVaccinatedTL[0] for num in dosesVaccinatedTH]
	relative_coverage_NL = [num/dosesVaccinatedNL[0] for num in dosesVaccinatedNL]
	relative_coverage_NH = [num/dosesVaccinatedNL[0] for num in dosesVaccinatedNH]
	
	#multiplication factor for final doses
	xFac_T = vacDoses[0]/(1.*(sum(relative_coverage_TL)+ sum(relative_coverage_TH)))
	xFac_N = vacDoses[1]/(1.*(sum(relative_coverage_NL)+ sum(relative_coverage_TH)))
	
	    
	doses_TL = [xFac_T * num for num in relative_coverage_TL]
	doses_TH = [xFac_T * num for num in relative_coverage_TH]
	doses_NL = [xFac_N * num for num in relative_coverage_NL]
	doses_NH = [xFac_N * num for num in relative_coverage_NH]
	
	    
        self.parameters.proportionVaccinatedTLPW.values = [(a/(1.*b)) for a,b in zip(doses_TL,  self.parameters.population_lowrisk[1:])]
	self.parameters.proportionVaccinatedTHPW.values = [(a/(1.*b)) for a,b in zip(doses_TH,  self.parameters.population_highrisk[1:])]
	self.parameters.proportionVaccinatedNLPW.values = [(a/(1.*b)) for a,b in zip(doses_NL,  self.parameters.population_lowrisk[1:])]
	self.parameters.proportionVaccinatedNHPW.values = [(a/(1.*b)) for a,b in zip(doses_NH,  self.parameters.population_highrisk[1:])]

	## extend to full ages groups. Proportions calculated by multiplying PVPWVal 
	##values with the matrix defined in S.160
	
	
	self.parameters.proportionVaccinatedTL = self.parameters.proportionVaccinatedTLPW.full(self.parameters.ages)
	self.parameters.proportionVaccinatedTH = self.parameters.proportionVaccinatedTHPW.full(self.parameters.ages)
	self.parameters.proportionVaccinatedNL = self.parameters.proportionVaccinatedNLPW.full(self.parameters.ages)
	self.parameters.proportionVaccinatedNH = self.parameters.proportionVaccinatedNHPW.full(self.parameters.ages)
	
	

	vacsUsedTypical = sum(doses_TL + doses_TH)
	vacsUsedUniversal = sum(doses_NL + doses_NH)

	   
        # Update initial condition for ODEs
        self.updateIC()

        return vacsUsedTypical, vacsUsedUniversal
    
    def vaccinated_output(self):
        return list(self.parameters.proportionVaccinatedTypical), list(self.parameters.proportionVaccinatedUniversal), [(a*b) for (a,b) in zip(self.parameters.proportionVaccinatedTypical, self.parameters.population)],[(a*b) for (a,b) in zip(self.parameters.proportionVaccinatedUniversal, self.parameters.population)]
    
    
        
    def simulateWithVaccine(self, PVPWVals, vacEfficacy, vacDoses):
	
        nVacTypes = len(vacEfficacy)

        self.resetSolution()

        # Vaccinate the population
	vacsUsedTypical, vacsUsedUniversal = self.updateProportionVaccinated(PVPWVals, nVacTypes, vacDoses)
	
	if vacsUsedUniversal <0 and min(list(PVPWVals)) >0 : print ("check!!!!"), PVPWVals, nVacTypes, vacsUsedUniversal
        tEnd = self.tMax
	tStart = self.tMin
	self.solve(tStart = tStart, tEnd = tEnd)
	self.updateStats()
	
	unvax  = ((self.NUL[0,:].sum() + self.NUH[0,:].sum() - self.SUL[-1,:].sum() - self.SUH[-1,:].sum()))/1e6
	typical = (self.NTL[0,:].sum() + self.NTH[0,:].sum() - self.STL[-1,:].sum() - self.STH[-1,:].sum())/1e6
	universal = ((self.NNL[0,:].sum() + self.NNH[0,:].sum() - self.SNL[-1,:].sum() - self.SNH[-1,:].sum()))/1e6
	print ("Unvax infected"), unvax
	print ("Typical infected"), typical
	print ("Universal infected"), universal
	print ("total"), unvax + typical + universal
	print ("check total"), (sum(self.infections_H1) + sum(self.infections_H3) + sum(self.infections_B))/1e6
	print ("infections H1"), (sum(self.infections_H1) *100.)/(1e6*(unvax + typical + universal))
	print ("infections H3"), (sum(self.infections_H3) *100.)/(1e6*(unvax + typical + universal))
	print ("infections B"), (sum(self.infections_B) *100.)/(1e6 * (unvax + typical + universal))

	#import matplotlib.pyplot as plt
	#times = [num for num in xrange(tEnd+1)]
	#plt.plot(times, (self.infections_H1).sum(axis=1), color = "red", linestyle= "-")
	#plt.plot(times, (self.infections_H3).sum(axis=1), color = "blue", linestyle= "-")
	#plt.plot(times, (self.infections_B).sum(axis=1), color = "green", linestyle = "-")

	
		 
	#plt.show()


        self.updateStats()

	return vacsUsedTypical, vacsUsedUniversal, self.parameters.proportionVaccinatedTypical, self.parameters.proportionVaccinatedUniversal 
    
    def outputInfo(self):
        print 'R0:\t\t\t %g' % self.parameters.R0
        print ('Infections (in millions):'),  self.totalInfections/1000000.
        print 'Deaths:\t\t\t %g' % self.totalDeaths
        print 'Hospitalizations:\t %g' % self.totalHospitalizations
        #print ('Age-specific infections:'), list(self.infections)
	print ("unvaccinated"),self.infectionsU.sum()
	print ("vaccinated"),self.infectionsV.sum()
	print ("total pop size"), self.parameters.population.sum()
	print ("total unvaccinated"),  ((1 - self.parameters.proportionVaccinatedTypical -  self.parameters.proportionVaccinatedUniversal) * self.parameters.population).sum()
	print ("total vaccinated Typical"),  ((self.parameters.proportionVaccinatedTypical) * self.parameters.population).sum()
	print ("total vaccinated Universal"),  ((self.parameters.proportionVaccinatedUniversal) * self.parameters.population).sum()


    def optimization_output(self):
	return self.parameters.proportionVaccinatedTypical ,  self.parameters.proportionVaccinatedUniversal


    def short_output(self):
	return list(self.infections), list(self.hospitalizations), list(self.deaths), list(self.DALY)

    def debug_info(self):
	return self.infectionsU.sum()
	

    #def vaccinated_output(self):
    #    return self.parameters.population.sum(),((1 - self.parameters.proportionVaccinated) * self.parameters.population).sum(), ((self.parameters.proportionVaccinated) * self.parameters.population).sum(), self.infectionsU.sum(), self.infectionsV.sum()


