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
        self.Y0 = numpy.zeros(84 * self.parameters.ages.size)

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
        return (self.SUL[-1, :], self.IUL_H1[-1, :], self.IUL_H3[-1, :], self.IUL_B[-1, :], self.RUL_H1[-1, :], self.RUL_H3[-1, :], self.RUL_B[-1, :], self.IUL_H1_H3[-1, :], self.IUL_H1_B[-1, :], self.IUL_H3_H1[-1, :], self.IUL_H3_B[-1, :], self.IUL_B_H3[-1, :], self.IUL_B_H1[-1, :], self.RUL[-1, :], 
		self.SUH[-1, :], self.IUH_H1[-1, :], self.IUH_H3[-1, :], self.IUH_B[-1, :], self.RUH_H1[-1, :], self.RUH_H3[-1, :], self.RUH_B[-1, :], self.IUH_H1_H3[-1, :], self.IUH_H1_B[-1, :], self.IUH_H3_H1[-1, :], self.IUH_H3_B[-1, :], self.IUH_B_H3[-1, :], self.IUH_B_H1[-1, :], self.RUH[-1, :], 
		self.STL[-1, :], self.ITL_H1[-1, :], self.ITL_H3[-1, :], self.ITL_B[-1, :], self.RTL_H1[-1, :], self.RTL_H3[-1, :], self.RTL_B[-1, :], self.ITL_H1_H3[-1, :], self.ITL_H1_B[-1, :], self.ITL_H3_H1[-1, :], self.ITL_H3_B[-1, :], self.ITL_B_H3[-1, :], self.ITL_B_H1[-1, :], self.RTL[-1, :], 
		self.STH[-1, :], self.ITH_H1[-1, :], self.ITH_H3[-1, :], self.ITH_B[-1, :], self.RTH_H1[-1, :], self.RTH_H3[-1, :], self.RTH_B[-1, :], self.ITH_H1_H3[-1, :], self.ITH_H1_B[-1, :], self.ITH_H3_H1[-1, :], self.ITH_H3_B[-1, :], self.ITH_B_H3[-1, :], self.ITH_B_H1[-1, :], self.RTH[-1, :], 
		self.SNL[-1, :], self.INL_H1[-1, :], self.INL_H3[-1, :], self.INL_B[-1, :], self.RNL_H1[-1, :], self.RNL_H3[-1, :], self.RNL_B[-1, :], self.INL_H1_H3[-1, :], self.INL_H1_B[-1, :], self.INL_H3_H1[-1, :], self.INL_H3_B[-1, :], self.INL_B_H3[-1, :], self.INL_B_H1[-1, :], self.RNL[-1, :], 
		self.SNH[-1, :], self.INH_H1[-1, :], self.INH_H3[-1, :], self.INH_B[-1, :], self.RNH_H1[-1, :], self.RNH_H3[-1, :], self.RNH_B[-1, :], self.INH_H1_H3[-1, :], self.INH_H1_B[-1, :], self.INH_H3_H1[-1, :], self.INH_H3_B[-1, :], self.INH_B_H3[-1, :], self.INH_B_H1[-1, :], self.RNH[-1, :])
    
    def updateIC(self):
        if not self.hasSolution:
            # S
	    ## SUL
            self.Y0[ 0: : 84] =  (1 - self.parameters.proportionVaccinatedTL -  self.parameters.proportionVaccinatedNL) * self.parameters.population_lowrisk 

	    ## SUH
            self.Y0[14: : 84] =  (1 - self.parameters.proportionVaccinatedTH -  self.parameters.proportionVaccinatedNH) * self.parameters.population_highrisk 
  
	    ## STL
            self.Y0[28: : 84] = self.parameters.proportionVaccinatedTL * self.parameters.population_lowrisk 
	    
	    ## STH
            self.Y0[42: : 84] =  self.parameters.proportionVaccinatedTH * self.parameters.population_highrisk

	    ## SNL
            self.Y0[56: : 84] = self.parameters.proportionVaccinatedNL * self.parameters.population_lowrisk 
	    ## SNH 
            self.Y0[70: : 84] = self.parameters.proportionVaccinatedNH  * self.parameters.population_highrisk
	
            # I: Add a single infectious person in each age
	    # IUL_H1, IUL_H3, IUL_B
	    self.Y0[ 1: : 84] = numpy.full(self.parameters.ages.size,1000)
	    self.Y0[ 2: : 84] = numpy.full(self.parameters.ages.size,1000)
	    self.Y0[ 3: : 84] = numpy.full(self.parameters.ages.size,1000)
	    
	    # IUH_H1, IUH_H3, IUH_B
            self.Y0[ 15: : 84] = numpy.full(self.parameters.ages.size, 1000)
	    self.Y0[ 16: : 84] = numpy.full(self.parameters.ages.size, 1000)
	    self.Y0[ 17: : 84] = numpy.full(self.parameters.ages.size, 1000)
	    
	    # ITL_H1, ITL_H3, ITL_B
	    self.Y0[ 29: : 84] = numpy.full(self.parameters.ages.size, 1000)
	    self.Y0[ 30: : 84] = numpy.full(self.parameters.ages.size, 1000)
	    self.Y0[ 31: : 84] = numpy.full(self.parameters.ages.size, 1000)
	    
	    # ITH_H1, ITH_H3, ITH_B
	    self.Y0[43: : 84] = numpy.full(self.parameters.ages.size, 1000)
	    self.Y0[44: : 84] = numpy.full(self.parameters.ages.size, 1000)
	    self.Y0[45: : 84] = numpy.full(self.parameters.ages.size, 1000)
	    
	    # INL_H1, INL_H3, INL_B
	    self.Y0[57: : 84] = numpy.full(self.parameters.ages.size, 1000)
	    self.Y0[58: : 84] = numpy.full(self.parameters.ages.size, 1000)
	    self.Y0[59: : 84] = numpy.full(self.parameters.ages.size, 1000)
	    
	    # INH_H1, INH_H3, INH_B
	    self.Y0[71: : 84] = numpy.full(self.parameters.ages.size, 1000)
	    self.Y0[72: : 84] = numpy.full(self.parameters.ages.size, 1000)
	    self.Y0[73: : 84] = numpy.full(self.parameters.ages.size, 1000)

            # S: Remove those new infectious people from the susceptibles
            self.Y0[ 0: : 84] -= (self.Y0[ 1: : 84] + self.Y0[ 2: : 84] + self.Y0[ 3: : 84])
	    self.Y0[ 14: : 84] -= (self.Y0[ 15: : 84] + self.Y0[ 16: : 84] + self.Y0[ 17: : 84])
	    self.Y0[ 28: : 84] -= (self.Y0[ 29: : 84] + self.Y0[ 30: : 84] + self.Y0[ 31: : 84])
	    self.Y0[ 42: : 84] -= (self.Y0[ 43: : 84] + self.Y0[ 44: : 84] + self.Y0[ 45: : 84])
	    self.Y0[ 56: : 84] -= (self.Y0[ 57: : 84] + self.Y0[ 58: : 84] + self.Y0[ 59: : 84])
	    self.Y0[ 70: : 84] -= (self.Y0[ 71: : 84] + self.Y0[ 72: : 84] + self.Y0[ 73: : 84])
 
            # R, I reinfections ans final R
	    for num in range(4,14):
		self.Y0[ num:  : 84] = 0.
	    
	    for num in range(18, 28):
		self.Y0[ num:  : 84] = 0.
	    
	    for num in range(32, 42):
		self.Y0[ num:  : 84] = 0.
		
	    for num in range(46, 56):
		self.Y0[ num:  : 84] = 0.
		
	    for num in range(60, 70):
		self.Y0[ num:  : 84] = 0.
		
	    for num in range(74, 84):
		self.Y0[ num:  : 84] = 0.

	   
        else:
	    SUL, IUL_H1, IUL_H3, IUL_B, RUL_H1, RUL_H3, RUL_B, IUL_H1_H3 , IUL_H1_B, IUL_H3_H1, IUL_H3_B, IUL_B_H3, IUL_B_H1, RUL, 
	    SUH, IUH_H1, IUH_H3, IUH_B, RUH, RUH_H3, RUH_B, IUH_H1_H3 , IUH_H1_B, IUH_H3_H1, IUH_H3_B, IUH_B_H3, IUH_B_H1, RUH, 
	    STL, ITL_H1, ITL_H3, ITL_B, RTL, RTL_H3, RTL_B,ITL_H1_H3 , ITL_H1_B, ITL_H3_H1, ITL_H3_B, ITL_B_H3, ITL_B_H1, RTL, 
	    STH, ITH_H1, ITH_H3, ITH_B, RTH, RTH_H3, RTH_B,ITH_H1_H3 , ITH_H1_B, ITH_H3_H1, ITH_H3_B, ITH_B_H3, ITH_B_H1, RTH, 
	    SNL, INL_H1, INL_H3, INL_B, RNL, RNL_H3, RNL_B,INL_H1_H3 , INL_H1_B, INL_H3_H1, INL_H3_B, INL_B_H3, INL_B_H1, RNL, 
	    SNH, INH_H1, INH_H3, INH_B, RNH, RNH_H3, RNH_B, INH_H1_H3 , INH_H1_B, INH_H3_H1, INH_H3_B, INH_B_H3, INH_B_H1, RNH   = self.getLastValues()
            
            self.Y0[ 0 : : 84] = (1 - self.parameters.proportionVaccinatedL) * SUL
	    self.Y0[ 14 : : 84] = (1 - self.parameters.proportionVaccinatedH) * SUH
	    self.Y0[ 28 : : 84] = STL + (1 - self.parameters.proportionVaccinatedTL) * SUL
	    self.Y0[ 42 : : 84] = STH + (1 - self.parameters.proportionVaccinatedTH) * SUH
	    self.Y0[ 56 : : 84] = SNL + (1 - self.parameters.proportionVaccinatedNL) * SUL
	    self.Y0[ 70 : : 84] = SNH + (1 - self.parameters.proportionVaccinatedNH) * SUH
	    
	    #I
	    self.Y0[ 1 : : 84] = IUL_H1
	    self.Y0[ 2 : : 84] = IUL_H3
	    self.Y0[ 3 : : 84] = IUL_B
	    
	    self.Y0[ 15 : : 84] = IUH_H1
	    self.Y0[ 16 : : 84] = IUH_H3
	    self.Y0[ 17 : : 84] = IUH_B
	    
	    self.Y0[ 29 : : 84] = ITL_H1
	    self.Y0[ 30 : : 84] = ITL_H3
	    self.Y0[ 31 : : 84] = ITL_B
	    
	    self.Y0[ 43 : : 84] = ITH_H1
	    self.Y0[ 44 : : 84] = ITH_H3
	    self.Y0[ 45 : : 84] = ITH_B
	    
	    self.Y0[ 57 : : 84] = INL_H1
	    self.Y0[ 58 : : 84] = INL_H3
	    self.Y0[ 59 : : 84] = INL_B
	    
	    self.Y0[ 71 : : 84] = INH_H1
	    self.Y0[ 72 : : 84] = INH_H3
	    self.Y0[ 73 : : 84] = INH_B
	    
	    #R
	    self.Y0[ 4 : : 84] = RUL_H1
	    self.Y0[ 5 : : 84] = RUL_H3
	    self.Y0[ 6 : : 84] = RUL_B
	    
	    self.Y0[ 18 : : 84] = RUH_H1
	    self.Y0[ 19 : : 84] = RUH_H3
	    self.Y0[ 20 : : 84] = RUH_B
	    
	    self.Y0[ 32 : : 84] = RTL_H1
	    self.Y0[ 33 : : 84] = RTL_H3
	    self.Y0[ 34 : : 84] = RTL_B
	    
	    self.Y0[ 46 : : 84] = RTH_H1
	    self.Y0[ 47 : : 84] = RTH_H3
	    self.Y0[ 48 : : 84] = RTH_B
	    
	    self.Y0[ 60 : : 84] = RNL_H1
	    self.Y0[ 61 : : 84] = RNL_H3
	    self.Y0[ 62 : : 84] = RNL_B
	    
	    self.Y0[ 74 : : 84] = RNH_H1
	    self.Y0[ 75 : : 84] = RNH_H3
	    self.Y0[ 76 : : 84] = RNH_B
	    
	    #re-infected T class
	    self.Y0[ 7 : : 84] = IUL_H1_H3
	    self.Y0[ 8 : : 84] = IUL_H1_B
	    self.Y0[ 9 : : 84] = IUL_H3_H1
	    self.Y0[ 10 : : 84] = IUL_H3_B
	    self.Y0[ 11 : : 84] = IUL_B_H3
	    self.Y0[ 12 : : 84] = IUL_B_H1
	    
	    self.Y0[ 7 : : 84] = IUH_H1_H3
	    self.Y0[ 8 : : 84] = IUH_H1_B
	    self.Y0[ 9 : : 84] = IUH_H3_H1
	    self.Y0[ 10 : : 84] = IUH_H3_B
	    self.Y0[ 11 : : 84] = IUH_B_H3
	    self.Y0[ 12 : : 84] = IUH_B_H1
	    
	    self.Y0[ 7 : : 84] = ITL_H1_H3
	    self.Y0[ 8 : : 84] = ITL_H1_B
	    self.Y0[ 9 : : 84] = ITL_H3_H1
	    self.Y0[ 10 : : 84] = ITL_H3_B
	    self.Y0[ 11 : : 84] = ITL_B_H3
	    self.Y0[ 12 : : 84] = ITL_B_H1
	    
	    self.Y0[ 7 : : 84] = ITH_H1_H3
	    self.Y0[ 8 : : 84] = ITH_H1_B
	    self.Y0[ 9 : : 84] = ITH_H3_H1
	    self.Y0[ 10 : : 84] = ITH_H3_B
	    self.Y0[ 11 : : 84] = ITH_B_H3
	    self.Y0[ 12 : : 84] = ITH_B_H1
	    
	    self.Y0[ 7 : : 84] = INL_H1_H3
	    self.Y0[ 8 : : 84] = INL_H1_B
	    self.Y0[ 9 : : 84] = INL_H3_H1
	    self.Y0[ 10 : : 84] = INL_H3_B
	    self.Y0[ 11 : : 84] = INL_B_H3
	    self.Y0[ 12 : : 84] = INL_B_H1
	    
	    self.Y0[ 7 : : 84] = INH_H1_H3
	    self.Y0[ 8 : : 84] = INH_H1_B
	    self.Y0[ 9 : : 84] = INH_H3_H1
	    self.Y0[ 10 : : 84] = INH_H3_B
	    self.Y0[ 11 : : 84] = INH_B_H3
	    self.Y0[ 12 : : 84] = INH_B_H1
	    
	    #R class
	    self.Y0[ 13 : : 84] = RUL
	    self.Y0[ 27 : : 84] = RUH
	    self.Y0[ 41 : : 84] = RTL
	    self.Y0[ 55 : : 84] = RTH
	    self.Y0[ 69 : : 84] = RNL
	    self.Y0[ 83 : : 84] = RNH


    def RHS(self, Y, t):
        '''
        SIR model with multiple host types.
        
        This function gives the right-hand sides of the ODEs.
        '''
        
        # Convert vector to meaningful component vectors

	SUL    = Y[ 0 : : 84]
        IUL_H1 = Y[ 1 : : 84]
	IUL_H3 = Y[ 2 : : 84]
	IUL_B  = Y[ 3 : : 84]
	RUL_H1 = Y[ 4 : : 84]
	RUL_H3 = Y[ 5 : : 84]
	RUL_B = Y[ 6 : : 84]
	IUL_H1_H3 = Y[ 7 : : 84]
	IUL_H1_B = Y[ 8 : : 84]
	IUL_H3_H1 = Y[ 9 : : 84]
	IUL_H3_B = Y[ 10 : : 84]
	IUL_B_H3 = Y[ 11 : : 84]
	IUL_B_H1 = Y[ 12 : : 84]
	RUL =  Y[13: : 84]
	
	SUH    = Y[ 14 : : 84]
        IUH_H1 = Y[ 15 : : 84]
	IUH_H3 = Y[ 16 : : 84]
	IUH_B  = Y[ 17 : : 84]
	RUH_H1 = Y[ 18 : : 84]
	RUH_H3 = Y[ 19 : : 84]
	RUH_B = Y[ 20 : : 84]
	IUH_H1_H3 = Y[ 21 : : 84]
	IUH_H1_B = Y[ 22 : : 84]
	IUH_H3_H1 = Y[ 23 : : 84]
	IUH_H3_B = Y[ 24 : : 84]
	IUH_B_H3 = Y[ 25 : : 84]
	IUH_B_H1 = Y[ 26 : : 84]
	RUH = Y[27: : 84]
	
	STL    = Y[ 28 : : 84]
        ITL_H1 = Y[ 29 : : 84]
	ITL_H3 = Y[ 30 : : 84]
	ITL_B  = Y[ 31 : : 84]
	RTL_H1 = Y[ 32 : : 84]
	RTL_H3 = Y[ 33 : : 84]
	RTL_B = Y[ 34 : : 84]
	ITL_H1_H3 = Y[ 35 : : 84]
	ITL_H1_B = Y[ 36 : : 84]
	ITL_H3_H1 = Y[ 37 : : 84]
	ITL_H3_B = Y[ 38 : : 84]
	ITL_B_H3 = Y[ 39 : : 84]
	ITL_B_H1 = Y[ 40 : : 84]
	RTL = Y[41: : 84]
	
	STH    = Y[ 42 : : 84]
        ITH_H1 = Y[ 43 : : 84]
	ITH_H3 = Y[ 44 : : 84]
	ITH_B  = Y[ 45 : : 84]
	RTH_H1 = Y[ 46 : : 84]
	RTH_H3 = Y[ 47 : : 84]
	RTH_B = Y[ 48 : : 84]
	ITH_H1_H3 = Y[ 49 : : 84]
	ITH_H1_B = Y[ 50 : : 84]
	ITH_H3_H1 = Y[ 51 : : 84]
	ITH_H3_B = Y[ 52 : : 84]
	ITH_B_H3 = Y[ 53 : : 84]
	ITH_B_H1 = Y[ 54 : : 84]
	RTH = Y[55: : 84]
	
	SNL    = Y[ 56 : : 84]
        INL_H1 = Y[ 57 : : 84]
	INL_H3 = Y[ 58 : : 84]
	INL_B  = Y[ 59 : : 84]
	RNL_H1 = Y[ 60 : : 84]
	RNL_H3 = Y[ 61 : : 84]
	RNL_B = Y[ 62 : : 84]
	INL_H1_H3 = Y[ 63 : : 84]
	INL_H1_B = Y[ 64 : : 84]
	INL_H3_H1 = Y[ 65 : : 84]
	INL_H3_B = Y[ 66 : : 84]
	INL_B_H3 = Y[ 67 : : 84]
	INL_B_H1 = Y[ 68 : : 84]
	RNL = Y[69: : 84]
	
	SNH    = Y[ 70 : : 84]
        INH_H1 = Y[ 71 : : 84]
	INH_H3 = Y[ 72 : : 84]
	INH_B  = Y[ 73 : : 84]
	RNH_H1 = Y[ 74 : : 84]
	RNH_H3 = Y[ 75 : : 84]
	RNH_B = Y[ 76 : : 84]
	INH_H1_H3 = Y[ 77 : : 84]
	INH_H1_B = Y[ 78 : : 84]
	INH_H3_H1 = Y[ 79 : : 84]
	INH_H3_B = Y[ 80 : : 84]
	INH_B_H3 = Y[ 81 : : 84]
	INH_B_H1 = Y[ 82 : : 84]
	RNH = Y[83: : 84]
	
        N =  sum(SUL+ IUL_H1+ IUL_H3+ IUL_B+ RUL_H1+ RUL_H3+ RUL_B+ IUL_H1_H3 + IUL_H1_B+ IUL_H3_H1+ IUL_H3_B+ IUL_B_H3+ IUL_B_H1+ RUL+ 
	    SUH+ IUH_H1+ IUH_H3+ IUH_B+ RUH+ RUH_H3+ RUH_B+ IUH_H1_H3 + IUH_H1_B+ IUH_H3_H1+ IUH_H3_B+ IUH_B_H3+ IUH_B_H1+ RUH+ 
	    STL+ ITL_H1+ ITL_H3+ ITL_B+ RTL+ RTL_H3+ RTL_B+ITL_H1_H3 + ITL_H1_B+ ITL_H3_H1+ ITL_H3_B+ ITL_B_H3+ ITL_B_H1+ RTL+ 
	    STH+ ITH_H1+ ITH_H3+ ITH_B+ RTH+ RTH_H3+ RTH_B+ITH_H1_H3 + ITH_H1_B+ ITH_H3_H1+ ITH_H3_B+ ITH_B_H3+ ITH_B_H1+ RTH+ 
	    SNL+ INL_H1+ INL_H3+ INL_B+ RNL+ RNL_H3+ RNL_B+INL_H1_H3 + INL_H1_B+ INL_H3_H1+ INL_H3_B+ INL_B_H3+ INL_B_H1+ RNL+ 
	    SNH+ INH_H1+ INH_H3+ INH_B+ RNH+ RNH_H3+ RNH_B+ INH_H1_H3 + INH_H1_B+ INH_H3_H1+ INH_H3_B+ INH_B_H3+ INH_B_H1+ RNH) 
      
        # The force of infection
        Lambda_H1 = self.parameters.transmissionScaling_H1 * self.parameters.susceptibility_H1\
		    * numpy.dot(self.parameters.contactMatrix, self.parameters.transmissibility * (IUL_H1 + IUL_H1_H3 + IUL_H1_B + IUH_H1 + IUH_H1_H3 + IUH_H1_B+  ITL_H1 + ITL_H1_H3 + ITL_H1_B+ ITH_H1+ ITH_H1_H3 + ITH_H1_B + INL_H1+ INL_H1_H3 + INL_H1_B + INH_H1+ INH_H1_H3 + INH_H1_B)) / N
	
	Lambda_H3 = self.parameters.transmissionScaling_H3 * self.parameters.susceptibility_H3 \
                 * numpy.dot(self.parameters.contactMatrix, 
                             self.parameters.transmissibility * (IUL_H3 + IUL_H3_H1 + IUL_H3_B + IUH_H3 + IUH_H3_H1 + IUH_H3_B+  ITL_H3 + ITL_H3_H1 + ITL_H3_B+ ITH_H3+ ITH_H3_H1 + ITH_H3_B + INL_H3+ INL_H3_H1 + INL_H3_B + INH_H3+ INH_H3_H1 + INH_H3_B)) / N
		
	Lambda_B = self.parameters.transmissionScaling_B * self.parameters.susceptibility_B \
                 * numpy.dot(self.parameters.contactMatrix,
                             self.parameters.transmissibility * (IUL_B + IUL_B_H3 + IUL_B_H1 + IUH_B + IUH_B_H3 + IUH_B_H1 +  ITL_B + ITL_B_H3 + ITL_B_H1+ ITH_B+ ITH_B_H3 + ITH_B_H1 + INL_B+ INL_B_H3 + INL_B_H1 + INH_B+ INH_B_H3 + INH_B_H1)) / N
	
	
        
        # The right-hand sides
	
	#UL
        dSUL    = - (Lambda_H1 + Lambda_H3 + Lambda_B) * SUL 
        dIUL_H1 = (Lambda_H1 * SUL) - (self.parameters.recoveryRate + self.parameters.deathRateUL_H1) * IUL_H1
	dIUL_H3 = (Lambda_H3 * SUL) - (self.parameters.recoveryRate + self.parameters.deathRateUL_H3) * IUL_H3
	dIUL_B  = (Lambda_B * SUL) - (self.parameters.recoveryRate + self.parameters.deathRateUL_B) * IUL_B
        dRUL_H1 = self.parameters.recoveryRate * IUL_H1 - ((1 - self.parameters.crossImmunity)*Lambda_H3 * RUL_H1) - ((1 - self.parameters.crossImmunity)*Lambda_B * RUL_H1)
	dRUL_H3 = self.parameters.recoveryRate * IUL_H3 - ((1 - self.parameters.crossImmunity)*Lambda_H1 * RUL_H3) - ((1 - self.parameters.crossImmunity)*Lambda_B * RUL_H3)
	dRUL_B = self.parameters.recoveryRate * IUL_B - ((1 - self.parameters.crossImmunity)*Lambda_H1 * RUL_B) - ((1 - self.parameters.crossImmunity)*Lambda_H3 * RUL_B)
	dIUL_H1_H3 = ((1 - self.parameters.crossImmunity)*Lambda_H1 * RUL_H3)  - (self.parameters.recoveryRate + self.parameters.deathRateUL_H1) * IUL_H1_H3
	dIUL_H1_B =  ((1 - self.parameters.crossImmunity)*Lambda_H1 * RUL_B) - (self.parameters.recoveryRate + self.parameters.deathRateUL_H1) * IUL_H1_B
	dIUL_H3_H1 = ((1 - self.parameters.crossImmunity)*Lambda_H3 * RUL_H1) - (self.parameters.recoveryRate + self.parameters.deathRateUL_H3) * IUL_H3_H1
	dIUL_H3_B = ((1 - self.parameters.crossImmunity)*Lambda_H3 * RUL_B) - (self.parameters.recoveryRate + self.parameters.deathRateUL_H3) * IUL_H3_B
	dIUL_B_H3 = ((1 - self.parameters.crossImmunity)*Lambda_B * RUL_H3) - (self.parameters.recoveryRate + self.parameters.deathRateUL_B) * IUL_B_H3
	dIUL_B_H1 = ((1 - self.parameters.crossImmunity)*Lambda_B * RUL_H1) - (self.parameters.recoveryRate + self.parameters.deathRateUL_B) * IUL_B_H1
	dRUL = self.parameters.recoveryRate *(IUL_H1_H3 + IUL_H1_B + IUL_H3_H1 + IUL_H3_B + IUL_B_H3 + IUL_B_H1)

	
	#UH
        dSUH    = - (Lambda_H1 + Lambda_H3 + Lambda_B) * SUH 
        dIUH_H1 = (Lambda_H1 * SUH) - (self.parameters.recoveryRate + self.parameters.deathRateUH_H1) * IUH_H1
	dIUH_H3 = (Lambda_H3 * SUH) - (self.parameters.recoveryRate + self.parameters.deathRateUH_H3) * IUH_H3
	dIUH_B  = (Lambda_B * SUH) - (self.parameters.recoveryRate + self.parameters.deathRateUH_B) * IUH_B
	dRUH_H1 = self.parameters.recoveryRate * IUH_H1 - ((1 - self.parameters.crossImmunity)*Lambda_H3 * RUH_H1) - ((1 - self.parameters.crossImmunity)*Lambda_B * RUH_H1)
	dRUH_H3 = self.parameters.recoveryRate * IUH_H3 - ((1 - self.parameters.crossImmunity)*Lambda_H1 * RUH_H3) - ((1 - self.parameters.crossImmunity)*Lambda_B * RUH_H3)
	dRUH_B = self.parameters.recoveryRate * IUH_B - ((1 - self.parameters.crossImmunity)*Lambda_H1 * RUH_B) - ((1 - self.parameters.crossImmunity)*Lambda_H3 * RUH_B)
	dIUH_H1_H3 = ((1 - self.parameters.crossImmunity)*Lambda_H1 * RUH_H3)  - (self.parameters.recoveryRate + self.parameters.deathRateUH_H1) * IUH_H1_H3
	dIUH_H1_B =  ((1 - self.parameters.crossImmunity)*Lambda_H1 * RUH_B) - (self.parameters.recoveryRate + self.parameters.deathRateUH_H1) * IUH_H1_B
	dIUH_H3_H1 = ((1 - self.parameters.crossImmunity)*Lambda_H3 * RUH_H1) - (self.parameters.recoveryRate + self.parameters.deathRateUH_H3) * IUH_H3_H1
	dIUH_H3_B = ((1 - self.parameters.crossImmunity)*Lambda_H3 * RUH_B) - (self.parameters.recoveryRate + self.parameters.deathRateUH_H3) * IUH_H3_B
	dIUH_B_H3 = ((1 - self.parameters.crossImmunity)*Lambda_B * RUH_H3) - (self.parameters.recoveryRate + self.parameters.deathRateUH_B) * IUH_B_H3
	dIUH_B_H1 = ((1 - self.parameters.crossImmunity)*Lambda_B * RUH_H1) - (self.parameters.recoveryRate + self.parameters.deathRateUH_B) * IUH_B_H1
	dRUH = self.parameters.recoveryRate *(IUH_H1_H3 + IUH_H1_B + IUH_H3_H1 + IUH_H3_B + IUH_B_H3 + IUH_B_H1)
	
	#TL
	dSTL = - ((1 - self.parameters.vaccineEfficacyVsInfectionTypical_H1) * Lambda_H1 *STL)
	- ((1 - self.parameters.vaccineEfficacyVsInfectionTypical_H3) * Lambda_H3 *STL)
	- ((1 - self.parameters.vaccineEfficacyVsInfectionTypical_B) * Lambda_B *STL)

        dITL_H1 = ((1 - self.parameters.vaccineEfficacyVsInfectionTypical_H1) * Lambda_H1 * STL) - (self.parameters.recoveryRate + self.parameters.deathRateVL_H1) * ITL_H1
	dITL_H3 = ((1 - self.parameters.vaccineEfficacyVsInfectionTypical_H3) * Lambda_H3 * STL) - (self.parameters.recoveryRate + self.parameters.deathRateVL_H3) * ITL_H3
	dITL_B = ((1 - self.parameters.vaccineEfficacyVsInfectionTypical_B) * Lambda_B * STL) - (self.parameters.recoveryRate + self.parameters.deathRateVL_B) * ITL_B
	dRTL_H1 = self.parameters.recoveryRate * ITL_H1 - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_H3) *Lambda_H3 * RTL_H1) - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_B) *Lambda_B * RTL_H1)
	dRTL_H3 = self.parameters.recoveryRate * ITL_H3 - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_H1) *Lambda_H1 * RTL_H3) - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_B) *Lambda_B * RTL_H3)
	dRTL_B = self.parameters.recoveryRate * ITL_B - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_H1) *Lambda_H1 * RTL_B) - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_H3) *Lambda_H3 * RTL_B)
	dITL_H1_H3 = ((1 - self.parameters.crossImmunity)*Lambda_H1 * RTL_H3)  - (self.parameters.recoveryRate + self.parameters.deathRateVL_H1) * ITL_H1_H3
	dITL_H1_B =  ((1 - self.parameters.crossImmunity)*Lambda_H1 * RTL_B) - (self.parameters.recoveryRate + self.parameters.deathRateVL_H1) * ITL_H1_B
	dITL_H3_H1 = ((1 - self.parameters.crossImmunity)*Lambda_H3 * RTL_H1) - (self.parameters.recoveryRate + self.parameters.deathRateVL_H3) * ITL_H3_H1
	dITL_H3_B = ((1 - self.parameters.crossImmunity)*Lambda_H3 * RTL_B) - (self.parameters.recoveryRate + self.parameters.deathRateVL_H3) * ITL_H3_B
	dITL_B_H3 = ((1 - self.parameters.crossImmunity)*Lambda_B * RTL_H3) - (self.parameters.recoveryRate + self.parameters.deathRateVL_B) * ITL_B_H3
	dITL_B_H1 = ((1 - self.parameters.crossImmunity)*Lambda_B * RTL_H1) - (self.parameters.recoveryRate + self.parameters.deathRateVL_B) * ITL_B_H1
	dRTL = self.parameters.recoveryRate *(ITL_H1_H3 + ITL_H1_B + ITL_H3_H1 + ITL_H3_B + ITL_B_H3 + ITL_B_H1)	
	
	#TH
	dSTH = - ((1 - self.parameters.vaccineEfficacyVsInfectionTypical_H1) * Lambda_H1 *STH)
	- ((1 - self.parameters.vaccineEfficacyVsInfectionTypical_H3) * Lambda_H3 *STH)
	- ((1 - self.parameters.vaccineEfficacyVsInfectionTypical_B) * Lambda_B *STH)
	
        dITH_H1 = ((1 - self.parameters.vaccineEfficacyVsInfectionTypical_H1) *Lambda_H1 * STH) - (self.parameters.recoveryRate + self.parameters.deathRateVH_H1) * ITH_H1
	dITH_H3 = ((1 - self.parameters.vaccineEfficacyVsInfectionTypical_H3) *Lambda_H3 * STH) - (self.parameters.recoveryRate + self.parameters.deathRateVH_H3) * ITH_H3
	dITH_B = ((1 - self.parameters.vaccineEfficacyVsInfectionTypical_B) * Lambda_B * STH) - (self.parameters.recoveryRate + self.parameters.deathRateVH_B) * ITH_B
	dRTH_H1 = self.parameters.recoveryRate * ITH_H1 - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_H3) *Lambda_H3 * RTH_H1) - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_B) *Lambda_B * RTH_H1)
	dRTH_H3 = self.parameters.recoveryRate * ITH_H3 - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_H1) *Lambda_H1 * RTH_H3) - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_B) *Lambda_B * RTH_H3)
	dRTH_B = self.parameters.recoveryRate * ITH_B - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_H1) *Lambda_H1 * RTH_B) - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_H3) *Lambda_H3 * RTH_B)
	dITH_H1_H3 = ((1 - self.parameters.crossImmunity)*Lambda_H1 * RTH_H3)  - (self.parameters.recoveryRate + self.parameters.deathRateVH_H1) * ITH_H1_H3
	dITH_H1_B =  ((1 - self.parameters.crossImmunity)*Lambda_H1 * RTH_B) - (self.parameters.recoveryRate + self.parameters.deathRateVH_H1) * ITH_H1_B
	dITH_H3_H1 = ((1 - self.parameters.crossImmunity)*Lambda_H3 * RTH_H1) - (self.parameters.recoveryRate + self.parameters.deathRateVH_H3) * ITH_H3_H1
	dITH_H3_B = ((1 - self.parameters.crossImmunity)*Lambda_H3 * RTH_B) - (self.parameters.recoveryRate + self.parameters.deathRateVH_H3) * ITH_H3_B
	dITH_B_H3 = ((1 - self.parameters.crossImmunity)*Lambda_B * RTH_H3) - (self.parameters.recoveryRate + self.parameters.deathRateVH_B) * ITH_B_H3
	dITH_B_H1 = ((1 - self.parameters.crossImmunity)*Lambda_B * RTH_H1) - (self.parameters.recoveryRate + self.parameters.deathRateVH_B) * ITH_B_H1
	dRTH = self.parameters.recoveryRate *(ITH_H1_H3 + ITH_H1_B + ITH_H3_H1 + ITH_H3_B + ITH_B_H3 + ITH_B_H1)
	
	#NL
	dSNL = - ((1 - self.parameters.vaccineEfficacyVsInfectionUniversal_H1) * Lambda_H1 *SNL)
	- ((1 - self.parameters.vaccineEfficacyVsInfectionUniversal_H3) * Lambda_H3 *SNL)
	- ((1 - self.parameters.vaccineEfficacyVsInfectionUniversal_B) * Lambda_B *SNL)
	
        dINL_H1 = ((1 - self.parameters.vaccineEfficacyVsInfectionUniversal_H1) *Lambda_H1 * SNL) - (self.parameters.recoveryRate + self.parameters.deathRateVL_H1) * INL_H1
	dINL_H3 = ((1 - self.parameters.vaccineEfficacyVsInfectionUniversal_H3) *Lambda_H3 * SNL) - (self.parameters.recoveryRate + self.parameters.deathRateVL_H3) * INL_H3
	dINL_B = ((1 - self.parameters.vaccineEfficacyVsInfectionUniversal_B) * Lambda_B * SNL) - (self.parameters.recoveryRate + self.parameters.deathRateVL_B) * INL_B
	dRNL_H1 = self.parameters.recoveryRate * INL_H1 - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_H3) *Lambda_H3 * RNL_H1) - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_B) *Lambda_B * RNL_H1)
	dRNL_H3 = self.parameters.recoveryRate * INL_H3 - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_H1) *Lambda_H1 * RNL_H3) - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_B) *Lambda_B * RNL_H3)
	dRNL_B = self.parameters.recoveryRate * INL_B - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_H1) *Lambda_H1 * RNL_B) - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_H3) *Lambda_H3 * RNL_B)
	dINL_H1_H3 = ((1 - self.parameters.crossImmunity)*Lambda_H1 * RNL_H3)  - (self.parameters.recoveryRate + self.parameters.deathRateVL_H1) * INL_H1_H3
	dINL_H1_B =  ((1 - self.parameters.crossImmunity)*Lambda_H1 * RNL_B) - (self.parameters.recoveryRate + self.parameters.deathRateVL_H1) * INL_H1_B
	dINL_H3_H1 = ((1 - self.parameters.crossImmunity)*Lambda_H3 * RNL_H1) - (self.parameters.recoveryRate + self.parameters.deathRateVL_H3) * INL_H3_H1
	dINL_H3_B = ((1 - self.parameters.crossImmunity)*Lambda_H3 * RNL_B) - (self.parameters.recoveryRate + self.parameters.deathRateVL_H3) * INL_H3_B
	dINL_B_H3 = ((1 - self.parameters.crossImmunity)*Lambda_B * RNL_H3) - (self.parameters.recoveryRate + self.parameters.deathRateVL_B) * INL_B_H3
	dINL_B_H1 = ((1 - self.parameters.crossImmunity)*Lambda_B * RNL_H1) - (self.parameters.recoveryRate + self.parameters.deathRateVL_B) * INL_B_H1
	dRNL = self.parameters.recoveryRate *(INL_H1_H3 + INL_H1_B + INL_H3_H1 + INL_H3_B + INL_B_H3 + INL_B_H1)
	
	#NH
	dSNH = - ((1 - self.parameters.vaccineEfficacyVsInfectionUniversal_H1) * Lambda_H1 *SNH)
	- ((1 - self.parameters.vaccineEfficacyVsInfectionUniversal_H3) * Lambda_H3 *SNH)
	- ((1 - self.parameters.vaccineEfficacyVsInfectionUniversal_B) * Lambda_B *SNH)
        dINH_H1 = ((1 - self.parameters.vaccineEfficacyVsInfectionUniversal_H1) *Lambda_H1 * SNH) - (self.parameters.recoveryRate + self.parameters.deathRateVH_H1) * INH_H1
	dINH_H3 = ((1 - self.parameters.vaccineEfficacyVsInfectionUniversal_H3) *Lambda_H3 * SNH) - (self.parameters.recoveryRate + self.parameters.deathRateVH_H3) * INH_H3
	dINH_B = ((1 - self.parameters.vaccineEfficacyVsInfectionUniversal_B) *Lambda_B * SNH) - (self.parameters.recoveryRate + self.parameters.deathRateVH_B) * INH_B
	dRNH_H1 = self.parameters.recoveryRate * INH_H1 - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_H3) *Lambda_H3 * RNH_H1) - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_B) *Lambda_B * RNH_H1)
	dRNH_H3 = self.parameters.recoveryRate * INH_H3 - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_H1) *Lambda_H1 * RNH_H3) - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_B) *Lambda_B * RNH_H3)
	dRNH_B = self.parameters.recoveryRate * INH_B - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_H1) *Lambda_H1 * RNH_B) - ((1 - self.parameters.crossImmunity)*(1 - self.parameters.vaccineEfficacyVsInfectionTypical_H3) *Lambda_H3 * RNH_B)
	dINH_H1_H3 = ((1 - self.parameters.crossImmunity)*Lambda_H1 * RNH_H3)  - (self.parameters.recoveryRate + self.parameters.deathRateVH_H1) * INH_H1_H3
	dINH_H1_B =  ((1 - self.parameters.crossImmunity)*Lambda_H1 * RNH_B) - (self.parameters.recoveryRate + self.parameters.deathRateVH_H1) * INH_H1_B
	dINH_H3_H1 = ((1 - self.parameters.crossImmunity)*Lambda_H3 * RNH_H1) - (self.parameters.recoveryRate + self.parameters.deathRateVH_H3) * INH_H3_H1
	dINH_H3_B = ((1 - self.parameters.crossImmunity)*Lambda_H3 * RNH_B) - (self.parameters.recoveryRate + self.parameters.deathRateVH_H3) * INH_H3_B
	dINH_B_H3 = ((1 - self.parameters.crossImmunity)*Lambda_B * RNH_H3) - (self.parameters.recoveryRate + self.parameters.deathRateVH_B) * INH_B_H3
	dINH_B_H1 = ((1 - self.parameters.crossImmunity)*Lambda_B * RNH_H1) - (self.parameters.recoveryRate + self.parameters.deathRateVH_B) * INH_B_H1
	dRNH = self.parameters.recoveryRate *(INH_H1_H3 + INH_H1_B + INH_H3_H1 + INH_H3_B + INH_B_H3 + INH_B_H1)
	
	
        # Convert meaningful component vectors into a single vector
        dY = numpy.empty(Y.size, dtype = float)
        dY[ 0 : : 84] = dSUL
        dY[ 1 : : 84] = dIUL_H1
        dY[ 2 : : 84] = dIUL_H3
	dY[ 3 : : 84] = dIUL_B
        dY[ 4 : : 84] = dRUL_H1
	dY[ 5 : : 84] = dRUL_H3
	dY[ 6 : : 84] = dRUL_B
	dY[ 7 : : 84] = dIUL_H1_H3
	dY[ 8 : : 84] = dIUL_H1_B
	dY[ 9 : : 84] = dIUL_H3_H1
	dY[ 10 : : 84] = dIUL_H3_B
	dY[ 11 : : 84] = dIUL_B_H3
	dY[ 12 : : 84] = dIUL_B_H1
	dY[ 13 : : 84] = dRUL
	
	dY[ 14 : : 84] = dSUH
        dY[ 15 : : 84] = dIUH_H1
        dY[ 16 : : 84] = dIUH_H3
	dY[ 17 : : 84] = dIUH_B
        dY[ 18 : : 84] = dRUH_H1
	dY[ 19 : : 84] = dRUH_H3
	dY[ 20 : : 84] = dRUH_B
	dY[ 21 : : 84] = dIUH_H1_H3
	dY[ 22 : : 84] = dIUH_H1_B
	dY[ 23 : : 84] = dIUH_H3_H1
	dY[ 24 : : 84] = dIUH_H3_B
	dY[ 25 : : 84] = dIUH_B_H3
	dY[ 26 : : 84] = dIUH_B_H1
	dY[ 27 : : 84] = dRUH
	
	dY[ 28 : : 84] = dSTL
        dY[ 29 : : 84] = dITL_H1
        dY[ 30 : : 84] = dITL_H3
	dY[ 31 : : 84] = dITL_B
        dY[ 32 : : 84] = dRTL_H1
	dY[ 33 : : 84] = dRTL_H3
	dY[ 34 : : 84] = dRTL_B
	dY[ 35 : : 84] = dITL_H1_H3
	dY[ 36 : : 84] = dITL_H1_B
	dY[ 37 : : 84] = dITL_H3_H1
	dY[ 38 : : 84] = dITL_H3_B
	dY[ 39 : : 84] = dITL_B_H3
	dY[ 40 : : 84] = dITL_B_H1
	dY[ 41 : : 84] = dRTL
	
	dY[ 42 : : 84] = dSTH
        dY[ 43 : : 84] = dITH_H1
        dY[ 44 : : 84] = dITH_H3
	dY[ 45 : : 84] = dITH_B
        dY[ 46 : : 84] = dRTH_H1
	dY[ 47 : : 84] = dRTH_H3
	dY[ 48 : : 84] = dRTH_B
	dY[ 49 : : 84] = dITH_H1_H3
	dY[ 50 : : 84] = dITH_H1_B
	dY[ 51 : : 84] = dITH_H3_H1
	dY[ 52 : : 84] = dITH_H3_B
	dY[ 53 : : 84] = dITH_B_H3
	dY[ 54 : : 84] = dITH_B_H1
	dY[ 55 : : 84] = dRTH
	
	dY[ 56 : : 84] = dSNL
        dY[ 57 : : 84] = dINL_H1
        dY[ 58 : : 84] = dINL_H3
	dY[ 59 : : 84] = dINL_B
        dY[ 60 : : 84] = dRNL_H1
	dY[ 61 : : 84] = dRNL_H3
	dY[ 62 : : 84] = dRNL_B
	dY[ 63 : : 84] = dINL_H1_H3
	dY[ 64 : : 84] = dINL_H1_B
	dY[ 65 : : 84] = dINL_H3_H1
	dY[ 66 : : 84] = dINL_H3_B
	dY[ 67 : : 84] = dINL_B_H3
	dY[ 68 : : 84] = dINL_B_H1
	dY[ 69 : : 84] = dRNL
	
	dY[ 70 : : 84] = dSNH
        dY[ 71 : : 84] = dINH_H1
        dY[ 72 : : 84] = dINH_H3
	dY[ 73 : : 84] = dINH_B
        dY[ 74 : : 84] = dRNH_H1
	dY[ 75 : : 84] = dRNH_H3
	dY[ 76 : : 84] = dRNH_B
	dY[ 77 : : 84] = dINH_H1_H3
	dY[ 78 : : 84] = dINH_H1_B
	dY[ 79 : : 84] = dINH_H3_H1
	dY[ 80 : : 84] = dINH_H3_B
	dY[ 81 : : 84] = dINH_B_H3
	dY[ 82 : : 84] = dINH_B_H1
	dY[ 83 : : 84] = dRNH
	
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
	    IUL_H1_H3_Old = self.IUL_H1_H3.copy()
	    IUL_H1_B_Old = self.IUL_H1_B.copy()
	    IUL_H3_H1_Old = self.IUL_H3_H1.copy()
	    IUL_H3_B_Old = self.IUL_H3_B.copy()
	    IUL_B_H3_Old = self.IUL_B_H3.copy()
	    IUL_B_H1_Old = self.IUL_B_H1.copy()
	    RUL_Old = self.RUL.copy()
	    
	    SUH_Old = self.SUH.copy()
            IUH_H1_Old = self.IUH_H1.copy()
	    IUH_H3_Old = self.IUH_H3.copy()
	    IUH_B_Old = self.IUH_B.copy()
            RUH_H1_Old = self.RUH_H1.copy()
	    RUH_H3_Old = self.RUH_H3.copy()
	    RUH_B_Old = self.RUH_B.copy()
	    IUH_H1_H3_Old = self.IUH_H1_H3.copy()
	    IUH_H1_B_Old = self.IUH_H1_B.copy()
	    IUH_H3_H1_Old = self.IUH_H3_H1.copy()
	    IUH_H3_B_Old = self.IUH_H3_B.copy()
	    IUH_B_H3_Old = self.IUH_B_H3.copy()
	    IUH_B_H1_Old = self.IUH_B_H1.copy()
	    RUH_Old = self.RUH.copy()
	    
	    STL_Old = self.STL.copy()
            ITL_H1_Old = self.ITL_H1.copy()
	    ITL_H3_Old = self.ITL_H3.copy()
	    ITL_B_Old = self.ITL_B.copy()
            RTL_H1_Old = self.RTL_H1.copy()
	    RTL_H3_Old = self.RTL_H3.copy()
	    RTL_B_Old = self.RTL_B.copy()
	    ITL_H1_H3_Old = self.ITL_H1_H3.copy()
	    ITL_H1_B_Old = self.ITL_H1_B.copy()
	    ITL_H3_H1_Old = self.ITL_H3_H1.copy()
	    ITL_H3_B_Old = self.ITL_H3_B.copy()
	    ITL_B_H3_Old = self.ITL_B_H3.copy()
	    ITL_B_H1_Old = self.ITL_B_H1.copy()
	    RTL_Old = self.RTL.copy()
	    
	    STH_Old = self.STH.copy()
            ITH_H1_Old = self.ITH_H1.copy()
	    ITH_H3_Old = self.ITH_H3.copy()
	    ITH_B_Old = self.ITH_B.copy()
            RTH_H1_Old = self.RTH_H1.copy()
	    RTH_H3_Old = self.RTH_H3.copy()
	    RTH_B_Old = self.RTH_B.copy()
	    ITH_H1_H3_Old = self.ITH_H1_H3.copy()
	    ITH_H1_B_Old = self.ITH_H1_B.copy()
	    ITH_H3_H1_Old = self.ITH_H3_H1.copy()
	    ITH_H3_B_Old = self.ITH_H3_B.copy()
	    ITH_B_H3_Old = self.ITH_B_H3.copy()
	    ITH_B_H1_Old = self.ITH_B_H1.copy()
	    RTH_Old = self.RTH.copy()
	    
	    SNL_Old = self.SNL.copy()
            INL_H1_Old = self.INL_H1.copy()
	    INL_H3_Old = self.INL_H3.copy()
	    INL_B_Old = self.INL_B.copy()
            RNL_H1_Old = self.RNL_H1.copy()
	    RNL_H3_Old = self.RNL_H3.copy()
	    RNL_B_Old = self.RNL_B.copy()
	    INL_H1_H3_Old = self.INL_H1_H3.copy()
	    INL_H1_B_Old = self.INL_H1_B.copy()
	    INL_H3_H1_Old = self.INL_H3_H1.copy()
	    INL_H3_B_Old = self.INL_H3_B.copy()
	    INL_B_H3_Old = self.INL_B_H3.copy()
	    INL_B_H1_Old = self.INL_B_H1.copy()
	    RNL_Old = self.RNL.copy()
	    
	    SNH_Old = self.SNH.copy()
            INH_H1_Old = self.INH_H1.copy()
	    INH_H3_Old = self.INH_H3.copy()
	    INH_B_Old = self.INH_B.copy()
            RNH_H1_Old = self.RNH_H1.copy()
	    RNH_H3_Old = self.RNH_H3.copy()
	    RNH_B_Old = self.RNH_B.copy()
	    INH_H1_H3_Old = self.INH_H1_H3.copy()
	    INH_H1_B_Old = self.INH_H1_B.copy()
	    INH_H3_H1_Old = self.INH_H3_H1.copy()
	    INH_H3_B_Old = self.INH_H3_B.copy()
	    INH_B_H3_Old = self.INH_B_H3.copy()
	    INH_B_H1_Old = self.INH_B_H1.copy()
	    RNH_Old = self.RNH.copy()
	
        # Time vector for solution
        self.T = numpy.hstack((numpy.arange(tStart, tEnd, tStep), tEnd))
        
        # Integrate the ODE
        from scipy.integrate import odeint
        self.Y = odeint(self.RHS,
                        self.Y0.copy(),
                        self.T,
                        mxstep = 1000)
        Z = self.Y.copy()
	
	self.SUL    = Z[:, 0 : : 84]
        self.IUL_H1 = Z[:, 1 : : 84]
	self.IUL_H3 = Z[:, 2 : : 84]
	self.IUL_B  = Z[:, 3 : : 84]
	self.RUL_H1 = Z[:, 4 : : 84]
	self.RUL_H3 = Z[:, 5 : : 84]
	self.RUL_B = Z[:, 6 : : 84]
	self.IUL_H1_H3 = Z[:, 7 : : 84]
	self.IUL_H1_B = Z[:, 8 : : 84]
	self.IUL_H3_H1 = Z[:, 9 : : 84]
	self.IUL_H3_B = Z[:, 10 : : 84]
	self.IUL_B_H3 = Z[:, 11 : : 84]
	self.IUL_B_H1 = Z[:, 12 : : 84]
	self.RUL =  Z[:,13: : 84]
	
	self.SUH    = Z[:, 14 : : 84]
        self.IUH_H1 = Z[:, 15 : : 84]
	self.IUH_H3 = Z[:, 16 : : 84]
	self.IUH_B  = Z[:, 17 : : 84]
	self.RUH_H1 = Z[:, 18 : : 84]
	self.RUH_H3 = Z[:, 19 : : 84]
	self.RUH_B = Z[:, 20 : : 84]
	self.IUH_H1_H3 = Z[:, 21 : : 84]
	self.IUH_H1_B = Z[:, 22 : : 84]
	self.IUH_H3_H1 = Z[:, 23 : : 84]
	self.IUH_H3_B = Z[:, 24 : : 84]
	self.IUH_B_H3 = Z[:, 25 : : 84]
	self.IUH_B_H1 = Z[:, 26 : : 84]
	self.RUH = Z[:,27: : 84]
	
	self.STL    = Z[:, 28 : : 84]
        self.ITL_H1 = Z[:, 29 : : 84]
	self.ITL_H3 = Z[:, 30 : : 84]
	self.ITL_B  = Z[:, 31 : : 84]
	self.RTL_H1 = Z[:, 32 : : 84]
	self.RTL_H3 = Z[:, 33 : : 84]
	self.RTL_B = Z[:, 34 : : 84]
	self.ITL_H1_H3 = Z[:, 35 : : 84]
	self.ITL_H1_B = Z[:, 36 : : 84]
	self.ITL_H3_H1 = Z[:, 37 : : 84]
	self.ITL_H3_B = Z[:, 38 : : 84]
	self.ITL_B_H3 = Z[:, 39 : : 84]
	self.ITL_B_H1 = Z[:, 40 : : 84]
	self.RTL = Z[:,41: : 84]
	
	self.STH    = Z[:, 42 : : 84]
        self.ITH_H1 = Z[:, 43 : : 84]
	self.ITH_H3 = Z[:, 44 : : 84]
	self.ITH_B  = Z[:, 45 : : 84]
	self.RTH_H1 = Z[:, 46 : : 84]
	self.RTH_H3 = Z[:, 47 : : 84]
	self.RTH_B = Z[:, 48 : : 84]
	self.ITH_H1_H3 = Z[:, 49 : : 84]
	self.ITH_H1_B = Z[:, 50 : : 84]
	self.ITH_H3_H1 = Z[:, 51 : : 84]
	self.ITH_H3_B = Z[:, 52 : : 84]
	self.ITH_B_H3 = Z[:, 53 : : 84]
	self.ITH_B_H1 = Z[:, 54 : : 84]
	self.RTH = Z[:,55: : 84]
	
	self.SNL    = Z[:, 56 : : 84]
        self.INL_H1 = Z[:, 57 : : 84]
	self.INL_H3 = Z[:, 58 : : 84]
	self.INL_B  = Z[:, 59 : : 84]
	self.RNL_H1 = Z[:, 60 : : 84]
	self.RNL_H3 = Z[:, 61 : : 84]
	self.RNL_B = Z[:, 62 : : 84]
	self.INL_H1_H3 = Z[:, 63 : : 84]
	self.INL_H1_B = Z[:, 64 : : 84]
	self.INL_H3_H1 = Z[:, 65 : : 84]
	self.INL_H3_B = Z[:, 66 : : 84]
	self.INL_B_H3 = Z[:, 67 : : 84]
	self.INL_B_H1 = Z[:, 68 : : 84]
	self.RNL = Z[:,69: : 84]
	
	self.SNH    = Z[:, 70 : : 84]
        self.INH_H1 = Z[:, 71 : : 84]
	self.INH_H3 = Z[:, 72 : : 84]
	self.INH_B  = Z[:, 73 : : 84]
	self.RNH_H1 = Z[:, 74 : : 84]
	self.RNH_H3 = Z[:, 75 : : 84]
	self.RNH_B = Z[:, 76 : : 84]
	self.INH_H1_H3 = Z[:, 77 : : 84]
	self.INH_H1_B = Z[:, 78 : : 84]
	self.INH_H3_H1 = Z[:, 79 : : 84]
	self.INH_H3_B = Z[:, 80 : : 84]
	self.INH_B_H3 = Z[:, 81 : : 84]
	self.INH_B_H1 = Z[:, 82 : : 84]
	self.RNH = Z[:,83: : 84]
	
	
        if self.hasSolution:
	   
            self.T = numpy.hstack((TOld, self.T))
	    #UL
            self.SUL = numpy.vstack((SUL_Old, self.SUL))
            self.IUL_H1 = numpy.vstack((IUL_H1_Old, self.IUL_H1))
	    self.IUL_H3 = numpy.vstack((IUL_H3_Old, self.IUL_H3))
	    self.IUL_B = numpy.vstack((IUL_B_Old, self.IUL_B))
            self.RUL_H1 = numpy.vstack((RUL_H1_Old, self.RUL_H1))
	    self.RUL_H3 = numpy.vstack((RUL_H3_Old, self.RUL_H3))
	    self.RUL_B = numpy.vstack((RUL_B_Old, self.RUL_B))
	    self.IUL_H1_H3 = numpy.vstack((IUL_H1_H3_Old, self.IUL_H1_H3))
	    self.IUL_H1_B = numpy.vstack((IUL_H1_B_Old, self.IUL_H1_B))
	    self.IUL_H3_H1 = numpy.vstack((IUL_H3_H1_Old, self.IUL_H3_H1))
	    self.IUL_H3_B = numpy.vstack((IUL_H3_B_Old, self.IUL_H3_B))
	    self.IUL_B_H3 = numpy.vstack((IUL_B_H3_Old, self.IUL_B_H3))
	    self.IUL_B_H1 = numpy.vstack((IUL_B_H1_Old, self.IUL_B_H1))
	    self.RUL =  numpy.vstack((RUL_Old, self.RUL))
	    
	    #UH
	    self.SUH = numpy.vstack((SUH_Old, self.SUH))
            self.IUH_H1 = numpy.vstack((IUH_H1_Old, self.IUH_H1))
	    self.IUH_H3 = numpy.vstack((IUH_H3_Old, self.IUH_H3))
	    self.IUH_B = numpy.vstack((IUH_B_Old, self.IUH_B))
            self.RUH_H1 = numpy.vstack((RUH_H1_Old, self.RUH_H1))
	    self.RUH_H3 = numpy.vstack((RUH_H3_Old, self.RUH_H3))
	    self.RUH_B = numpy.vstack((RUH_B_Old, self.RUH_B))
	    self.IUH_H1_H3 = numpy.vstack((IUH_H1_H3_Old, self.IUH_H1_H3))
	    self.IUH_H1_B = numpy.vstack((IUH_H1_B_Old, self.IUH_H1_B))
	    self.IUH_H3_H1 = numpy.vstack((IUH_H3_H1_Old, self.IUH_H3_H1))
	    self.IUH_H3_B = numpy.vstack((IUH_H3_B_Old, self.IUH_H3_B))
	    self.IUH_B_H3 = numpy.vstack((IUH_B_H3_Old, self.IUH_B_H3))
	    self.IUH_B_H1 = numpy.vstack((IUH_B_H1_Old, self.IUH_B_H1))
	    self.RUH =  numpy.vstack((RUH_Old, self.RUH))
	    
	    #TL
	    self.STL = numpy.vstack((STL_Old, self.STL))
            self.ITL_H1 = numpy.vstack((ITL_H1_Old, self.ITL_H1))
	    self.ITL_H3 = numpy.vstack((ITL_H3_Old, self.ITL_H3))
	    self.ITL_B = numpy.vstack((ITL_B_Old, self.ITL_B))
            self.RTL_H1 = numpy.vstack((RTL_H1_Old, self.RTL_H1))
	    self.RTL_H3 = numpy.vstack((RTL_H3_Old, self.RTL_H3))
	    self.RTL_B = numpy.vstack((RTL_B_Old, self.RTL_B))
	    self.ITL_H1_H3 = numpy.vstack((ITL_H1_H3_Old, self.ITL_H1_H3))
	    self.ITL_H1_B = numpy.vstack((ITL_H1_B_Old, self.ITL_H1_B))
	    self.ITL_H3_H1 = numpy.vstack((ITL_H3_H1_Old, self.ITL_H3_H1))
	    self.ITL_H3_B = numpy.vstack((ITL_H3_B_Old, self.ITL_H3_B))
	    self.ITL_B_H3 = numpy.vstack((ITL_B_H3_Old, self.ITL_B_H3))
	    self.ITL_B_H1 = numpy.vstack((ITL_B_H1_Old, self.ITL_B_H1))
	    self.RTL =  numpy.vstack((RTL_Old, self.RTL))
	    
	    #TH
	    self.STH = numpy.vstack((STH_Old, self.STH))
            self.ITH_H1 = numpy.vstack((ITH_H1_Old, self.ITH_H1))
	    self.ITH_H3 = numpy.vstack((ITH_H3_Old, self.ITH_H3))
	    self.ITH_B = numpy.vstack((ITH_B_Old, self.ITH_B))
            self.RTH_H1 = numpy.vstack((RTH_H1_Old, self.RTH_H1))
	    self.RTH_H3 = numpy.vstack((RTH_H3_Old, self.RTH_H3))
	    self.RTH_B = numpy.vstack((RTH_B_Old, self.RTH_B))
	    self.ITH_H1_H3 = numpy.vstack((ITH_H1_H3_Old, self.ITH_H1_H3))
	    self.ITH_H1_B = numpy.vstack((ITH_H1_B_Old, self.ITH_H1_B))
	    self.ITH_H3_H1 = numpy.vstack((ITH_H3_H1_Old, self.ITH_H3_H1))
	    self.ITH_H3_B = numpy.vstack((ITH_H3_B_Old, self.ITH_H3_B))
	    self.ITH_B_H3 = numpy.vstack((ITH_B_H3_Old, self.ITH_B_H3))
	    self.ITH_B_H1 = numpy.vstack((ITH_B_H1_Old, self.ITH_B_H1))
	    self.RTH =  numpy.vstack((RTH_Old, self.RTH))
	    
	    #NL
	    self.SNL = numpy.vstack((SNL_Old, self.SNL))
            self.INL_H1 = numpy.vstack((INL_H1_Old, self.INL_H1))
	    self.INL_H3 = numpy.vstack((INL_H3_Old, self.INL_H3))
	    self.INL_B = numpy.vstack((INL_B_Old, self.INL_B))
            self.RNL_H1 = numpy.vstack((RNL_H1_Old, self.RNL_H1))
	    self.RNL_H3 = numpy.vstack((RNL_H3_Old, self.RNL_H3))
	    self.RNL_B = numpy.vstack((RNL_B_Old, self.RNL_B))
	    self.INL_H1_H3 = numpy.vstack((INL_H1_H3_Old, self.INL_H1_H3))
	    self.INL_H1_B = numpy.vstack((INL_H1_B_Old, self.INL_H1_B))
	    self.INL_H3_H1 = numpy.vstack((INL_H3_H1_Old, self.INL_H3_H1))
	    self.INL_H3_B = numpy.vstack((INL_H3_B_Old, self.INL_H3_B))
	    self.INL_B_H3 = numpy.vstack((INL_B_H3_Old, self.INL_B_H3))
	    self.INL_B_H1 = numpy.vstack((INL_B_H1_Old, self.INL_B_H1))
	    self.RNL =  numpy.vstack((RNL_Old, self.RNL))
	    
	    #NH
	    self.SNH = numpy.vstack((SNH_Old, self.SNH))
            self.INH_H1 = numpy.vstack((INH_H1_Old, self.INH_H1))
	    self.INH_H3 = numpy.vstack((INH_H3_Old, self.INH_H3))
	    self.INH_B = numpy.vstack((INH_B_Old, self.INH_B))
            self.RNH_H1 = numpy.vstack((RNH_H1_Old, self.RNH_H1))
	    self.RNH_H3 = numpy.vstack((RNH_H3_Old, self.RNH_H3))
	    self.RNH_B = numpy.vstack((RNH_B_Old, self.RNH_B))
	    self.INH_H1_H3 = numpy.vstack((INH_H1_H3_Old, self.INH_H1_H3))
	    self.INH_H1_B = numpy.vstack((INH_H1_B_Old, self.INH_H1_B))
	    self.INH_H3_H1 = numpy.vstack((INH_H3_H1_Old, self.INH_H3_H1))
	    self.INH_H3_B = numpy.vstack((INH_H3_B_Old, self.INH_H3_B))
	    self.INH_B_H3 = numpy.vstack((INH_B_H3_Old, self.INH_B_H3))
	    self.INH_B_H1 = numpy.vstack((INH_B_H1_Old, self.INH_B_H1))
	    self.RNH =  numpy.vstack((RNH_Old, self.RNH))
	    
	    
        self.hasSolution = True

    def updateStats(self):
	    
	self.NUL = self.SUL +  self.IUL_H1 + self.IUL_H3 + self.IUL_B + self.RUL_H1 + self.RUL_H3 + self.RUL_B + self.IUL_H1_H3 + self.IUL_H1_B+ self.IUL_H3_H1+ self.IUL_H3_B+ self.IUL_B_H3+ self.IUL_B_H1+ self.RUL
	self.NUH = self.SUH +  self.IUH_H1 + self.IUH_H3 + self.IUH_B + self.RUH_H1 + self.RUH_H3 + self.RUH_B + self.IUH_H1_H3 + self.IUH_H1_B+ self.IUH_H3_H1+ self.IUH_H3_B+ self.IUH_B_H3+ self.IUH_B_H1+ self.RUH
	self.NTL = self.STL +  self.ITL_H1 + self.ITL_H3 + self.ITL_B + self.RTL_H1 + self.RTL_H3 + self.RTL_B + self.ITL_H1_H3 + self.ITL_H1_B+ self.ITL_H3_H1+ self.ITL_H3_B+ self.ITL_B_H3+ self.ITL_B_H1+ self.RTL
	self.NTH = self.STH +  self.ITH_H1 + self.ITH_H3 + self.ITH_B + self.RTH_H1 + self.RTH_H3 + self.RTH_B + self.ITH_H1_H3 + self.ITH_H1_B+ self.ITH_H3_H1+ self.ITH_H3_B+ self.ITH_B_H3+ self.ITH_B_H1+ self.RTH
	self.NNL = self.SNL +  self.INL_H1 + self.INL_H3 + self.INL_B + self.RNL_H1 + self.RNL_H3 + self.RNL_B + self.INL_H1_H3 + self.INL_H1_B+ self.INL_H3_H1+ self.INL_H3_B+ self.INL_B_H3+ self.INL_B_H1+ self.RNL
	self.NNH = self.SNH +  self.INH_H1 + self.INH_H3 + self.INH_B + self.RNH_H1 + self.RNH_H3 + self.RNH_B + self.INH_H1_H3 + self.INH_H1_B+ self.INH_H3_H1+ self.INH_H3_B+ self.INH_B_H3+ self.INH_B_H1+ self.RNH
	
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
	
	self.infectionsL_H1 = 

        # Find duplicate times: these are where vaccination occurs
        #for i in numpy.arange(len(self.T)).compress(numpy.diff(self.T) == 0):
            # Update for vaccinated people
        #    self.infectionsUL += self.SUL[i + 1, :] - self.SUL[i, :]
	#    self.infectionsUH += self.SUH[i + 1, :] - self.SUH[i, :]
        #    self.infectionsTL += self.STL[i + 1, :] - self.STL[i, :]
	#    self.infectionsTH += self.STH[i + 1, :] - self.STH[i, :]
	#    self.infectionsNL += self.SNL[i + 1, :] - self.SNL[i, :]
	#    self.infectionsNH += self.SNH[i + 1, :] - self.SNH[i, :]

	self.infectionsL  = self.infectionsUL + self.infectionsTL + self.infectionsNL
        self.infectionsH  = self.infectionsUH + self.infectionsTH + self.infectionsNH
	self.infectionsU =  self.infectionsUL + self.infectionsUH
	self.infectionsV  = self.infectionsTL + self.infectionsTH + self.infectionsNL + self.infectionsNH        
	self.infections  = self.infectionsU + self.infectionsV
        self.totalInfections = self.infections.sum()
	
	"""
	self.hospitalizationsL_H1 = self.infectionsL_H1 * self.parameters.lowRiskcaseHospitalization_H1
	self.hospitalizationsL_H3 = self.infectionsL_H3 * self.parameters.lowRiskcaseHospitalization_H3
	self.hospitalizationsL_B = self.infectionsL_B * self.parameters.lowRiskcaseHospitalization_B
	
	self.hospitalizationsH_H1 = self.infectionsH_H1 * self.parameters.highRiskcaseHospitalization_H1
	self.hospitalizationsH_H3 = self.infectionsH_H3 * self.parameters.highRiskcaseHospitalization_H3
	self.hospitalizationsH_B = self.infectionsH_B * self.parameters.highRiskcaseHospitalization_B
	
	self.hospitalizationsL  = self.hospitalizationsL_H1 + self.hospitalizationsL_H3 + self.hospitalizationsL_B
	self.hospitalizationsH  = self.hospitalizationsH_H1 + self.hospitalizationsH_H3 + self.hospitalizationsH_B
        self.hospitalizations = self.hospitalizationsL + self.hospitalizationsH
	self.totalHospitalizations = self.hospitalizations.sum()
        
        self.deathsUL = self.NUL[0, :] - self.NUL[-1, :]
	self.deathsUH = self.NUH[0, :] - self.NUH[-1, :]
	self.deathsTL = self.NTL[0, :] - self.NTL[-1, :]
	self.deathsTH = self.NTH[0, :] - self.NTH[-1, :]
	self.deathsNL = self.NNL[0, :] - self.NNL[-1, :]
	self.deathsNH = self.NNH[0, :] - self.NNH[-1, :]
	
        # Find duplicate times: these are where vaccination occurs
        #for i in numpy.arange(len(self.T)).compress(numpy.diff(self.T) == 0):
            # Update for vaccinated people
        #    self.deathsUL += self.NUL[i + 1, :] - self.NUL[i, :]
	#    self.deathsUH += self.NUH[i + 1, :] - self.NUH[i, :]
        #    self.deathsTL += self.NTL[i + 1, :] - self.NTL[i, :]
	#    self.deathsTH += self.NTH[i + 1, :] - self.NTH[i, :]
	#    self.deathsNL += self.NNL[i + 1, :] - self.NNL[i, :]
	#    self.deathsNH += self.NNH[i + 1, :] - self.NNH[i, :]


	self.deathsL = self.deathsUL + self.deathsTL + self.deathsNL
	self.deathsH = self.deathsUH + self.deathsTH + self.deathsNH
	self.deathsU = self.deathsUL + self.deathsUH
	self.deathsV = self.deathsTL + self.deathsTH + self.deathsNL + self.deathsNH
        self.deaths  = self.deathsL + self.deathsH
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
	
	"""
        
        
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


	#remove age group zero from all sections here
	population_highrisk = [(a*b) for (a,b) in zip(self.parameters.population[1:] , self.parameters.proportionHighRisk[1:])]
        population_lowrisk = self.parameters.population[1:] - population_highrisk
	
	## raw coverage among age groups
	proportionVaccinatedTL = PVPWVal[0][: self.parameters.proportionVaccinatedLLength]
	proportionVaccinatedTH = list(PVPWVal[0][self.parameters.proportionVaccinatedLLength:])*len(proportionVaccinatedTL)
	proportionVaccinatedNL = PVPWVal[1][: self.parameters.proportionVaccinatedLLength]
	proportionVaccinatedNH = list(PVPWVal[1][self.parameters.proportionVaccinatedLLength:])*len(proportionVaccinatedNL)
	
	
	#raw dose uptake among age groups
	dosesVaccinatedTL =  [(a*b) for a, b in zip(proportionVaccinatedTL, population_lowrisk)]
	dosesVaccinatedTH =  [(a*b) for a, b in zip(proportionVaccinatedTH, population_highrisk)]
	dosesVaccinatedNL =  [(a*b) for a, b in zip(proportionVaccinatedNL, population_lowrisk)]
	dosesVaccinatedNH =  [(a*b) for a, b in zip(proportionVaccinatedNH, population_highrisk)]
	
	
	
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
	##values with the matrix defined in S.142
	
	
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
	
	unvax  = ((self.SUL[0,:].sum() + self.SUH[0,:].sum() - self.SUL[-1,:].sum() - self.SUH[-1,:].sum()))/1e6
	typical = (self.STL[0,:].sum() + self.STH[0,:].sum() - self.STL[-1,:].sum() - self.STH[-1,:].sum())/1e6
	universal = ((self.SNL[0,:].sum() + self.SNH[0,:].sum() - self.SNL[-1,:].sum() - self.SNH[-1,:].sum()))/1e6
	print ("Unvax infected"), unvax
	print ("Typical infected"), typical
	print ("Universal infected"), universal
	print ("total"), unvax + typical + universal

	import matplotlib.pyplot as plt
	times = [num for num in xrange(tEnd+1)]
	plt.plot(times, (self.IUL_H3+ self.IUL_H1+ self.IUL_B).sum(axis=1), color = "red", linestyle= "-")
	plt.plot(times, (self.IUH_H3+ self.IUH_H1+ self.IUH_B).sum(axis=1), color = "red", linestyle= "--")
	plt.plot(times, (self.ITL_H3+ self.ITL_H1+ self.ITL_B).sum(axis=1), color = "blue", linestyle = "-")
	plt.plot(times, (self.ITH_H3+ self.ITH_H1+ self.ITH_B).sum(axis=1), color = "blue", linestyle = "--")
	plt.plot(times, (self.INL_H3+ self.INL_H1+ self.INL_B).sum(axis=1), color = "green", linestyle = "-")
	plt.plot(times, (self.INH_H3+ self.INH_H1+ self.INH_B).sum(axis=1), color = "green", linestyle = "--")

	
		 
	plt.show()


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


