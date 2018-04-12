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
        self.Y0 = numpy.zeros(54 * self.parameters.ages.size)

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
	
	    self.vacNumbers = doses
	    
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
	    typical_vaccination = [ 0.70,  0.56,  0.56, 0.5152, 0.34, 0.34,0.34,0.34,0.34,0.34,
				   0.46,  0.46, 0.46, 0.653, 0.653, 0.653, 0.7]
	    
	    typical_proportion = self.vacNumbers[0]/(1. *sum(self.vacNumbers))
	    universal_proportion = 1 - typical_proportion
	    
	    typical_vax_coverage = [typical_proportion * num for num in typical_vaccination]
	    universal_vax_coverage = [(a-b) for (a,b) in zip(typical_vaccination, typical_vax_coverage)]
	    
	    #print ("check vaccine numbers"), sum([(a+b)*c for (a,b,c) in zip(typical_vax_coverage, universal_vax_coverage, self.parameters.population[1:])])		   

	    
	    return typical_vax_coverage + universal_vax_coverage	
	    

    def computeR0(self):
        return self.parameters.computeR0()

    def getLastValues(self):
        return (self.SUL_H3[-1, :], self.IUL_H1[-1, :], self.RUL_H1[-1, :],
		self.SUL_H3[-1, :], self.IUL_H3[-1, :], self.RUL_H3[-1, :],
		self.SUL_B[-1, :], self.IUL_HB[-1, :], self.RUL_HB[-1, :],
		
		self.SUH_H1[-1, :], self.IUH_H1[-1, :], self.RUH_H1[-1, :],
		self.SUH_H3[-1, :], self.IUH_H3[-1, :], self.RUH_H3[-1, :],
		self.SUH_B[-1, :], self.IUH_HB[-1, :], self.RUH_HB[-1, :],
		
		self.STL_H1[-1, :], self.ITL_H1[-1, :], self.RTL_H1[-1, :],
		self.STL_H3[-1, :], self.ITL_H3[-1, :], self.RTL_H3[-1, :],
		self.STL_B[-1, :], self.ITL_HB[-1, :], self.RTL_HB[-1, :],
		
		self.STH_H1[-1, :], self.ITH_H1[-1, :], self.RTH_H1[-1, :],
		self.STH_H3[-1, :], self.ITH_H3[-1, :], self.RTH_H3[-1, :],
		self.STH_B[-1, :], self.ITH_HB[-1, :], self.RTH_HB[-1, :],
		
		self.SNL_H1[-1, :], self.INL_H1[-1, :], self.RNL_H1[-1, :],
		self.SNL_H3[-1, :], self.INL_H3[-1, :], self.RNL_H3[-1, :],
		self.SNL_B[-1, :], self.INL_HB[-1, :], self.RNL_HB[-1, :],
		
		self.SNH_H1[-1, :], self.INH_H1[-1, :], self.RNH_H1[-1, :],
		self.SNH_H3[-1, :], self.INH_H3[-1, :], self.RNH_H3[-1, :],
		self.SNH_B[-1, :], self.INH_HB[-1, :], self.RNH_HB[-1, :])
    
    def updateIC(self):
        if not self.hasSolution:
            # S
	    ## SUL
            self.Y0[ 0 : : 54] = \
                     (1 - self.parameters.proportionVaccinatedL - self.parameters.proportionVaccinatedH) * self.parameters.population * (1 - self.parameters.proportionHighRisk)
	    self.Y0[ 3 : : 54] = \
                     (1 - self.parameters.proportionVaccinatedL - self.parameters.proportionVaccinatedH) * self.parameters.population * (1 - self.parameters.proportionHighRisk)
	    self.Y0[ 6 : : 54] = \
                     (1 - self.parameters.proportionVaccinatedL - self.parameters.proportionVaccinatedH) * self.parameters.population * (1 - self.parameters.proportionHighRisk)
	    
	    
	    ## SUH
            self.Y0[ 9 : : 54] = \
                     (1 - self.parameters.proportionVaccinatedL - self.parameters.proportionVaccinatedH) * self.parameters.population * self.parameters.proportionHighRisk 
	    self.Y0[ 12 : : 54] = \
                     (1 - self.parameters.proportionVaccinatedL - self.parameters.proportionVaccinatedH) * self.parameters.population * self.parameters.proportionHighRisk 
	    self.Y0[ 15 : : 54] = \
                     (1 - self.parameters.proportionVaccinatedL - self.parameters.proportionVaccinatedH) * self.parameters.population * self.parameters.proportionHighRisk 
	    
	    ## STL
            self.Y0[ 18 : : 54] = self.parameters.proportionVaccinatedTL * self.parameters.population * (1 - self.parameters.proportionHighRisk)
	    self.Y0[ 21 : : 54] = self.parameters.proportionVaccinatedTL * self.parameters.population * (1 - self.parameters.proportionHighRisk)
	    self.Y0[ 24 : : 54] = self.parameters.proportionVaccinatedTL * self.parameters.population * (1 - self.parameters.proportionHighRisk)
	    
	    
	    ## STH
            self.Y0[ 27 : : 54] =  self.parameters.proportionVaccinatedTH * self.parameters.population * self.parameters.proportionHighRisk
	    self.Y0[ 30 : : 54] =  self.parameters.proportionVaccinatedTH * self.parameters.population * self.parameters.proportionHighRisk
	    self.Y0[ 33 : : 54] =  self.parameters.proportionVaccinatedTH * self.parameters.population * self.parameters.proportionHighRisk
	    
	    ## SNL
            self.Y0[ 36 : : 54] = self.parameters.proportionVaccinatedNL * self.parameters.population  * (1 - self.parameters.proportionHighRisk)
	    self.Y0[ 39 : : 54] = self.parameters.proportionVaccinatedNL * self.parameters.population  * (1 - self.parameters.proportionHighRisk)
	    self.Y0[ 42 : : 54] = self.parameters.proportionVaccinatedNL * self.parameters.population  * (1 - self.parameters.proportionHighRisk)
	    
	    
	    ## SNH 
            self.Y0[ 45 : : 54] = self.parameters.proportionVaccinatedNH * self.parameters.population * self.parameters.proportionHighRisk
	    self.Y0[ 48 : : 54] = self.parameters.proportionVaccinatedNH * self.parameters.population * self.parameters.proportionHighRisk
	    self.Y0[ 51 : : 54] = self.parameters.proportionVaccinatedNH * self.parameters.population * self.parameters.proportionHighRisk


            # I: Add a single infectious person in each age
	    
	    self.Y0[ 1 : : 54] = numpy.full(self.parameters.ages.size,1)
	    self.Y0[ 4 : : 54] = numpy.full(self.parameters.ages.size,1)
	    self.Y0[ 7 : : 54] = numpy.full(self.parameters.ages.size,1)
	    
            self.Y0[ 10 : : 54] = numpy.full(self.parameters.ages.size, 1)
	    self.Y0[ 13 : : 54] = numpy.full(self.parameters.ages.size, 1)
	    self.Y0[ 16 : : 54] = numpy.full(self.parameters.ages.size, 1)
	    
	    self.Y0[ 19 : : 54] = numpy.full(self.parameters.ages.size, 1)
	    self.Y0[ 22 : : 54] = numpy.full(self.parameters.ages.size, 1)
	    self.Y0[ 25 : : 54] = numpy.full(self.parameters.ages.size, 1)
	    
	    self.Y0[28 : : 54] = numpy.full(self.parameters.ages.size, 1)
	    self.Y0[31 : : 54] = numpy.full(self.parameters.ages.size, 1)
	    self.Y0[34 : : 54] = numpy.full(self.parameters.ages.size, 1)
	    
	    self.Y0[37 : : 54] = numpy.full(self.parameters.ages.size, 1)
	    self.Y0[40 : : 54] = numpy.full(self.parameters.ages.size, 1)
	    self.Y0[43 : : 54] = numpy.full(self.parameters.ages.size, 1)
	    
	    self.Y0[46 : : 54] = numpy.full(self.parameters.ages.size, 1)
	    self.Y0[49 : : 54] = numpy.full(self.parameters.ages.size, 1)
	    self.Y0[52 : : 54] = numpy.full(self.parameters.ages.size, 1)
	    

            # S: Remove those new infectious people from the susceptibles
            self.Y0[ 0 : : 54] -= self.Y0[ 1 : : 54]
            self.Y0[ 3 : : 54] -= self.Y0[ 4 : : 54]
	    self.Y0[ 6 : : 54] -= self.Y0[ 7 : : 54]
	    self.Y0[ 9 : : 54] -= self.Y0[ 10 : : 54]
	    self.Y0[ 12 : : 54] -= self.Y0[13 : : 54]
	    self.Y0[ 15 : : 54] -= self.Y0[16 : : 54]
	    self.Y0[ 18 : : 54] -= self.Y0[19 : : 54]
	    self.Y0[ 21 : : 54] -= self.Y0[22 : : 54]
	    self.Y0[ 24 : : 54] -= self.Y0[25 : : 54]
	    self.Y0[ 27 : : 54] -= self.Y0[28 : : 54]
	    self.Y0[ 30 : : 54] -= self.Y0[31 : : 54]
	    self.Y0[ 33 : : 54] -= self.Y0[34 : : 54]
	    self.Y0[ 36 : : 54] -= self.Y0[37 : : 54]
	    self.Y0[ 39 : : 54] -= self.Y0[40 : : 54]
	    self.Y0[ 42 : : 54] -= self.Y0[43 : : 54]
	    self.Y0[ 45 : : 54] -= self.Y0[46 : : 54]
	    self.Y0[ 48 : : 54] -= self.Y0[49 : : 54]
	    self.Y0[ 51 : : 54] -= self.Y0[52 : : 54]
	    
            # R (RU, RVL, RVH)
            self.Y0[ 2 : : 54] = 0.
            self.Y0[ 5 : : 54] = 0.
	    self.Y0[ 8 : : 54] = 0.
	    self.Y0[ 11 : : 54] = 0.
	    self.Y0[ 14 : : 54] = 0.
	    self.Y0[ 17 : : 54] = 0.
	    self.Y0[ 20 : : 54] = 0.
	    self.Y0[ 23 : : 54] = 0.
	    self.Y0[ 26 : : 54] = 0.
	    self.Y0[ 29 : : 54] = 0.
	    self.Y0[ 32 : : 54] = 0.
	    self.Y0[ 35 : : 54] = 0.
	    self.Y0[ 38 : : 54] = 0.
	    self.Y0[ 41 : : 54] = 0.
	    self.Y0[ 44 : : 54] = 0.
	    self.Y0[ 47 : : 54] = 0.
	    self.Y0[ 50 : : 54] = 0.
	    self.Y0[ 53 : : 54] = 0.
	   
        else:
	    
            SUL_H1, IUL_H1, RUL_H1, SUL_H3, IUL_H3, RUL_H3, SUL_B, IUL_B, RUL_B,
	    SUH_H1, IUH_H1, RUH_H1, SUH_H3, IUH_H3, RUH_H3, SUH_B, IUH_B, RUH_B,
	    STL_H1, ITL_H1, RTL_H1, STL_H3, ITL_H3, RTL_H3, STL_B, ITL_B, RTL_B,
	    STH_H1, ITH_H1, RTH_H1, STH_H3, ITH_H3, RTH_H3, STH_B, ITH_B, RTH_B,
	    SNL_H1, INL_H1, RNL_H1, SNL_H3, INL_H3, RNL_H3, SNL_B, INL_B, RNL_B,
	    SNH_H1, INH_H1, RNH_H1, SNH_H3, INH_H3, RNH_H3, SNH_B, INH_B, RNH_B  = self.getLastValues()
            
            self.Y0[ 0 : : 54] = (1 - self.parameters.proportionVaccinatedL) * SUL_H1
	    self.Y0[ 3 : : 54] = (1 - self.parameters.proportionVaccinatedL) * SUL_H3
	    self.Y0[ 6 : : 54] = (1 - self.parameters.proportionVaccinatedL) * SUL_B
	    
	    self.Y0[ 9 : : 54] = (1 - self.parameters.proportionVaccinatedH) * SUH_H1
	    self.Y0[ 12 : : 54] = (1 - self.parameters.proportionVaccinatedH) * SUH_H3
	    self.Y0[ 15 : : 54] = (1 - self.parameters.proportionVaccinatedH) * SUH_B
	    
	    self.Y0[ 18 : : 54] = (1 - self.parameters.proportionVaccinatedH) * STL_H1
	    self.Y0[ 21 : : 54] = (1 - self.parameters.proportionVaccinatedH) * STL_H3
	    self.Y0[ 24 : : 54] = (1 - self.parameters.proportionVaccinatedH) * STL_B
	    
	    self.Y0[ 27 : : 54] = (1 - self.parameters.proportionVaccinatedH) * STH_H1
	    self.Y0[ 30 : : 54] = (1 - self.parameters.proportionVaccinatedH) * STH_H3
	    self.Y0[ 33 : : 54] = (1 - self.parameters.proportionVaccinatedH) * STH_B
	    
	    self.Y0[ 36 : : 54] = (1 - self.parameters.proportionVaccinatedH) * SNL_H1
	    self.Y0[ 39 : : 54] = (1 - self.parameters.proportionVaccinatedH) * SNL_H3
	    self.Y0[ 42 : : 54] = (1 - self.parameters.proportionVaccinatedH) * SNL_B
	    
	    self.Y0[ 45 : : 54] = (1 - self.parameters.proportionVaccinatedH) * SNH_H1
	    self.Y0[ 48 : : 54] = (1 - self.parameters.proportionVaccinatedH) * SNH_H3
	    self.Y0[ 51 : : 54] = (1 - self.parameters.proportionVaccinatedH) * SNH_B
	    
	    #I
	    self.Y0[ 1 : : 54] = IUL_H1
	    self.Y0[ 4 : : 54] = IUL_H3
	    self.Y0[ 7 : : 54] = IUL_B
	    
	    self.Y0[ 10 : : 54] = IUH_H1
	    self.Y0[ 13 : : 54] = IUH_H3
	    self.Y0[ 16 : : 54] = IUH_B
	    
	    self.Y0[ 19 : : 54] = ITL_H1
	    self.Y0[ 22 : : 54] = ITL_H3
	    self.Y0[ 25 : : 54] = ITL_B
	    
	    self.Y0[ 28 : : 54] = ITH_H1
	    self.Y0[ 31 : : 54] = ITH_H3
	    self.Y0[ 34 : : 54] = ITH_B
	    
	    self.Y0[ 37 : : 54] = INL_H1
	    self.Y0[ 40 : : 54] = INL_H3
	    self.Y0[ 43 : : 54] = INL_B
	    
	    self.Y0[ 46 : : 54] = INH_H1
	    self.Y0[ 48 : : 54] = INH_H3
	    self.Y0[ 52 : : 54] = INH_B
	    
	    #R
	    self.Y0[ 2 : : 54] = RUL_H1
	    self.Y0[ 5 : : 54] = RUL_H3
	    self.Y0[ 8 : : 54] = RUL_B
	    
	    self.Y0[ 11 : : 54] = RUH_H1
	    self.Y0[ 14 : : 54] = RUH_H3
	    self.Y0[ 17 : : 54] = RUH_B
	    
	    self.Y0[ 20 : : 54] = RTL_H1
	    self.Y0[ 23 : : 54] = RTL_H3
	    self.Y0[ 26 : : 54] = RTL_B
	    
	    self.Y0[ 29 : : 54] = RTH_H1
	    self.Y0[ 32 : : 54] = RTH_H3
	    self.Y0[ 35 : : 54] = RTH_B
	    
	    self.Y0[ 38 : : 54] = RNL_H1
	    self.Y0[ 41 : : 54] = RNL_H3
	    self.Y0[ 44 : : 54] = RNL_B
	    
	    self.Y0[ 47 : : 54] = RNH_H1
	    self.Y0[ 50 : : 54] = RNH_H3
	    self.Y0[ 53 : : 54] = RNH_B
	
	    
       

    def RHS(self, Y, t):
        '''
        SEIR model with multiple host types.
        
        This function gives the right-hand sides of the ODEs.
        '''
        
        # Convert vector to meaningful component vectors

	SUL_H1 = Y[ 0 : : 54]
        IUL_H1 = Y[ 1 : : 54]
        RUL_H1 = Y[ 2 : : 54]
	
	SUL_H3 = Y[ 0 : : 54]
        IUL_H3 = Y[ 1 : : 54]
        RUL_H3 = Y[ 2 : : 54]
	
	SUL_B = Y[ 0 : : 54]
        IUL_B = Y[ 1 : : 54]
        RUL_B = Y[ 2 : : 54]
	
	#UH
	SUH_H1 = Y[ 0 : : 54]
        IUH_H1 = Y[ 1 : : 54]
        RUH_H1 = Y[ 2 : : 54]
	
	SUH_H3 = Y[ 0 : : 54]
        IUH_H3 = Y[ 1 : : 54]
        RUH_H3 = Y[ 2 : : 54]
	
	SUH_B = Y[ 0 : : 54]
        IUH_B = Y[ 1 : : 54]
        RUH_B = Y[ 2 : : 54]
	
	#TL
	STL_H1 = Y[ 0 : : 54]
        ITL_H1 = Y[ 1 : : 54]
        RTL_H1 = Y[ 2 : : 54]
	
	STL_H3 = Y[ 0 : : 54]
        ITL_H3 = Y[ 1 : : 54]
        RTL_H3 = Y[ 2 : : 54]
	
	STL_B = Y[ 0 : : 54]
        ITL_B = Y[ 1 : : 54]
        RTL_B = Y[ 2 : : 54]
	
	#TH
	STH_H1 = Y[ 0 : : 54]
        ITH_H1 = Y[ 1 : : 54]
        RTH_H1 = Y[ 2 : : 54]
	
	STH_H3 = Y[ 0 : : 54]
        ITH_H3 = Y[ 1 : : 54]
        RTH_H3 = Y[ 2 : : 54]
	
	STH_B = Y[ 0 : : 54]
        ITH_B = Y[ 1 : : 54]
        RTH_B = Y[ 2 : : 54]
	
	#NL
	SNL_H1 = Y[ 0 : : 54]
        INL_H1 = Y[ 1 : : 54]
        RNL_H1 = Y[ 2 : : 54]
	
	SNL_H3 = Y[ 0 : : 54]
        INL_H3 = Y[ 1 : : 54]
        RNL_H3 = Y[ 2 : : 54]
	
	SNL_B = Y[ 0 : : 54]
        INL_B = Y[ 1 : : 54]
        RNL_B = Y[ 2 : : 54]
	
	#NH
	SNH_H1 = Y[ 0 : : 54]
        INH_H1 = Y[ 1 : : 54]
        RNH_H1 = Y[ 2 : : 54]
	
	SNH_H3 = Y[ 0 : : 54]
        INH_H3 = Y[ 1 : : 54]
        RNH_H3 = Y[ 2 : : 54]
	
	SNH_B = Y[ 0 : : 54]
        INH_B = Y[ 1 : : 54]
        RNH_B = Y[ 2 : : 54]
	
        
        N = sum(SUL_H1+ IUL_H1+ RUL_H1+ SUL_H3+ IUL_H3+ RUL_H3+ SUL_B+ IUL_B+ RUL_B+
	    SUH_H1+ IUH_H1+ RUH_H1+ SUH_H3+ IUH_H3+ RUH_H3+ SUH_B+ IUH_B+ RUH_B+
	    STL_H1+ ITL_H1+ RTL_H1+ STL_H3+ ITL_H3+ RTL_H3+ STL_B+ ITL_B+ RTL_B+
	    STH_H1+ ITH_H1+ RTH_H1+ STH_H3+ ITH_H3+ RTH_H3+ STH_B+ ITH_B+ RTH_B+
	    SNL_H1+ INL_H1+ RNL_H1+ SNL_H3+ INL_H3+ RNL_H3+ SNL_B+ INL_B+ RNL_B+
	    SNH_H1+ INH_H1+ RNH_H1+ SNH_H3+ INH_H3+ RNH_H3+ SNH_B+ INH_B+ RNH_B)
	    
      
        # The force of infection
        Lambda_H1 = self.parameters.transmissionScaling * self.parameters.susceptibility_H1\
		    * numpy.dot(self.parameters.contactMatrix, self.parameters.transmissibility * (IUL_H1 + IUH_H1 + ITL_H1 + ITH_H1 + INL_H1 + INH_H1)) / N
	
	Lambda_H3 = self.parameters.transmissionScaling * self.parameters.susceptibility_H3 \
                 * numpy.dot(self.parameters.contactMatrix, 
                             self.parameters.transmissibility * (IUL_H3 + IUH_H3 + ITL_H3 + ITH_H3 + INL_H3 + INH_H3)) / N
		
	Lambda_B = self.parameters.transmissionScaling * self.parameters.susceptibility_B \
                 * numpy.dot(self.parameters.contactMatrix,
                             self.parameters.transmissibility * (IUL_B + IUH_B + ITL_B + ITH_B + INL_B + INH_B)) / N
        
        # The right-hand sides
	
	#UL
        dSUL_H1 = - Lambda_H1 * SUL_H1
        dIUL_H1 = Lambda_H1 * SUL_H1 - (self.parameters.recoveryRate + self.parameters.deathRateUL_H1) * IUL_H1
        dRUL_H1 = self.parameters.recoveryRate * IUL_H1
	
	dSUL_H3 = - Lambda_H3 * SUL_H3
        dIUL_H3 = Lambda_H3 * SUL_H3 - (self.parameters.recoveryRate + self.parameters.deathRateUL_H3) * IUL_H3
        dRUL_H3 = self.parameters.recoveryRate * IUL_H3
	
	dSUL_B = - Lambda_B * SUL_B
        dIUL_B = Lambda_B * SUL_B - (self.parameters.recoveryRate + self.parameters.deathRateUL_B) * IUL_B
        dRUL_B = self.parameters.recoveryRate * IUL_B
	
	#UH
	dSUH_H1 = - Lambda_H1 * SUH_H1
        dIUH_H1 = Lambda_H1 * SUH_H1 - (self.parameters.recoveryRate + self.parameters.deathRateUH_H1) * IUH_H1
        dRUH_H1 = self.parameters.recoveryRate * IUH_H1
	
	dSUH_H3 = - Lambda_H3 * SUH_H3
        dIUH_H3 = Lambda_H3 * SUH_H3 - (self.parameters.recoveryRate + self.parameters.deathRateUH_H3) * IUH_H3
        dRUH_H3 = self.parameters.recoveryRate * IUH_H3
	
	dSUH_B = - Lambda_B * SUH_B
        dIUH_B = Lambda_B * SUH_B - (self.parameters.recoveryRate + self.parameters.deathRateUH_B) * IUH_B
        dRUH_B = self.parameters.recoveryRate * IUH_B
	
	#TL
	dSTL_H1 = - (1 - self.parameters.vaccineEfficacyVsInfectionTypical_H1) * Lambda_H1 * STL_H1
        dITL_H1 = Lambda_H1 * STL_H1 - (self.parameters.recoveryRate + self.parameters.deathRateVL_H1) * ITL_H1
        dRTL_H1 = self.parameters.recoveryRate * ITL_H1
	
	dSTL_H3 = - (1 - self.parameters.vaccineEfficacyVsInfectionTypical_H3) * Lambda_H3 * STL_H3
        dITL_H3 = Lambda_H3 * STL_H3 - (self.parameters.recoveryRate + self.parameters.deathRateVL_H3) * ITL_H3
        dRTL_H3 = self.parameters.recoveryRate * ITL_H3
	
	dSTL_B = - (1 - self.parameters.vaccineEfficacyVsInfectionTypical_B) * Lambda_B * STL_B
        dITL_B = Lambda_B * STL_B - (self.parameters.recoveryRate + self.parameters.deathRateVL_B) * ITL_B
        dRTL_B = self.parameters.recoveryRate * ITL_B
	
	#TH
	dSTH_H1 = - (1 - self.parameters.vaccineEfficacyVsInfectionTypical_H1) * Lambda_H1 * STH_H1
        dITH_H1 = Lambda_H1 * STH_H1 - (self.parameters.recoveryRate + self.parameters.deathRateVH_H1) * ITH_H1
        dRTH_H1 = self.parameters.recoveryRate * ITH_H1
	
	dSTH_H3 = - (1 - self.parameters.vaccineEfficacyVsInfectionTypical_H3) * Lambda_H3 * STH_H3
        dITH_H3 = Lambda_H3 * STH_H3 - (self.parameters.recoveryRate + self.parameters.deathRateVH_H3) * ITH_H3
        dRTH_H3 = self.parameters.recoveryRate * ITH_H3
	
	dSTH_B = - (1 - self.parameters.vaccineEfficacyVsInfectionTypical_B) * Lambda_B * STH_B
        dITH_B = Lambda_B * STH_B - (self.parameters.recoveryRate + self.parameters.deathRateVH_B) * ITH_B
        dRTH_B = self.parameters.recoveryRate * ITH_B
	
	#NL
	dSNL_H1 = - (1 - self.parameters.vaccineEfficacyVsInfectionUniversal_H1) * Lambda_H1 * SNL_H1
        dINL_H1 = Lambda_H1 * SNL_H1 - (self.parameters.recoveryRate + self.parameters.deathRateVL_H1) * INL_H1
        dRNL_H1 = self.parameters.recoveryRate * INL_H1
	
	dSNL_H3 = - (1 - self.parameters.vaccineEfficacyVsInfectionUniversal_H3) * Lambda_H3 * SNL_H3
        dINL_H3 = Lambda_H3 * SNL_H3 - (self.parameters.recoveryRate + self.parameters.deathRateVL_H3) * INL_H3
        dRNL_H3 = self.parameters.recoveryRate * INL_H3
	
	dSNL_B = - (1 - self.parameters.vaccineEfficacyVsInfectionUniversal_B) * Lambda_B * SNL_B
        dINL_B = Lambda_B * SNL_B - (self.parameters.recoveryRate + self.parameters.deathRateVL_B) * INL_B
        dRNL_B = self.parameters.recoveryRate * INL_B
	
	#NH
	dSNH_H1 = - (1 - self.parameters.vaccineEfficacyVsInfectionUniversal_H1) * Lambda_H1 * SNH_H1
        dINH_H1 = Lambda_H1 * SNH_H1 - (self.parameters.recoveryRate + self.parameters.deathRateVH_H1) * INH_H1
        dRNH_H1 = self.parameters.recoveryRate * INH_H1
	
	dSNH_H3 = - (1 - self.parameters.vaccineEfficacyVsInfectionUniversal_H3) * Lambda_H3 * SNH_H3
        dINH_H3 = Lambda_H3 * SNH_H3 - (self.parameters.recoveryRate + self.parameters.deathRateVH_H3) * INH_H3
        dRNH_H3 = self.parameters.recoveryRate * INH_H3
	
	dSNH_B = - (1 - self.parameters.vaccineEfficacyVsInfectionUniversal_B) * Lambda_B * SNH_B
        dINH_B = Lambda_B * SNH_B - (self.parameters.recoveryRate + self.parameters.deathRateVH_B) * INH_B
        dRNH_B = self.parameters.recoveryRate * INH_B
	
        
        # Convert meaningful component vectors into a single vector
        dY = numpy.empty(Y.size, dtype = float)
        dY[ 0 : : 54] = dSUL_H1
        dY[ 1 : : 54] = dIUL_H1
        dY[ 2 : : 54] = dRUL_H1
	dY[ 3 : : 54] = dSUL_H3
        dY[ 4 : : 54] = dIUL_H3
        dY[ 5 : : 54] = dRUL_H3
	dY[ 6 : : 54] = dSUL_B
        dY[ 7 : : 54] = dIUL_B
        dY[ 8 : : 54] = dRUL_B
	
	dY[ 9 : : 54] = dSUH_H1
        dY[ 10 : : 54] = dIUH_H1
        dY[ 11 : : 54] = dRUH_H1
	dY[ 12 : : 54] = dSUH_H3
        dY[ 13 : : 54] = dIUH_H3
        dY[ 14 : : 54] = dRUH_H3
	dY[ 15 : : 54] = dSUH_B
        dY[ 16 : : 54] = dIUH_B
        dY[ 17 : : 54] = dRUH_B
	
	dY[ 18 : : 54] = dSTL_H1
        dY[ 19 : : 54] = dITL_H1
        dY[ 20 : : 54] = dRTL_H1
	dY[ 21 : : 54] = dSTL_H3
        dY[ 22 : : 54] = dITL_H3
        dY[ 23 : : 54] = dRTL_H3
	dY[ 24 : : 54] = dSTL_B
        dY[ 25 : : 54] = dITL_B
        dY[ 26 : : 54] = dRTL_B
	
	dY[ 27 : : 54] = dSTH_H1
        dY[ 28 : : 54] = dITH_H1
        dY[ 29 : : 54] = dRTH_H1
	dY[ 30 : : 54] = dSTH_H3
        dY[ 31 : : 54] = dITH_H3
        dY[ 32 : : 54] = dRTH_H3
	dY[ 33 : : 54] = dSTH_B
        dY[ 34 : : 54] = dITH_B
        dY[ 35 : : 54] = dRTH_B
	
	dY[ 36 : : 54] = dSNL_H1
        dY[ 37 : : 54] = dINL_H1
        dY[ 38 : : 54] = dRNL_H1
	dY[ 39 : : 54] = dSNL_H3
        dY[ 40 : : 54] = dINL_H3
        dY[ 41 : : 54] = dRNL_H3
	dY[ 42 : : 54] = dSNL_B
        dY[ 43 : : 54] = dINL_B
        dY[ 44 : : 54] = dRNL_B
	
	dY[ 45 : : 54] = dSNH_H1
        dY[ 46 : : 54] = dINH_H1
        dY[ 47 : : 54] = dRNH_H1
	dY[ 48 : : 54] = dSNH_H3
        dY[ 49 : : 54] = dINH_H3
        dY[ 50 : : 54] = dRNH_H3
	dY[ 51 : : 54] = dSNH_B
        dY[ 52 : : 54] = dINH_B	

        return dY
    
    def resetSolution(self):
        self.hasSolution = False

    def solve(self, tStart = 0., tEnd = None, tStep = 1.):
        if tEnd == None:
            tEnd = self.tMax

        if self.hasSolution:
            TOld  = self.T.copy()
            SUL_H1_Old = self.SUL_H1.copy()
            IUL_H1_Old = self.IUL_H1.copy()
            RUL_H1_Old = self.RUL_H1.copy()
	    SUL_H3_Old = self.SUL_H3.copy()
            IUL_H3_Old = self.IUL_H3.copy()
            RUL_H3_Old = self.RUL_H3.copy()
	    SUL_B_Old = self.SUL_B.copy()
            IUL_B_Old = self.IUL_B.copy()
            RUL_B_Old = self.RUL_B.copy()
	    
	    SUH_H1_Old = self.SUH_H1.copy()
            IUH_H1_Old = self.IUH_H1.copy()
            RUH_H1_Old = self.RUH_H1.copy()
	    SUH_H3_Old = self.SUH_H3.copy()
            IUH_H3_Old = self.IUH_H3.copy()
            RUH_H3_Old = self.RUH_H3.copy()
	    SUH_B_Old = self.SUH_B.copy()
            IUH_B_Old = self.IUH_B.copy()
            RUH_B_Old = self.RUH_B.copy()
	    
	    STL_H1_Old = self.STL_H1.copy()
            ITL_H1_Old = self.ITL_H1.copy()
            RTL_H1_Old = self.RTL_H1.copy()
	    STL_H3_Old = self.STL_H3.copy()
            ITL_H3_Old = self.ITL_H3.copy()
            RTL_H3_Old = self.RTL_H3.copy()
	    STL_B_Old = self.STL_B.copy()
            ITL_B_Old = self.ITL_B.copy()
            RTL_B_Old = self.RTL_B.copy()
	    
	    STH_H1_Old = self.STH_H1.copy()
            ITH_H1_Old = self.ITH_H1.copy()
            RTH_H1_Old = self.RTH_H1.copy()
	    STH_H3_Old = self.STH_H3.copy()
            ITH_H3_Old = self.ITH_H3.copy()
            RTH_H3_Old = self.RTH_H3.copy()
	    STH_B_Old = self.STH_B.copy()
            ITH_B_Old = self.ITH_B.copy()
            RTH_B_Old = self.RTH_B.copy()
	    
	    SNL_H1_Old = self.SNL_H1.copy()
            INL_H1_Old = self.INL_H1.copy()
            RNL_H1_Old = self.RNL_H1.copy()
	    SNL_H3_Old = self.SNL_H3.copy()
            INL_H3_Old = self.INL_H3.copy()
            RNL_H3_Old = self.RNL_H3.copy()
	    SNL_B_Old = self.SNL_B.copy()
            INL_B_Old = self.INL_B.copy()
            RNL_B_Old = self.RNL_B.copy()
	    
	    SNH_H1_Old = self.SNH_H1.copy()
            INH_H1_Old = self.INH_H1.copy()
            RNH_H1_Old = self.RNH_H1.copy()
	    SNH_H3_Old = self.SNH_H3.copy()
            INH_H3_Old = self.INH_H3.copy()
            RNH_H3_Old = self.RNH_H3.copy()
	    SNH_B_Old = self.SNH_B.copy()
            INH_B_Old = self.INH_B.copy()
            RNH_B_Old = self.RNH_B.copy()
            
        # Time vector for solution
        self.T = numpy.hstack((numpy.arange(tStart, tEnd, tStep), tEnd))
        
        # Integrate the ODE
        from scipy.integrate import odeint
        self.Y = odeint(self.RHS,
                        self.Y0.copy(),
                        self.T,
                        mxstep = 1000)
        Z = self.Y.copy()
	
	self.SUL_H1 = Z[:, 0 : : 54]
        self.IUL_H1 = Z[:, 1 : : 54]
        self.RUL_H1 = Z[:, 2 : : 54]
	
	self.SUL_H3 = Z[:, 0 : : 54]
        self.IUL_H3 = Z[:, 1 : : 54]
        self.RUL_H3 = Z[:, 2 : : 54]
	
	self.SUL_B = Z[:, 0 : : 54]
        self.IUL_B = Z[:, 1 : : 54]
        self.RUL_B = Z[:, 2 : : 54]
	
	#UH
	self.SUH_H1 = Z[:, 0 : : 54]
        self.IUH_H1 = Z[:, 1 : : 54]
        self.RUH_H1 = Z[:, 2 : : 54]
	
	self.SUH_H3 = Z[:, 0 : : 54]
        self.IUH_H3 = Z[:, 1 : : 54]
        self.RUH_H3 = Z[:, 2 : : 54]
	
	self.SUH_B = Z[:, 0 : : 54]
        self.IUH_B = Z[:, 1 : : 54]
        self.RUH_B = Z[:, 2 : : 54]
	
	#TL
	self.STL_H1 = Z[:, 0 : : 54]
        self.ITL_H1 = Z[:, 1 : : 54]
        self.RTL_H1 = Z[:, 2 : : 54]
	
	self.STL_H3 = Z[:, 0 : : 54]
        self.ITL_H3 = Z[:, 1 : : 54]
        self.RTL_H3 = Z[:, 2 : : 54]
	
	self.STL_B = Z[:, 0 : : 54]
        self.ITL_B = Z[:, 1 : : 54]
        self.RTL_B = Z[:, 2 : : 54]
	
	#TH
	self.STH_H1 = Z[:, 0 : : 54]
        self.ITH_H1 = Z[:, 1 : : 54]
        self.RTH_H1 = Z[:, 2 : : 54]
	
	self.STH_H3 = Z[:, 0 : : 54]
        self.ITH_H3 = Z[:, 1 : : 54]
        self.RTH_H3 = Z[:, 2 : : 54]
	
	self.STH_B = Z[:, 0 : : 54]
        self.ITH_B = Z[:, 1 : : 54]
        self.RTH_B = Z[:, 2 : : 54]
	
	#NL
	self.SNL_H1 = Z[:, 0 : : 54]
        self.INL_H1 = Z[:, 1 : : 54]
        self.RNL_H1 = Z[:, 2 : : 54]
	
	self.SNL_H3 = Z[:, 0 : : 54]
        self.INL_H3 = Z[:, 1 : : 54]
        self.RNL_H3 = Z[:, 2 : : 54]
	
	self.SNL_B = Z[:, 0 : : 54]
        self.INL_B = Z[:, 1 : : 54]
        self.RNL_B = Z[:, 2 : : 54]
	
	#NH
	self.SNH_H1 = Z[:, 0 : : 54]
        self.INH_H1 = Z[:, 1 : : 54]
        self.RNH_H1 = Z[:, 2 : : 54]
	
	self.SNH_H3 = Z[:, 0 : : 54]
        self.INH_H3 = Z[:, 1 : : 54]
        self.RNH_H3 = Z[:, 2 : : 54]
	
	self.SNH_B = Z[:, 0 : : 54]
        self.INH_B = Z[:, 1 : : 54]
        self.RNH_B = Z[:, 2 : : 54]
	

        if self.hasSolution:
	   
            self.T = numpy.hstack((TOld, self.T))
	    #UL
            self.SUL_H1 = numpy.vstack((SUL_H1_Old, self.SUL_H1))
            self.IUL_H1 = numpy.vstack((IUL_H1_Old, self.IUL_H1))
            self.RUL_H1 = numpy.vstack((RUL_H1_Old, self.RUL_H1))
	    
	    self.SUL_H3 = numpy.vstack((SUL_H3_Old, self.SUL_H3))
            self.IUL_H3 = numpy.vstack((IUL_H3_Old, self.IUL_H3))
            self.RUL_H3 = numpy.vstack((RUL_H3_Old, self.RUL_H3))
	    
	    self.SUL_B = numpy.vstack((SUL_B_Old, self.SUL_B))
            self.IUL_B = numpy.vstack((IUL_B_Old, self.IUL_B))
            self.RUL_B = numpy.vstack((RUL_B_Old, self.RUL_B))
	    
	    #UH
	    self.SUH_H1 = numpy.vstack((SUH_H1_Old, self.SUH_H1))
            self.IUH_H1 = numpy.vstack((IUH_H1_Old, self.IUH_H1))
            self.RUH_H1 = numpy.vstack((RUH_H1_Old, self.RUH_H1))
	    
	    self.SUH_H3 = numpy.vstack((SUH_H3_Old, self.SUH_H3))
            self.IUH_H3 = numpy.vstack((IUH_H3_Old, self.IUH_H3))
            self.RUH_H3 = numpy.vstack((RUH_H3_Old, self.RUH_H3))
	    
	    self.SUH_B = numpy.vstack((SUH_B_Old, self.SUH_B))
            self.IUH_B = numpy.vstack((IUH_B_Old, self.IUH_B))
            self.RUH_B = numpy.vstack((RUH_B_Old, self.RUH_B))
	    
	    #TL
	    self.STL_H1 = numpy.vstack((STL_H1_Old, self.STL_H1))
            self.ITL_H1 = numpy.vstack((ITL_H1_Old, self.ITL_H1))
            self.RTL_H1 = numpy.vstack((RTL_H1_Old, self.RTL_H1))
	    
	    self.STL_H3 = numpy.vstack((STL_H3_Old, self.STL_H3))
            self.ITL_H3 = numpy.vstack((ITL_H3_Old, self.ITL_H3))
            self.RTL_H3 = numpy.vstack((RTL_H3_Old, self.RTL_H3))
	    
	    self.STL_B = numpy.vstack((STL_B_Old, self.STL_B))
            self.ITL_B = numpy.vstack((ITL_B_Old, self.ITL_B))
            self.RTL_B = numpy.vstack((RTL_B_Old, self.RTL_B))
	    
	    #TH
	    self.STH_H1 = numpy.vstack((STH_H1_Old, self.STH_H1))
            self.ITH_H1 = numpy.vstack((ITH_H1_Old, self.ITH_H1))
            self.RTH_H1 = numpy.vstack((RTH_H1_Old, self.RTH_H1))
	    
	    self.STH_H3 = numpy.vstack((STH_H3_Old, self.STH_H3))
            self.ITH_H3 = numpy.vstack((ITH_H3_Old, self.ITH_H3))
            self.RTH_H3 = numpy.vstack((RTH_H3_Old, self.RTH_H3))
	    
	    self.STH_B = numpy.vstack((STH_B_Old, self.STH_B))
            self.ITH_B = numpy.vstack((ITH_B_Old, self.ITH_B))
            self.RTH_B = numpy.vstack((RTH_B_Old, self.RTH_B))
	    
	    #NL
	    self.SNL_H1 = numpy.vstack((SNL_H1_Old, self.SNL_H1))
            self.INL_H1 = numpy.vstack((INL_H1_Old, self.INL_H1))
            self.RNL_H1 = numpy.vstack((RNL_H1_Old, self.RNL_H1))
	    
	    self.SNL_H3 = numpy.vstack((SNL_H3_Old, self.SNL_H3))
            self.INL_H3 = numpy.vstack((INL_H3_Old, self.INL_H3))
            self.RNL_H3 = numpy.vstack((RNL_H3_Old, self.RNL_H3))
	    
	    self.SNL_B = numpy.vstack((SNL_B_Old, self.SNL_B))
            self.INL_B = numpy.vstack((INL_B_Old, self.INL_B))
            self.RNL_B = numpy.vstack((RNL_B_Old, self.RNL_B))
	    
	    #NH
	    self.SNH_H1 = numpy.vstack((SNH_H1_Old, self.SNH_H1))
            self.INH_H1 = numpy.vstack((INH_H1_Old, self.INH_H1))
            self.RNH_H1 = numpy.vstack((RNH_H1_Old, self.RNH_H1))
	    
	    self.SNH_H3 = numpy.vstack((SNH_H3_Old, self.SNH_H3))
            self.INH_H3 = numpy.vstack((INH_H3_Old, self.INH_H3))
            self.RNH_H3 = numpy.vstack((RNH_H3_Old, self.RNH_H3))
	    
	    self.SNH_B = numpy.vstack((SNH_B_Old, self.SNH_B))
            self.INH_B = numpy.vstack((INH_B_Old, self.INH_B))
            self.RNH_B = numpy.vstack((RNH_B_Old, self.RNH_B))
	    
	    
        self.hasSolution = True

    def updateStats(self):
	self.NUL_H1 = self.SUL_H1 +  self.IUL_H1 + self.RUL_H1
	self.NUL_H3 = self.SUL_H3 +  self.IUL_H3 + self.RUL_H3
        self.NUL_B = self.SUL_B +  self.IUL_B + self.RUL_B
	
	self.NUH_H1 = self.SUH_H1 +  self.IUH_H1 + self.RUH_H1
	self.NUH_H3 = self.SUH_H3 +  self.IUH_H3 + self.RUH_H3
        self.NUH_B = self.SUH_B +  self.IUH_B + self.RUH_B
	
	self.NTL_H1 = self.STL_H1 +  self.ITL_H1 + self.RTL_H1
	self.NTL_H3 = self.STL_H3 +  self.ITL_H3 + self.RTL_H3
        self.NTL_B = self.STL_B +  self.ITL_B + self.RTL_B
	
	self.NTH_H1 = self.STH_H1 +  self.ITH_H1 + self.RTH_H1
	self.NTH_H3 = self.STH_H3 +  self.ITH_H3 + self.RTH_H3
        self.NTH_B = self.STH_B +  self.ITH_B + self.RTH_B
	
	self.NNL_H1 = self.SNL_H1 +  self.INL_H1 + self.RNL_H1
	self.NNL_H3 = self.SNL_H3 +  self.INL_H3 + self.RNL_H3
        self.NNL_B = self.SNL_B +  self.INL_B + self.RNL_B
	
	self.NNH_H1 = self.SNH_H1 +  self.INH_H1 + self.RNH_H1
	self.NNH_H3 = self.SNH_H3 +  self.INH_H3 + self.RNH_H3
        self.NNH_B = self.SNH_B +  self.INH_B + self.RNH_B
	
	self.NU = self.NUL_H1 + self.NUL_H3 + self.NUL_B + self.NUH_H1 + self.NUH_H3 + self.NUH_B
	self.NT = self.NTL_H1 + self.NTL_H3 + self.NTL_B + self.NTH_H1 + self.NTH_H3 + self.NTH_B
	self.NN = self.NNL_H1 + self.NNL_H3 + self.NNL_B + self.NNH_H1 + self.NNH_H3 + self.NNH_B
	
        self.N  = self.NU + self.NT + self.NN
	
	self.infectionsUL_H1 = self.NUL_H1[0, :] - self.SUL_H1[-1, :]
	self.infectionsUL_H3 = self.NUL_H3[0, :] - self.SUL_H3[-1, :]
	self.infectionsUL_B = self.NUL_B[0, :] - self.SUL_B[-1, :]
	
	self.infectionsUH_H1 = self.NUH_H1[0, :] - self.SUH_H1[-1, :]
	self.infectionsUH_H3 = self.NUH_H3[0, :] - self.SUH_H3[-1, :]
	self.infectionsUH_B = self.NUH_B[0, :] - self.SUH_B[-1, :]
	
	self.infectionsTL_H1 = self.NTL_H1[0, :] - self.STL_H1[-1, :]
	self.infectionsTL_H3 = self.NTL_H3[0, :] - self.STL_H3[-1, :]
	self.infectionsTL_B = self.NTL_B[0, :] - self.STL_B[-1, :]
	
	self.infectionsTH_H1 = self.NTH_H1[0, :] - self.STH_H1[-1, :]
	self.infectionsTH_H3 = self.NTH_H3[0, :] - self.STH_H3[-1, :]
	self.infectionsTH_B = self.NTH_B[0, :] - self.STH_B[-1, :]
	
	self.infectionsNL_H1 = self.NNL_H1[0, :] - self.SNL_H1[-1, :]
	self.infectionsNL_H3 = self.NNL_H3[0, :] - self.SNL_H3[-1, :]
	self.infectionsNL_B = self.NNL_B[0, :] - self.SNL_B[-1, :]
	
	self.infectionsNH_H1 = self.NNH_H1[0, :] - self.SNH_H1[-1, :]
	self.infectionsNH_H3 = self.NNH_H3[0, :] - self.SNH_H3[-1, :]
	self.infectionsNH_B = self.NNH_B[0, :] - self.SNH_B[-1, :]
	

        self.infectionsUL = self.infectionsUL_H1 + self.infectionsUL_H3 + self.infectionsUL_B
	self.infectionsUH = self.infectionsUH_H1 + self.infectionsUH_H3 + self.infectionsUH_B
	self.infectionsTL = self.infectionsTL_H1 + self.infectionsTL_H3 + self.infectionsTL_B
	self.infectionsTH = self.infectionsTH_H1 + self.infectionsTH_H3 + self.infectionsTH_B
	self.infectionsNL = self.infectionsNL_H1 + self.infectionsNL_H3 + self.infectionsNL_B
	self.infectionsNH = self.infectionsNH_H1 + self.infectionsNH_H3 + self.infectionsNH_B
	
	self.infectionsL_H1 = self.infectionsUL_H1 + self.infectionsTL_H1 + self.infectionsNL_H1
	self.infectionsL_H3 = self.infectionsUL_H3 + self.infectionsTL_H3 + self.infectionsNL_H3
	self.infectionsL_B = self.infectionsUL_B + self.infectionsTL_B + self.infectionsNL_B
	
	self.infectionsH_H1 = self.infectionsUH_H1 + self.infectionsTH_H1 + self.infectionsNH_H1
	self.infectionsH_H3 = self.infectionsUH_H3 + self.infectionsTH_H3 + self.infectionsNH_H3
	self.infectionsH_B = self.infectionsUH_B + self.infectionsTH_B + self.infectionsNH_B

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
        
        self.deathsUL_H1 = self.NUL_H1[0, :] - self.NUL_H1[-1, :]
	self.deathsUL_H3 = self.NUL_H3[0, :] - self.NUL_H3[-1, :]
	self.deathsUL_B = self.NUL_B[0, :] - self.NUL_B[-1, :]
	
	self.deathsUH_H1 = self.NUH_H1[0, :] - self.NUH_H1[-1, :]
	self.deathsUH_H3 = self.NUH_H3[0, :] - self.NUH_H3[-1, :]
	self.deathsUH_B = self.NUH_B[0, :] - self.NUH_B[-1, :]
	
	self.deathsTL_H1 = self.NTL_H1[0, :] - self.NTL_H1[-1, :]
	self.deathsTL_H3 = self.NTL_H3[0, :] - self.NTL_H3[-1, :]
	self.deathsTL_B = self.NTL_B[0, :] - self.NTL_B[-1, :]
	
	self.deathsTH_H1 = self.NTH_H1[0, :] - self.NTH_H1[-1, :]
	self.deathsTH_H3 = self.NTH_H3[0, :] - self.NTH_H3[-1, :]
	self.deathsTH_B = self.NTH_B[0, :] - self.NTH_B[-1, :]
	
	self.deathsNL_H1 = self.NNL_H1[0, :] - self.NNL_H1[-1, :]
	self.deathsNL_H3 = self.NNL_H3[0, :] - self.NNL_H3[-1, :]
	self.deathsNL_B = self.NNL_B[0, :] - self.NNL_B[-1, :]
	
	self.deathsNH_H1 = self.NNH_H1[0, :] - self.NNH_H1[-1, :]
	self.deathsNH_H3 = self.NNH_H3[0, :] - self.NNH_H3[-1, :]
	self.deathsNH_B = self.NNH_B[0, :] - self.NNH_B[-1, :]
	
	
	self.deathsUL = self.deathsUL_H1 + self.deathsUL_H3 + self.deathsUL_B
	self.deathsUH = self.deathsUH_H1 + self.deathsUH_H3 + self.deathsUH_B
	self.deathsTL = self.deathsTL_H1 + self.deathsTL_H3 + self.deathsTL_B
	self.deathsTH = self.deathsTH_H1 + self.deathsTH_H3 + self.deathsTH_B
	self.deathsNL = self.deathsNL_H1 + self.deathsNL_H3 + self.deathsNL_B
	self.deathsNH = self.deathsNH_H1 + self.deathsNH_H3 + self.deathsNH_B

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
        
        
    def simulate(self):
        self.updateIC()
        self.solve()
        self.updateStats()

    def updateProportionVaccinated(self, PVPWVal, nVacTypes):
        # Update propotion vaccinated

	 # Convert flat vector to 2-D array
        if numpy.ndim(PVPWVal) != 2:
	  
            PVPWVal = numpy.asarray(PVPWVal).reshape(
                (nVacTypes,
                 self.parameters.proportionVaccinatedLLength + self.parameters.proportionVaccinatedHLength))

   
	
        self.parameters.proportionVaccinatedTLPW.values = PVPWVal[0][: self.parameters.proportionVaccinatedLLength]
	self.parameters.proportionVaccinatedTHPW.values = PVPWVal[0][self.parameters.proportionVaccinatedLLength:]
	self.parameters.proportionVaccinatedNLPW.values = PVPWVal[1][: self.parameters.proportionVaccinatedLLength]
	self.parameters.proportionVaccinatedNHPW.values = PVPWVal[1][self.parameters.proportionVaccinatedLLength:]

	## extend to full ages groups. Proportions calculated by multiplying PVPWVal 
	##values with the matrix defined in S.130
	
	self.parameters.proportionVaccinatedTL = self.parameters.proportionVaccinatedTLPW.full(self.parameters.ages)
	self.parameters.proportionVaccinatedTH = self.parameters.proportionVaccinatedTHPW.full(self.parameters.ages)
	self.parameters.proportionVaccinatedNL = self.parameters.proportionVaccinatedNLPW.full(self.parameters.ages)
	self.parameters.proportionVaccinatedNH = self.parameters.proportionVaccinatedNHPW.full(self.parameters.ages)
	
	
	if self.hasSolution:
	   
            vacsUsedTL = self.parameters.proportionVaccinatedTL* IC[0]
	    vacsUsedTH = self.parameters.proportionVaccinatedTH* IC[0]
	    vacsUsedTypical = vacsUsedTL + vacsUsedTH
	    vacsUsedNL = self.parameters.proportionVaccinatedNL* IC[0]
	    vacsUsedNH = self.parameters.proportionVaccinatedNH* IC[0]
	    vacsUsedUniversal = vacsUsedNL + vacsUsedNH

        else:
	    vacsUsedTypical = (self.parameters.proportionVaccinatedTL * self.parameters.population) + (self.parameters.proportionVaccinatedTH * self.parameters.population)
            vacsUsedUniversal = (self.parameters.proportionVaccinatedNL * self.parameters.population) +  (self.parameters.proportionVaccinatedNH * self.parameters.population)
	
	   
        # Update initial condition for ODEs
        self.updateIC()

        return vacsUsedTypical, vacsUsedUniversal
    
    def vaccinated_output(self):
        return list(self.parameters.proportionVaccinatedTypical), list(self.parameters.proportionVaccinatedUniversal), [(a*b) for (a,b) in zip(self.parameters.proportionVaccinatedTypical, self.parameters.population)],[(a*b) for (a,b) in zip(self.parameters.proportionVaccinatedUniversal, self.parameters.population)]
    
    
        
    def simulateWithVaccine(self, PVPWVals, vacEfficacy):
	
        nVacTypes = len(vacEfficacy)

        self.resetSolution()

        # Vaccinate the population
	vacsUsedTypical, vacsUsedUniversal = self.updateProportionVaccinated(PVPWVals, nVacTypes)
	
	if min(list(vacsUsedUniversal)) <0 and min(list(PVPWVals)) >0 : print ("check!!!!"), PVPWVals, nVacTypes, vacsUsedUniversal
        tEnd = self.tMax
	tStart = self.tMin
	self.solve(tStart = tStart, tEnd = tEnd)
	
	#print self.IU.sum(axis=1)
	import matplotlib.pyplot as plt
	times = [num for num in xrange(tEnd+1)]
	plt.plot(times, self.IUL_H3.sum(axis=1), color = "red", linestyle= "-")
	plt.plot(times, self.IUH_H3.sum(axis=1), color = "red", linestyle= "--")
	plt.plot(times, self.ITL_H3.sum(axis=1), color = "blue", linestyle = "-")
	plt.plot(times, self.ITH_H3.sum(axis=1), color = "blue", linestyle = "--")
	plt.plot(times, self.INL_H3.sum(axis=1), color = "green", linestyle = "-")
	plt.plot(times, self.INH_H3.sum(axis=1), color = "green", linestyle = "--")
		 
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


