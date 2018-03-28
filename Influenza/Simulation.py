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
        self.Y0 = numpy.zeros(18 * self.parameters.ages.size)

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
	    typical_vaccination = [ 0.70,  0.56,  0.56, 0.5152, 0.34, 0.34,0.34,0.34,0.34,0.34,
				   0.46,  0.46, 0.46, 0.653, 0.653, 0.653]
	    
	    typical_proportion = self.vacNumbers[0]/(1. *sum(self.vacNumbers))
	    universal_proportion = 1 - typical_proportion
	    
	    typical_vax_coverage = [typical_proportion * num for num in typical_vaccination]
	    universal_vax_coverage = [(a-b) for (a,b) in zip(typical_vaccination, typical_vax_coverage)]
	    
	    #print ("check vaccine numbers"), sum([(a+b)*c for (a,b,c) in zip(typical_vax_coverage, universal_vax_coverage, self.parameters.population[1:])])		   

	    
	    return typical_vax_coverage + universal_vax_coverage	
	    

    def computeR0(self):
        return self.parameters.computeR0()

    def getLastValues(self):
        return (self.SUL[-1, :], self.IUL[-1, :], self.RUL[-1, :],
		self.SUH[-1, :], self.IUH[-1, :], self.RUH[-1, :],
		self.STL[-1, :], self.ITL[-1, :], self.RTL[-1, :],
		self.STH[-1, :], self.ITH[-1, :], self.RTH[-1, :],
		self.SNL[-1, :], self.INL[-1, :], self.RNL[-1, :],
                self.SNH[-1, :], self.INH[-1, :], self.RNH[-1, :])

    def updateIC(self):
        if not self.hasSolution:
            # S
	    ## SUL
            self.Y0[ 0 : : 18] = \
                     (1 - self.parameters.proportionVaccinatedL - self.parameters.proportionVaccinatedH) * self.parameters.population * (1 - self.parameters.proportionHighRisk)
	    ## SUH
            self.Y0[ 3 : : 18] = \
                     (1 - self.parameters.proportionVaccinatedL - self.parameters.proportionVaccinatedH) * self.parameters.population * self.parameters.proportionHighRisk 

	    ## STL
            self.Y0[ 6 : : 18] = \
                     self.parameters.proportionVaccinatedTL * self.parameters.population * (1 - self.parameters.proportionHighRisk)
	    ## STH
            self.Y0[ 9 : : 18] = \
                     self.parameters.proportionVaccinatedTH * self.parameters.population * self.parameters.proportionHighRisk
	    
	    ## SNL
            self.Y0[ 12 : : 18] = \
                     self.parameters.proportionVaccinatedNL * self.parameters.population  * (1 - self.parameters.proportionHighRisk)
	    ## SNH
            self.Y0[ 15 : : 18] = \
                     self.parameters.proportionVaccinatedNH * self.parameters.population * self.parameters.proportionHighRisk


            # I: Add a single infectious person in each age
	    
	    self.Y0[ 1 : : 18] = numpy.full(self.parameters.ages.size,1)
            self.Y0[ 4 : : 18] = numpy.full(self.parameters.ages.size, 1)
	    self.Y0[ 7 : : 18] = numpy.full(self.parameters.ages.size, 1)
	    self.Y0[ 10 : : 18] = numpy.full(self.parameters.ages.size, 1)
	    self.Y0[ 13 : : 18] = numpy.full(self.parameters.ages.size, 1)
	    self.Y0[ 16 : : 18] = numpy.full(self.parameters.ages.size, 1)
	    

            # S: Remove those new infectious people from the susceptibles
            self.Y0[ 0 : : 18] -= self.Y0[ 1 : : 18]
            self.Y0[ 3 : : 18] -= self.Y0[ 4 : : 18]
	    self.Y0[ 6 : : 18] -= self.Y0[ 7 : : 18]
	    self.Y0[ 9 : : 18] -= self.Y0[ 10 : : 18]
	    self.Y0[ 12 : : 18] -= self.Y0[13 : : 18]
	    self.Y0[ 15 : : 18] -= self.Y0[16 : : 18]

            # R (RU, RVL, RVH)
            self.Y0[ 2 : : 18] = 0.
            self.Y0[ 5 : : 18] = 0.
	    self.Y0[ 8 : : 18] = 0.
	    self.Y0[ 11 : : 18] = 0.
	    self.Y0[ 14 : : 18] = 0.
	    self.Y0[ 17 : : 18] = 0.
	    print ("check!!!!"),  self.Y0[ 1 : : 18],  self.Y0[ 4 : : 18],  self.Y0[ 7 : : 18],  self.Y0[ 10 : : 18],  self.Y0[ 13 : : 18],  self.Y0[ 16 : : 18]

        else:
	    
            SUL, IUL, RUL, SUH, IUH, RUH, STL,  ITL, RTL, STH,  ITH, RTH, SNL, INL, RNL, SNH, INH, RNH = self.getLastValues()
            
            self.Y0[ 0 : : 18] = (1 - self.parameters.proportionVaccinatedL) * SUL
	    self.Y0[ 3 : : 18] = (1 - self.parameters.proportionVaccinatedH) * SUH
            self.Y0[ 6 : : 18] = STL + self.parameters.proportionVaccinatedTL * SUL
	    self.Y0[ 9 : : 18] = STH + self.parameters.proportionVaccinatedTH * SUH
	    self.Y0[ 12 : : 18] = SNL + self.parameters.proportionVaccinatedNL * SUL
	    self.Y0[ 15 : : 18] = SNH + self.parameters.proportionVaccinatedNH * SUH
	    
            self.Y0[ 1 : : 18] = IUL
            self.Y0[ 4 : : 18] = IVH
	    self.Y0[ 7 : : 18] = ITL
	    self.Y0[ 10 : : 18] = ITH
	    self.Y0[ 13 : : 18] = INL
	    self.Y0[ 16 : : 18] = INH
	    
            self.Y0[ 2 : : 18] = RUL
            self.Y0[ 5 : : 18] = RUH
	    self.Y0[ 8 : : 18] = RTL
	    self.Y0[ 11 : : 18] = RTH
	    self.Y0[ 14 : : 18] = RNL
	    self.Y0[ 17 : : 18] = RNH

    def RHS(self, Y, t):
        '''
        SEIR model with multiple host types.
        
        This function gives the right-hand sides of the ODEs.
        '''
        
        # Convert vector to meaningful component vectors

	SUL = Y[ 0 : : 18]
        IUL = Y[ 1 : : 18]
        RUL = Y[ 2 : : 18]
	SUH = Y[ 3 : : 18]
        IUH = Y[ 4 : : 18]
        RUH = Y[ 5 : : 18]
        STL = Y[ 6 : : 18]
        ITL = Y[ 7 : : 18]
        RTL = Y[ 8 : : 18]
	STH = Y[ 9 : : 18]
        ITH = Y[ 10 : : 18]
        RTH = Y[ 11: : 18]
	SNL = Y[ 12 : : 18]
        INL = Y[ 13 : : 18]
        RNL = Y[ 14 : : 18]
	SNH = Y[ 15 : : 18]
        INH = Y[ 16 : : 18]
        RNH = Y[ 17: : 18]
        
        N = sum(SUL + IUL + RUL + SUH + IUH + RUH + STL + ITL + RTL+ STH + ITH + RTH + SNL + INL + RNL+ SNH + INH + RNH)
      
        # The force of infection
        Lambda = self.parameters.transmissionScaling * self.parameters.susceptibility \
                 * numpy.dot(self.parameters.contactMatrix,
                             self.parameters.transmissibility * (IUL + IUH+ ITL + ITH + INL + INH)) / N
        
        # The right-hand sides
        dSUL = - Lambda * SUL
        dIUL = Lambda * SUL - (self.parameters.recoveryRate + self.parameters.deathRateUL) * IUL
        dRUL = self.parameters.recoveryRate * IUL
	
	dSUH = - Lambda * SUH
        dIUH = Lambda * SUH - (self.parameters.recoveryRate + self.parameters.deathRateUH) * IUH
        dRUH = self.parameters.recoveryRate * IUH
        
	dSTL = - (1 - self.parameters.vaccineEfficacyVsInfectionTypical) * Lambda * STL
        dITL = Lambda * STL  - (self.parameters.recoveryRate + self.parameters.deathRateVL) * ITL
        dRTL = self.parameters.recoveryRate * ITL
	
	dSTH = - (1 - self.parameters.vaccineEfficacyVsInfectionTypical) * Lambda * STH
        dITH = Lambda * STH  - (self.parameters.recoveryRate + self.parameters.deathRateVH) * ITH
        dRTH = self.parameters.recoveryRate * ITH
	
	dSNL = - (1 - self.parameters.vaccineEfficacyVsInfectionUniversal) * Lambda * SNL
        dINL = Lambda * SNL  - (self.parameters.recoveryRate + self.parameters.deathRateVL) * INL
        dRNL = self.parameters.recoveryRate * INL
	
	dSNH = - (1 - self.parameters.vaccineEfficacyVsInfectionUniversal) * Lambda * SNH
        dINH = Lambda * SNH  - (self.parameters.recoveryRate + self.parameters.deathRateVH) * INH
        dRNH = self.parameters.recoveryRate * INH
        
        
        # Convert meaningful component vectors into a single vector
        dY = numpy.empty(Y.size, dtype = float)
        dY[ 0 : : 18] = dSUL
        dY[ 1 : : 18] = dIUL
        dY[ 2 : : 18] = dRUL
	dY[ 3 : : 18] = dSUH
        dY[ 4 : : 18] = dIUH
        dY[ 5 : : 18] = dRUH
        dY[ 6 : : 18] = dSTL
        dY[ 7 : : 18] = dITL
        dY[ 8 : : 18] = dRTL
        dY[ 9 : : 18] = dSTH
        dY[ 10 : : 18] = dITH
        dY[ 11 : : 18] = dRTH
	dY[ 12 : : 18] = dSNL
        dY[ 13 : : 18] = dINL
        dY[ 14 : : 18] = dRNL
        dY[ 15 : : 18] = dSNH
        dY[ 16 : : 18] = dINH
        dY[ 17 : : 18] = dRNH
        return dY
    
    def resetSolution(self):
        self.hasSolution = False

    def solve(self, tStart = 0., tEnd = None, tStep = 1.):
        if tEnd == None:
            tEnd = self.tMax

        if self.hasSolution:
            TOld  = self.T.copy()
            SULOld = self.SUL.copy()
            IULOld = self.IUL.copy()
            RULOld = self.RUL.copy()
	    SUHOld = self.SUH.copy()
            IUHOld = self.IUH.copy()
            RUHOld = self.RUH.copy()
            STLOld = self.STL.copy()
            ITLOld = self.ITL.copy()
            RTLOld = self.RTL.copy()
            STHOld = self.STH.copy()
            ITHOld = self.ITH.copy()
            RTHOld = self.RTH.copy()
	    SNLOld = self.SNL.copy()
            INLOld = self.INL.copy()
            RNLOld = self.RNL.copy()
            SNHOld = self.SNH.copy()
            INHOld = self.INH.copy()
            RNHOld = self.RNH.copy()
            
            
        # Time vector for solution
        self.T = numpy.hstack((numpy.arange(tStart, tEnd, tStep), tEnd))
        
        # Integrate the ODE
        from scipy.integrate import odeint
        self.Y = odeint(self.RHS,
                        self.Y0.copy(),
                        self.T,
                        mxstep = 1000)
        Z = self.Y.copy()
	self.SUL = Z[:,  0 : : 18]
        self.IUL = Z[:,  1 : : 18]
        self.RUL = Z[:,  2 : : 18]
	self.SUH = Z[:,  3 : : 18]
        self.IUH = Z[:,  4 : : 18]
        self.RUH = Z[:,  5 : : 18]
        self.STL = Z[:,  6 : : 18]
        self.ITL = Z[:,  7 : : 18]
        self.RTL = Z[:,  8 : : 18]
        self.STH = Z[:,  9 : : 18]
        self.ITH = Z[:,  10 : : 18]
        self.RTH = Z[:,  11 : : 18]
	self.SNL = Z[:,  12 : : 18]
        self.INL = Z[:,  13 : : 18]
        self.RNL = Z[:,  14 : : 18]
        self.SNH = Z[:,  15 : : 18]
        self.INH = Z[:,  16 : : 18]
        self.RNH = Z[:,  17 : : 18]

        if self.hasSolution:
	   
            self.T = numpy.hstack((TOld, self.T))

            self.SUL = numpy.vstack((SULOld, self.SUL))
            self.IUL = numpy.vstack((IULOld, self.IUL))
            self.RUL = numpy.vstack((RULOld, self.RUL))
	    self.SUH = numpy.vstack((SUHOld, self.SUL))
            self.IUH = numpy.vstack((IUHOld, self.IUL))
            self.RUH = numpy.vstack((RUHOld, self.RUL))
            self.STL = numpy.vstack((STLOld, self.STL))
            self.ITL = numpy.vstack((ITLOld, self.ITL))
            self.RTL = numpy.vstack((RTLOld, self.RTL))
            self.STH = numpy.vstack((STHOld, self.STH))
            self.ITH = numpy.vstack((ITHOld, self.ITH))
            self.RTH = numpy.vstack((RTHOld, self.RTH))
	    self.SNL = numpy.vstack((SNLOld, self.SNL))
            self.INL = numpy.vstack((INLOld, self.INL))
            self.RNL = numpy.vstack((RNLOld, self.RNL))
            self.SNH = numpy.vstack((SNHOld, self.SNH))
            self.INH = numpy.vstack((INHOld, self.INH))
            self.RNH = numpy.vstack((RNHOld, self.RNH))

        self.hasSolution = True

    def updateStats(self):
        self.NUL = self.SUL +  self.IUL + self.RUL
	self.NUH = self.SUH +  self.IUH + self.RUH
	self.NTL = self.STL +  self.ITL + self.RTL
	self.NTH = self.STH +  self.ITH + self.RTH
	self.NNL = self.SNL +  self.INL + self.RNL
	self.NNH = self.SNH +  self.INH + self.RNH
	
	self.NU = self.NUL + self.NUH
	self.NVL = self.NUL + self.NTL+ self.NNL
        self.NVH = self.NUH + self.NTH + self.NNH
	self.NV = self.NVL + self.NVH
        self.N  = self.NU + self.NV

        self.infectionsUL = self.NUL[0, :] - self.SUL[-1, :]
	self.infectionsUH = self.NUH[0, :] - self.SUH[-1, :]
        self.infectionsTL = self.NTL[0, :] - self.STL[-1, :]
	self.infectionsTH = self.NTH[0, :] - self.STH[-1, :]
	self.infectionsNL = self.NNL[0, :] - self.SNL[-1, :]
	self.infectionsNH = self.NNH[0, :] - self.SNH[-1, :]

        # Find duplicate times: these are where vaccination occurs
        for i in numpy.arange(len(self.T)).compress(numpy.diff(self.T) == 0):
            # Update for vaccinated people
            self.infectionsUL += self.SUL[i + 1, :] - self.SUL[i, :]
	    self.infectionsUH += self.SUH[i + 1, :] - self.SUH[i, :]
            self.infectionsTL += self.STL[i + 1, :] - self.STL[i, :]
	    self.infectionsTH += self.STH[i + 1, :] - self.STH[i, :]
	    self.infectionsNL += self.SNL[i + 1, :] - self.SNL[i, :]
	    self.infectionsNH += self.SNH[i + 1, :] - self.SNH[i, :]

	self.infectionsL  = self.infectionsUL + self.infectionsTL + self.infectionsNL
        self.infectionsH  = self.infectionsUH + self.infectionsTH + self.infectionsNH
	self.infectionsU =  self.infectionsUL + self.infectionsUH
	self.infectionsV  = self.infectionsTL + self.infectionsTH + self.infectionsNL + self.infectionsNH        
	self.infections  = self.infectionsU + self.infectionsV
        self.totalInfections = self.infections.sum()
        
	
	
	self.hospitalizationsL = self.infectionsL * self.parameters.caseHospitalizationL
	self.hospitalizationsH = self.infectionsH * self.parameters.caseHospitalizationH
        self.hospitalizations = self.hospitalizationsL + self.hospitalizationsH
	self.totalHospitalizations = self.hospitalizations.sum()
        
        self.deathsUL = self.NUL[0, :] - self.NUL[-1, :]
	self.deathsUH = self.NUH[0, :] - self.NUH[-1, :]
        self.deathsTL = self.NTL[0, :] - self.NTL[-1, :]
	self.deathsTH = self.NTH[0, :] - self.NTH[-1, :]
	self.deathsNL = self.NNL[0, :] - self.NNL[-1, :]
	self.deathsNH = self.NNH[0, :] - self.NNH[-1, :]

        # Find duplicate times: these are where vaccination occurs
        for i in numpy.arange(len(self.T)).compress(numpy.diff(self.T) == 0):
            # Update for vaccinated people
            self.deathsUL += self.NUL[i + 1, :] - self.NUL[i, :]
	    self.deathsUH += self.NUH[i + 1, :] - self.NUH[i, :]
            self.deathsTL += self.NTL[i + 1, :] - self.NTL[i, :]
	    self.deathsTH += self.NTH[i + 1, :] - self.NTH[i, :]
	    self.deathsNL += self.NNL[i + 1, :] - self.NNL[i, :]
	    self.deathsNH += self.NNH[i + 1, :] - self.NNH[i, :]


	self.deathsL = self.dealthsUL + self.deathsTL + self.deathsNL
	self.deathsH = self.dealthsUH + self.deathsTH + self.deathsNH
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
                 self.parameters.proportionVaccinatedLength))

   
	
        self.parameters.proportionVaccinatedTypicalPW.values = PVPWVal[0]
	self.parameters.proportionVaccinatedUniversalPW.values = PVPWVal[1]

	## extend to full ages groups. Proportions calculated by multiplying PVPWVal 
	##values with the matrix defined in S.130
	
	self.parameters.proportionVaccinatedTypical = self.parameters.proportionVaccinatedTypicalPW.full(self.parameters.ages)
	self.parameters.proportionVaccinatedUniversal = self.parameters.proportionVaccinatedUniversalPW.full(self.parameters.ages)
	
	
	
	if self.hasSolution:
	   
            vacsUsedTypical = (self.parameters.proportionVaccinatedTypical * IC[0])
	    vacsUsedUniversal = (self.parameters.proportionVaccinatedUniversal * IC[0])

        else:
	    vacsUsedTypical = (self.parameters.proportionVaccinatedTypical
                        * self.parameters.population)

            vacsUsedUniversal = (self.parameters.proportionVaccinatedUniversal
                        * self.parameters.population)
	
	   
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
	
	print self.IU.sum(axis=1)
	import matplotlib.pyplot as plt
	times = [num for num in xrange(tEnd+1)]
	plt.plot(times, self.IUL.sum(axis=1), color = "red", linestyle= "-")
	plt.plot(times, self.IUH.sum(axis=1), color = "red", linestyle= "--")
	plt.plot(times, self.ITL.sum(axis=1), color = "blue", linestyle = "-")
	plt.plot(times, self.ITH.sum(axis=1), color = "blue", linestyle = "--")
	plt.plot(times, self.INL.sum(axis=1), color = "green", linestyle = "-")
	plt.plot(times, self.INH.sum(axis=1), color = "green", linestyle = "--")
		 
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


