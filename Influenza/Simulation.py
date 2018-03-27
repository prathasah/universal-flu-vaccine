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
        self.Y0 = numpy.zeros(9 * self.parameters.ages.size)

        self.hasSolution = False
	
    def comupte_random_vaccination(self, doses):
	
	    self.vacNumbers = doses
	    low_vax_doses_condition = False
	    high_vax_doses_condition = False
	    ## chose a random between between 0 and min of (population size of subgroup and doses remaining). Reqd for when max doses = 0
	    low_vax_remaining = self.vacNumbers[0]
	    vacsUsedLow = {}
	    ##randomize age groups that go first to avoid high values chosen for initial age groups
	    random_age_group_list = [num for num in range(1,17)]
	    random.shuffle(random_age_group_list)
	    for num in random_age_group_list:
		vacsUsedLow[num] = random.randint(0, min(self.parameters.population[num], low_vax_remaining))
		##find the last value by substracting the sum of vacsUsedLow
		low_vax_remaining = self.vacNumbers[0] - sum(vacsUsedLow.values())
	    population_unvaccinated = [self.parameters.population[0]] + [(a-b) for (a,b) in zip(self.parameters.population[1:], [vacsUsedLow[age] for age in range(1,17)])]
	   
	    
	    high_vax_remaining = self.vacNumbers[1]
	    vacsUsedHigh = {}
	    random.shuffle(random_age_group_list)
	    for num in random_age_group_list:
		vacsUsedHigh[num] = random.randint(0, min(population_unvaccinated[num],high_vax_remaining))
		##find the last value by substracting the sum of vacsUsedLow
		high_vax_remaining = self.vacNumbers[1] - sum(vacsUsedHigh.values())
		
	    
	    return [a/(1.*b) for (a,b) in zip([vacsUsedLow[age] for age in range(1,17)], self.parameters.population[1:])] +  [a/(1.*b) for (a,b) in zip([vacsUsedHigh[age] for age in range(1,17)], self.parameters.population[1:])]
	    
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
	    
	    low_proportion = self.vacNumbers[0]/(1. *sum(self.vacNumbers))
	    high_proportion = 1 - low_proportion
	    
	    low_vax_coverage = [low_proportion * num for num in typical_vaccination]
	    high_vax_coverage = [(a-b) for (a,b) in zip(typical_vaccination, low_vax_coverage)]
	    
	    #print ("check vaccine numbers"), sum([(a+b)*c for (a,b,c) in zip(low_vax_coverage, high_vax_coverage, self.parameters.population[1:])])		   

	    
	    return low_vax_coverage + high_vax_coverage	
	    

    def computeR0(self):
        return self.parameters.computeR0()

    def getLastValues(self):
        return (self.SU[-1, :], self.IU[-1, :], self.RU[-1, :],
		self.SVL[-1, :], self.IVL[-1, :], self.RVL[-1, :],
                self.SVH[-1, :], self.IVH[-1, :], self.RVH[-1, :])

    def updateIC(self):
        if not self.hasSolution:
            # S
	    ## SU
	   
            self.Y0[ 0 : : 9] = \
                     (1 - self.parameters.proportionVaccinatedLow - self.parameters.proportionVaccinatedHigh) * self.parameters.population 

	    ## SVL
            self.Y0[ 3 : : 9] = \
                     self.parameters.proportionVaccinatedLow * self.parameters.population 
	    ## SVH
            self.Y0[ 6 : : 9] = \
                     self.parameters.proportionVaccinatedHigh * self.parameters.population 


            # I: Add a single infectious person in each age
	    
	    self.Y0[ 1 : : 9] = numpy.full(self.parameters.ages.size,1)
            self.Y0[ 4 : : 9] = numpy.full(self.parameters.ages.size, 1)
	    self.Y0[ 7 : : 9] = numpy.full(self.parameters.ages.size, 1)
	    
            #self.Y0[ 2 : : 12] = (1 - self.parameters.proportionVaccinatedLow - self.parameters.proportionVaccinatedHigh)
            #self.Y0[ 6 : : 12] = self.parameters.proportionVaccinatedLow
	    #self.Y0[ 10 : : 12] = self.parameters.proportionVaccinatedHigh

            # S: Remove those new infectious people from the susceptibles
            self.Y0[ 0 : : 9] -= self.Y0[ 1 : : 9]
            self.Y0[ 3 : : 9] -= self.Y0[ 4 : : 9]
	    self.Y0[ 6 : : 9] -= self.Y0[ 7 : : 9]

            # R (RU, RVL, RVH)
            self.Y0[ 2 : : 9] = 0.
            self.Y0[ 5 : : 9] = 0.
	    self.Y0[ 8 : : 9] = 0.
	    print ("check!!!!"),  self.Y0[ 1 : : 9],  self.Y0[ 4 : : 9],  self.Y0[ 7 : : 9]

        else:
	    
            SU, IU, RU, SVL,  IVL, RVL, SVH, IVH, RVH = self.getLastValues()
            
            self.Y0[ 0 : : 9] = (1 - self.parameters.proportionVaccinatedLow - self.parameters.proportionVaccinatedHigh)* SU
            self.Y0[ 3 : : 9] = SVL + self.parameters.proportionVaccinatedLow * SU
	    self.Y0[ 6 : : 9] = SVH + self.parameters.proportionVaccinatedHigh * SU
            self.Y0[ 1 : : 9] = IU
            self.Y0[ 4 : : 9] = IVL
	    self.Y0[ 7 : : 9] = IVH
            self.Y0[ 2 : : 9] = RU
            self.Y0[ 5 : : 9] = RVL
	    self.Y0[ 8 : : 9] = RVH

    def RHS(self, Y, t):
        '''
        SEIR model with multiple host types.
        
        This function gives the right-hand sides of the ODEs.
        '''
        
        # Convert vector to meaningful component vectors

	SU = Y[ 0 : : 9]
        IU = Y[ 1 : : 9]
        RU = Y[ 2 : : 9]
        SVL = Y[ 3 : : 9]
        IVL = Y[ 4 : : 9]
        RVL = Y[ 5 : : 9]
	SVH = Y[ 6 : : 9]
        IVH = Y[ 7 : : 9]
        RVH = Y[ 8: : 9]
        
        N = sum(SU + IU + RU + SVL + IVL + RVL+ SVH + IVH + RVH)
      
        # The force of infection
        Lambda = self.parameters.transmissionScaling * self.parameters.susceptibility \
                 * numpy.dot(self.parameters.contactMatrix,
                             self.parameters.transmissibility * (IU + IVL + IVH)) / N
        
        # The right-hand sides
        dSU = - Lambda * SU
        dIU = Lambda * SU - (self.parameters.recoveryRate + self.parameters.deathRateU) * IU
        dRU = self.parameters.recoveryRate * IU
        
	dSVL = - (1 - self.parameters.vaccineEfficacyVsInfectionLow) * Lambda * SVL
        dIVL = Lambda * SVL  - (self.parameters.recoveryRate + self.parameters.deathRateV) * IVL
        dRVL = self.parameters.recoveryRate * IVL

	dSVH = - (1 - self.parameters.vaccineEfficacyVsInfectionHigh) * Lambda * SVH
        dIVH = Lambda * SVH - (self.parameters.recoveryRate + self.parameters.deathRateV) * IVH
        dRVH = self.parameters.recoveryRate * IVH
        
        
        # Convert meaningful component vectors into a single vector
        dY = numpy.empty(Y.size, dtype = float)
        dY[ 0 : : 9] = dSU
        dY[ 1 : : 9] = dIU
        dY[ 2 : : 9] = dRU
        dY[ 3 : : 9] = dSVL
        dY[ 4 : : 9] = dIVL
        dY[ 5 : : 9] = dRVL
        dY[ 6 : : 9] = dSVH
        dY[ 7 : : 9] = dIVH
        dY[ 8 : : 9] = dRVH
      
        return dY
    
    def resetSolution(self):
        self.hasSolution = False

    def solve(self, tStart = 0., tEnd = None, tStep = 1.):
        if tEnd == None:
            tEnd = self.tMax

        if self.hasSolution:
            TOld  = self.T.copy()
            SUOld = self.SU.copy()
            IUOld = self.IU.copy()
            RUOld = self.RU.copy()
            SVLOld = self.SVL.copy()
            IVLOld = self.IVL.copy()
            RVLOld = self.RVL.copy()
            SVHOld = self.SVH.copy()
            IVHOld = self.IVH.copy()
            RVHOld = self.RVH.copy()
            
            
        # Time vector for solution
        self.T = numpy.hstack((numpy.arange(tStart, tEnd, tStep), tEnd))
        
        # Integrate the ODE
        from scipy.integrate import odeint
        self.Y = odeint(self.RHS,
                        self.Y0.copy(),
                        self.T,
                        mxstep = 1000)
        Z = self.Y.copy()
	self.SU = Z[:,  0 : : 9]
        self.IU = Z[:,  1 : : 9]
        self.RU = Z[:,  2 : : 9]
        self.SVL = Z[:,  3 : : 9]
        self.IVL = Z[:,  4 : : 9]
        self.RVL = Z[:,  5 : : 9]
        self.SVH = Z[:,  6 : : 9]
        self.IVH = Z[:,  7 : : 9]
        self.RVH = Z[:,  8 : : 9]

        if self.hasSolution:
	   
            self.T = numpy.hstack((TOld, self.T))

            self.SU = numpy.vstack((SUOld, self.SU))
            self.IU = numpy.vstack((IUOld, self.IU))
            self.RU = numpy.vstack((RUOld, self.RU))
            self.SVL = numpy.vstack((SVLOld, self.SVL))
            self.IVL = numpy.vstack((IVLOld, self.IVL))
            self.RVL = numpy.vstack((RVLOld, self.RVL))
            self.SVH = numpy.vstack((SVHOld, self.SVH))
            self.IVH = numpy.vstack((IVHOld, self.IVH))
            self.RVH = numpy.vstack((RVHOld, self.RVH))

        self.hasSolution = True

    def updateStats(self):
        self.NU = self.SU +  self.IU + self.RU
	self.NVL = self.SVL +  self.IVL + self.RVL
	self.NVH = self.SVH +  self.IVH + self.RVH
        self.NV = self.NVL + self.NVH
        self.N  = self.NU + self.NV

        self.infectionsU = self.NU[0, :] - self.SU[-1, :]
        self.infectionsVL = self.NVL[0, :] - self.SVL[-1, :]
	self.infectionsVH = self.NVH[0, :] - self.SVH[-1, :]

        # Find duplicate times: these are where vaccination occurs
        for i in numpy.arange(len(self.T)).compress(numpy.diff(self.T) == 0):
            # Update for vaccinated people
            self.infectionsU += self.SU[i + 1, :] - self.SU[i, :]
            self.infectionsVL += self.SVL[i + 1, :] - self.SVL[i, :]
	    self.infectionsVH += self.SVH[i + 1, :] - self.SVH[i, :]

	
	self.infectionsV  = self.infectionsVL + self.infectionsVH        
	self.infections  = self.infectionsU + self.infectionsV
        self.totalInfections = self.infections.sum()
        
	
	# Hospitalization of *vaccinated* risk individuals
        self.hospitalizationV = self.infectionsV* self.parameters.caseHospitalization * (1 - self.parameters.vaccineEfficacyVsHospitalization)
	##Hospitalization rate of unvaccinated individuals
        self.hospitalizationU = self.infectionsU* self.parameters.caseHospitalization 
	##Hospitalization rate of vaccinated individuals
       
        self.hospitalizations = self.hospitalizationV + self.hospitalizationU
        self.totalHospitalizations = self.hospitalizations.sum()
        
        self.deathsU = self.NU[0, :] - self.NU[-1, :]
        self.deathsVL = self.NVL[0, :] - self.NVL[-1, :]
	self.deathsVH = self.NVH[0, :] - self.NVH[-1, :]

        # Find duplicate times: these are where vaccination occurs
        for i in numpy.arange(len(self.T)).compress(numpy.diff(self.T) == 0):
            # Update for vaccinated people
            self.deathsU += self.NU[i + 1, :] - self.NU[i, :]
            self.deathsVL += self.NVL[i + 1, :] - self.NVL[i, :]
	    self.deathsVH += self.NVH[i + 1, :] - self.NVH[i, :]

        self.deaths  = self.deathsU + self.deathsVL + self.deathsVH
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

   
	
        self.parameters.proportionVaccinatedLowPW.values = PVPWVal[0]
	self.parameters.proportionVaccinatedHighPW.values = PVPWVal[1]

	## extend to full ages groups. Proportions calculated by multiplying PVPWVal 
	##values with the matrix defined in S.130
	
	self.parameters.proportionVaccinatedLow = self.parameters.proportionVaccinatedLowPW.full(self.parameters.ages)
	self.parameters.proportionVaccinatedHigh = self.parameters.proportionVaccinatedHighPW.full(self.parameters.ages)
	
	
	
	if self.hasSolution:
	   
            vacsUsedLow = (self.parameters.proportionVaccinatedLow * IC[0])
	    vacsUsedHigh = (self.parameters.proportionVaccinatedHigh * IC[0])

        else:
	    vacsUsedLow = (self.parameters.proportionVaccinatedLow
                        * self.parameters.population)

            vacsUsedHigh = (self.parameters.proportionVaccinatedHigh
                        * self.parameters.population)
	
	   
        # Update initial condition for ODEs
        self.updateIC()

        return vacsUsedLow, vacsUsedHigh
    
    def vaccinated_output(self):
        return list(self.parameters.proportionVaccinatedLow), list(self.parameters.proportionVaccinatedHigh), [(a*b) for (a,b) in zip(self.parameters.proportionVaccinatedLow, self.parameters.population)],[(a*b) for (a,b) in zip(self.parameters.proportionVaccinatedHigh, self.parameters.population)]
    
    
        
    def simulateWithVaccine(self, PVPWVals, vacEfficacy):
	
        nVacTypes = len(vacEfficacy)

        self.resetSolution()

        # Vaccinate the population
	vacsUsedLow, vacsUsedHigh = self.updateProportionVaccinated(PVPWVals, nVacTypes)
	
	if min(list(vacsUsedHigh)) <0 and min(list(PVPWVals)) >0 : print ("check!!!!"), PVPWVals, nVacTypes, vacsUsedHigh
        tEnd = self.tMax
	tStart = self.tMin
	self.solve(tStart = tStart, tEnd = tEnd)
	
	print self.IU.sum(axis=1)
	import matplotlib.pyplot as plt
	times = [num for num in xrange(tEnd+1)]
	plt.plot(times, self.IU.sum(axis=1), color = "red")
	plt.plot(times, self.IVL.sum(axis=1), color = "blue")
	plt.plot(times, self.IVH.sum(axis=1), color = "green")
		 
	plt.show()


        self.updateStats()

	return vacsUsedLow, vacsUsedHigh, self.parameters.proportionVaccinatedLow, self.parameters.proportionVaccinatedHigh 

    def outputSolution(self):
        for (i, t) in enumerate(self.T):
            print '%g' % t,
            print '\t'.join(map(lambda f: ('%g' % f),
                                 sum(self.IU[i, :]),
                                 sum(self.RU[i, :]),
                                 sum(self.SV[i, :]),
                                 sum(self.IV[i, :]),
                                 sum(self.RV[i, :])))
    
    def outputInfo(self):
        print 'R0:\t\t\t %g' % self.parameters.R0
        print ('Infections (in millions):'),  self.totalInfections/1000000.
        print 'Deaths:\t\t\t %g' % self.totalDeaths
        print 'Hospitalizations:\t %g' % self.totalHospitalizations
        #print ('Age-specific infections:'), list(self.infections)
	print ("unvaccinated"),self.infectionsU.sum()
	print ("vaccinated"),self.infectionsV.sum()
	print ("total pop size"), self.parameters.population.sum()
	print ("total unvaccinated"),  ((1 - self.parameters.proportionVaccinatedLow -  self.parameters.proportionVaccinatedHigh) * self.parameters.population).sum()
	print ("total vaccinated Low"),  ((self.parameters.proportionVaccinatedLow) * self.parameters.population).sum()
	print ("total vaccinated High"),  ((self.parameters.proportionVaccinatedHigh) * self.parameters.population).sum()


    def optimization_output(self):
	return self.parameters.proportionVaccinatedLow ,  self.parameters.proportionVaccinatedHigh


    def short_output(self):
	return list(self.infections), list(self.hospitalizations), list(self.deaths), list(self.DALY)

    def debug_info(self):
	return self.infectionsU.sum()
	

    #def vaccinated_output(self):
    #    return self.parameters.population.sum(),((1 - self.parameters.proportionVaccinated) * self.parameters.population).sum(), ((self.parameters.proportionVaccinated) * self.parameters.population).sum(), self.infectionsU.sum(), self.infectionsV.sum()


