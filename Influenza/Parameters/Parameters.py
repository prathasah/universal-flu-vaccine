from PiecewiseAgeParameter import PiecewiseAgeParameter, PiecewiseAgeRate
from ages import ages, vaccinationLowRiskAges, vaccinationHighRiskAges
import demography
import epidemiology
import costs
from .. import fileIO

import numpy

import UserDict

class ParamDict(UserDict.UserDict):
    def valueOrAttrFromOther(self, key, other):
        '''
        Return value if key is in paramValues dict or return
        default as attribute from object other.
        '''
        return self.get(key, getattr(other, key))

class Parameters:
    def setPWAttr(self, namePW, value):
        '''
        Set a piecewise attribute as self.namePW
        and its expanded value as self.name.
        '''
        assert isinstance(value, PiecewiseAgeParameter)
        assert namePW.endswith('PW')
        name = namePW[ : -2]
        setattr(self, namePW, value)
        setattr(self, name, value.full(self.ages))
    
    def setPWAttrFromPassedOrOther(self, other, namePW):
        '''
        Set a piecewise attribute as self.namePW
        and its expanded value as self.name,
        taking values from passed paramValues
        or from attributes of object other.
        '''
	
        self.setPWAttr(namePW,
                       self.passedParamValues.valueOrAttrFromOther(namePW,
                                                                   other))

    def setAttrFromPassedOrOther(self, other, name):
        '''
        Set an attribute as self.name,
        taking value from passed paramValues
        or from attributes of object other.
        '''
	
        setattr(self, name,
                self.passedParamValues.valueOrAttrFromOther(name, other))

    def __init__(self, **paramValues):
        self.passedParamValues = ParamDict(paramValues)
        self.ages = numpy.array(ages)

        # Load in parameters and expand as necessary
	# Go thrrough each files
        for m in (demography, epidemiology, costs):
	    # list all the modules
            for p in dir(m):
		#if module returns a numbers, then..
                if isinstance(getattr(m, p),(float, int)):
		    
		    self.setAttrFromPassedOrOther(m, p)
		##if it is an agespecific parameter then..
                elif isinstance(getattr(m, p),
                                PiecewiseAgeParameter):
		    self.setPWAttrFromPassedOrOther(m, p)

        # Compute mortality rates
        # from case mortalities

	##equatiion S1.5: Mortality of low risk individuals
        self.caseMortalityL = self.caseMortality \
                              / (1 +
                                 (self.highRiskRelativeCaseMortality - 1)
                                 * self.proportionHighRisk)
	##Mortality of high risk individuals
        self.caseMortalityH = self.caseMortalityL \
                              * self.highRiskRelativeCaseMortality
	# Mortality of low risk *vaccinated* risk individuals
        self.caseMortalityVL = self.caseMortalityL \
                               * (1 - self.vaccineEfficacyVsDeath)
	##Mortality of high risk *vaccinated* risk individuals
        self.caseMortalityVH = self.caseMortalityVL \
                               * self.highRiskRelativeCaseMortality
	
	##equation S1.7. Death rate of low-risk unvaccinated individuals
        self.deathRateUL = self.recoveryRate \
                           * self.caseMortalityL \
                           / (1 - self.caseMortalityL)
	##Death rate of low-risk vaccinated individuals
        self.deathRateVL = self.recoveryRate \
                           * self.caseMortalityVL \
                           / (1 - self.caseMortalityVL)
	##Death rate of high-risk unvaccinated individuals
        self.deathRateUH = self.recoveryRate \
                           * self.caseMortalityH \
                           / (1 - self.caseMortalityH)
	##Death rate of high-risk vaccinated individuals
        self.deathRateVH = self.recoveryRate \
                           * self.caseMortalityVH \
                           / (1 - self.caseMortalityVH)

        # Compute specific case hospitalizations
        self.caseHospitalizationL = self.caseHospitalization \
           / (1 + (self.highRiskRelativeCaseHospitalization - 1)
              * self.proportionHighRisk)
        self.caseHospitalizationH = self.highRiskRelativeCaseHospitalization \
                                    * self.caseHospitalizationL

        # Set up proportion vaccinated vectors
        self.proportionVaccinatedLPW = PiecewiseAgeRate(
            [0.0] * len(vaccinationLowRiskAges),
            vaccinationLowRiskAges)
        self.proportionVaccinatedHPW = PiecewiseAgeRate(
            [0.0] * len(vaccinationHighRiskAges),
            vaccinationHighRiskAges)
        self.proportionVaccinatedL = \
                                   self.proportionVaccinatedLPW.full(self.ages)
        self.proportionVaccinatedH = \
                                   self.proportionVaccinatedHPW.full(self.ages)
        self.proportionVaccinatedLLength = len(vaccinationLowRiskAges)
        self.proportionVaccinatedHLength = len(vaccinationHighRiskAges)
        self.proportionVaccinatedLength = self.proportionVaccinatedLLength \
                                          + self.proportionVaccinatedHLength

        # Get contact matrix
        if 'contactMatrix' in paramValues:
            self.contactMatrix = paramValues.get('contactMatrix')
        else:
            from sys import modules
            import os.path
            import cPickle
            modulePath = os.path.dirname(modules[self.__module__].__file__)
            contactMatrixFile = os.path.join(modulePath, 'contactMatrix.p')
            self.contactMatrix = cPickle.load(open(contactMatrixFile))

        # One last parameter to fit to R0
        self.transmissionScaling = 1.0
	#print ("Check1"), self.transmissionScaling, self.computeR0()
        self.transmissionScaling *=  self.R0 / self.computeR0()
	#print ("Check1"), self.transmissionScaling
        
    def computeR0(self):
        s0 = self.population / sum(self.population)
        sL0 = s0 * (1 - self.proportionHighRisk)
        sH0 = s0 * self.proportionHighRisk
        sUL0 = sL0 * (1 - self.proportionVaccinatedL)
        sVL0 = sL0 * self.proportionVaccinatedL
        sUH0 = sH0 * (1 - self.proportionVaccinatedH)
        sVH0 = sH0 * self.proportionVaccinatedH

        FUL = self.transmissionScaling \
              * numpy.outer(self.susceptibility * sUL0,
                            self.transmissibility) * self.contactMatrix
        FUH = self.transmissionScaling \
              * numpy.outer(self.susceptibility * sUH0,
                            self.transmissibility) * self.contactMatrix
        FVL = self.transmissionScaling \
              * numpy.outer((1 - self.vaccineEfficacyVsInfection)
                            * self.susceptibility * sVL0,
                            self.transmissibility) * self.contactMatrix
        FVH = self.transmissionScaling \
              * numpy.outer((1 - self.vaccineEfficacyVsInfection)
                            * self.susceptibility * sVH0,
                            self.transmissibility) * self.contactMatrix
        
        F = numpy.vstack((numpy.hstack((FUL, FUL, FUL, FUL)),
                          numpy.hstack((FUH, FUH, FUH, FUH)),
                          numpy.hstack((FVL, FVL, FVL, FVL)),
                          numpy.hstack((FVH, FVH, FVH, FVH))))

        V = numpy.diag(numpy.hstack(
            (self.recoveryRate + self.deathRateUL,
             self.recoveryRate + self.deathRateUH,
             self.recoveryRate + self.deathRateVL,
             self.recoveryRate + self.deathRateVH)))

        G = numpy.dot(F, numpy.linalg.inv(V))

        (Lambda, Nu) = numpy.linalg.eig(G)

        i = numpy.argmax(Lambda.real)
        LambdaMax = numpy.real_if_close(Lambda[i])
        
        assert numpy.isreal(LambdaMax), \
               "Complex maximal eigenvalue %s!" % LambdaMax

        return LambdaMax

    def getDumpData(self, runData = True):
        dumpData = fileIO.dumpContainer()

        #dumpData.ages = self.ages
        #dumpData.population = self.population
        #dumpData.proportionHighRisk = self.proportionHighRiskPW
        dumpData.R0 = self.R0
        dumpData.recoveryRate = self.recoveryRatePW
        dumpData.latencyRate = self.latencyRatePW
        dumpData.susceptibility = self.susceptibilityPW
        dumpData.transmissibility = self.transmissibilityPW
        dumpData.vaccineEfficacyVsInfection = \
            self.vaccineEfficacyVsInfectionPW
        dumpData.vaccineEfficacyVsDeath = \
            self.vaccineEfficacyVsDeathPW
        dumpData.caseMortality = self.caseMortalityPW
        dumpData.highRiskRelativeCaseMortality = \
            self.highRiskRelativeCaseMortalityPW
        dumpData.contactMatrix = self.contactMatrix
        dumpData.expectationOfLife = self.expectationOfLife
        dumpData.contingentValue = self.contingentValue
        dumpData.caseHospitalization = \
            self.caseHospitalizationPW
        dumpData.highRiskRelativeCaseHospitalization = \
            self.highRiskRelativeCaseHospitalizationPW

        return dumpData

    dump = fileIO.dump
