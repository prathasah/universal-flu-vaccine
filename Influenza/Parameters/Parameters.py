from PiecewiseAgeParameter import PiecewiseAgeParameter, PiecewiseAgeRate
from ages import ages, vaccinationLowRiskAgesTypical, vaccinationLowRiskAgesUniversal, vaccinationHighRiskAgesTypical, vaccinationHighRiskAgesUniversal
import demography
import epidemiology
import costs
import os
import numpy as np
import types
#from .. import fileIO


import numpy

import UserDict

class ParamDict(UserDict.UserDict):
    def valueOrAttrFromOther(self, key, other):
        '''
        Return value if key is in paramValues dict or return
        default as attribute from object other.
        '''
	return getattr(other, key)


    def epi_valueOrAttrFromOther(self, key, other, index):
        '''
        Return value if key is in paramValues dict or return
        default as attribute from object other.
        '''
	return (getattr(other, key))(index)
    
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
	self.setPWAttr(namePW,self.passedParamValues.valueOrAttrFromOther(namePW, other))
	

    def setAttrFromPassedOrOther(self, other, name):
        '''
        Set an attribute as self.name,
        taking value from passed paramValues
        or from attributes of object other.
        '''
	setattr(self, name, 
                self.passedParamValues.valueOrAttrFromOther(name, other))
	
	
	
    def epi_setPWAttr(self, namePW, value):
        '''
        Set a piecewise attribute as self.namePW
        and its expanded value as self.name.
        '''

        assert isinstance(value, PiecewiseAgeParameter)
        assert namePW.endswith('PW')
        name = namePW[ : -2]
        setattr(self, namePW, value)
        setattr(self, name, value.full(self.ages))
    
    def epi_setPWAttrFromPassedOrOther(self, other, namePW, index):
        '''
        Set a piecewise attribute as self.namePW
        and its expanded value as self.name,
        taking values from passed paramValues
        or from attributes of object other.
        '''
	self.epi_setPWAttr(namePW,self.passedParamValues.epi_valueOrAttrFromOther(namePW, other, index))
	

    def epi_setAttrFromPassedOrOther(self, other, name, index):
        '''
        Set an attribute as self.name,
        taking value from passed paramValues
        or from attributes of object other.
        '''
	setattr(self, name, 
                self.passedParamValues.epi_valueOrAttrFromOther(name, other, index))

##################################################################

    def __init__(self, index, **paramValues):

	self.passedParamValues = ParamDict(paramValues)	

	self.ages = numpy.array(ages)


        # Load in parameters and expand as necessary
	# Go through each files
        for m in (demography, costs):
	    # list all the modules
            for p in dir(m):
		#if module returns a numbers, then..
                if isinstance(getattr(m, p),(float, int)):
		    self.setAttrFromPassedOrOther(m, p)
		##if it is an agespecific parameter then..
                elif isinstance(getattr(m, p),
                                PiecewiseAgeParameter): 
		    self.setPWAttrFromPassedOrOther(m, p)
		    
		    
        for p in dir(epidemiology):
	    #if module returns a numbers, then..
	    func = getattr(epidemiology, p)
	    if isinstance(func,types.FunctionType):  
		if isinstance(func(index),(float, int)):
		    self.epi_setAttrFromPassedOrOther(epidemiology,p,index)
		##if it is an agespecific parameter then..
		elif isinstance(func(index),
				    PiecewiseAgeParameter): 
			self.epi_setPWAttrFromPassedOrOther(epidemiology,p, index)

	
	self.vaccineEfficacyVsInfectionTypical = self.passedParamValues["vacEfficacy"] [0] * self.relative_vaccineEfficacyVsInfection
	self.vaccineEfficacyVsInfectionUniversal = self.passedParamValues["vacEfficacy"] [1] * self.relative_vaccineEfficacyVsInfection



        # Compute mortality rates
        # from case mortalities

	## Death probability for low risk **unvacccinated** individuals
        self.caseMortalityL = self.caseMortality/ ((1 - self.proportionHighRisk) + (self.HighRiskRelativeCaseMortality * self.proportionHighRisk))
	
	## Death probability for high risk **unvacccinated** individuals
        self.caseMortalityH = self.caseMortalityL * self.HighRiskRelativeCaseMortality
	
	# Death probability for low risk **vacccinated** individuals
        self.caseMortalityVL = self.caseMortalityL * (1 - self.vaccineEfficacyVsDeath)
	
	# Death probability for high risk **vacccinated** individuals
        self.caseMortalityVH = self.caseMortalityVL * self.HighRiskRelativeCaseMortality
	
	##equation S1.7. Death rate of low-risk unvaccinated individuals
        self.deathRateUL = self.recoveryRate  * self.caseMortalityL/ (1 - self.caseMortalityL)
	##Death rate of low-risk vaccinated individuals
        self.deathRateVL = self.recoveryRate  * self.caseMortalityVL / (1 - self.caseMortalityVL)
	##Death rate of high-risk unvaccinated individuals
        self.deathRateUH = self.recoveryRate  * self.caseMortalityH / (1 - self.caseMortalityH)
	##Death rate of high-risk vaccinated individuals
        self.deathRateVH = self.recoveryRate  * self.caseMortalityVH  / (1 - self.caseMortalityVH)
	
	# Compute specific case hospitalizations
        self.caseHospitalizationL = self.caseHospitalization / ((1 - self.proportionHighRisk) + (self.HighRiskRelativeCaseHospitalization * self.proportionHighRisk))
        self.caseHospitalizationH = self.HighRiskRelativeCaseHospitalization * self.caseHospitalizationL
	
	
        # Set up proportion vaccinated vectors
        self.proportionVaccinatedTypicalLPW = PiecewiseAgeRate([0.0] * len(vaccinationLowRiskAgesTypical),
            vaccinationLowRiskAgesTypical)
	self.proportionVaccinatedTypicalHPW = PiecewiseAgeRate([0.0] * len(vaccinationHighRiskAgesTypical),
            vaccinationHighRiskAgesTypical)
	self.proportionVaccinatedUniversalLPW = PiecewiseAgeRate([0.0] * len(vaccinationLowRiskAgesUniversal),
            vaccinationLowRiskAgesUniversal)
        self.proportionVaccinatedUniversalHPW = PiecewiseAgeRate([0.0] * len(vaccinationHighRiskAgesUniversal),
            vaccinationHighRiskAgesUniversal)
	
	if len(vaccinationLowRiskAgesTypical) != len(vaccinationLowRiskAgesUniversal):
		raise ValueError, "The number of age group bins for low and high efficacy vaccine should be the same!"
	if len(vaccinationHighRiskAgesTypical) != len(vaccinationHighRiskAgesUniversal):
		raise ValueError, "The number of age group bins for low and high efficacy vaccine should be the same!"
	self.proportionVaccinatedTLPW = PiecewiseAgeRate([0.0] * len(vaccinationLowRiskAgesTypical),
            vaccinationLowRiskAgesTypical)
	self.proportionVaccinatedNLPW = PiecewiseAgeRate([0.0] * len(vaccinationLowRiskAgesUniversal),
            vaccinationLowRiskAgesUniversal)
	self.proportionVaccinatedTHPW = PiecewiseAgeRate([0.0] * len(vaccinationHighRiskAgesTypical),
            vaccinationHighRiskAgesTypical)
	self.proportionVaccinatedNHPW = PiecewiseAgeRate([0.0] * len(vaccinationHighRiskAgesUniversal),
            vaccinationHighRiskAgesUniversal)
		
	
	self.proportionVaccinatedTL = self.proportionVaccinatedTypicalLPW.full(self.ages)
	self.proportionVaccinatedTH = self.proportionVaccinatedTypicalHPW.full(self.ages)
	self.proportionVaccinatedNL = self.proportionVaccinatedUniversalLPW.full(self.ages)
	self.proportionVaccinatedNH = self.proportionVaccinatedUniversalHPW.full(self.ages)
	self.proportionVaccinatedL = self.proportionVaccinatedTL + self.proportionVaccinatedNL
	self.proportionVaccinatedH = self.proportionVaccinatedTH + self.proportionVaccinatedNH
	

	self.proportionVaccinatedLLength = len(vaccinationLowRiskAgesTypical)
	self.proportionVaccinatedHLength = len(vaccinationHighRiskAgesTypical)
        

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
	self.transmissionScaling *=  self.R0 / self.computeR0()
	
        
    def computeR0(self):
	#normalized population size for each age groups
        s0 = self.population / sum(self.population)
	sL0 = s0 * (1 - self.proportionHighRisk)
        sH0 = s0 * self.proportionHighRisk
	
	sUL0 = sL0 * (1 - self.proportionVaccinatedTL - self.proportionVaccinatedNL)
        sTL0 = sL0 * self.proportionVaccinatedTL
	sNL0 = sL0 * self.proportionVaccinatedNL
        sUH0 = sH0 * (1 - self.proportionVaccinatedTH - self.proportionVaccinatedNH)
        sTH0 = sH0 * self.proportionVaccinatedTH
	sNH0 = sH0 * self.proportionVaccinatedNH
	

        FUL = self.transmissionScaling  * numpy.outer(self.susceptibility * sUL0,
                            self.transmissibility) * self.contactMatrix
	FUH = self.transmissionScaling  * numpy.outer(self.susceptibility * sUH0,
                            self.transmissibility) * self.contactMatrix
	FTL = self.transmissionScaling  * numpy.outer((1 - self.vaccineEfficacyVsInfectionTypical)
                            * self.susceptibility * sTL0,
                            self.transmissibility) * self.contactMatrix
	FTH = self.transmissionScaling  * numpy.outer((1 - self.vaccineEfficacyVsInfectionTypical)
                            * self.susceptibility * sTH0,
                            self.transmissibility) * self.contactMatrix
	FNL = self.transmissionScaling  * numpy.outer((1 - self.vaccineEfficacyVsInfectionUniversal)
                            * self.susceptibility * sNL0,
                            self.transmissibility) * self.contactMatrix
	FNH = self.transmissionScaling  * numpy.outer((1 - self.vaccineEfficacyVsInfectionUniversal)
                            * self.susceptibility * sNH0,
                            self.transmissibility) * self.contactMatrix
        
        F = numpy.vstack((numpy.hstack((FUL, FUL, FUL, FUL, FUL, FUL)),
			  numpy.hstack((FUH, FUH, FUH, FUH, FUH, FUH)),
			  numpy.hstack((FTL, FTL, FTL, FTL, FTL, FTL)),
			  numpy.hstack((FTH, FTH, FTH, FTH, FTH, FTH)),
			  numpy.hstack((FNL, FNL, FNL, FNL, FNL, FNL)),
			  numpy.hstack((FNH, FNH, FNH, FNH, FNH, FNH))))
                    

	##death rate of typical and universal vaccine is assumed to be the same.
        V = numpy.diag(numpy.hstack(
            (self.recoveryRate + self.deathRateUL,
	     self.recoveryRate + self.deathRateUH,
	     self.recoveryRate + self.deathRateVL,
	     self.recoveryRate + self.deathRateVH,
             self.recoveryRate + self.deathRateVL,
	     self.recoveryRate + self.deathRateVH)))


        G = numpy.dot(F, numpy.linalg.inv(V))

        (Lambda, Nu) = numpy.linalg.eig(G)

        i = numpy.argmax(Lambda.real)
        LambdaMax = numpy.real_if_close(Lambda[i])
        
        assert numpy.isreal(LambdaMax), \
               "Complex maximal eigenvalue %s!" % LambdaMax

        return LambdaMax

