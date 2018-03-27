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

	
	self.vaccineEfficacyVsInfectionLow = self.passedParamValues["vacEfficacy"] [0] * self.relative_vaccineEfficacyVsInfection
	self.vaccineEfficacyVsInfectionHigh = self.passedParamValues["vacEfficacy"] [1] * self.relative_vaccineEfficacyVsInfection



        # Compute mortality rates
        # from case mortalities

	# Mortality of *vaccinated* risk individuals
        self.caseMortalityV = self.caseMortality \
                               * (1 - self.vaccineEfficacyVsDeath)
	
	##Death rate of unvaccinated individuals
        self.deathRateU = (self.recoveryRate * self.caseMortality) / (1 - self.caseMortality)
	##Death rate of vaccinated individuals
        self.deathRateV = (self.recoveryRate * self.caseMortalityV)/ (1 - self.caseMortalityV)
	
        # Set up proportion vaccinated vectors
        self.proportionVaccinatedLowPW = PiecewiseAgeRate([0.0] * len(vaccinationLowRiskAgesTypical),
            vaccinationLowRiskAgesTypical)
        self.proportionVaccinatedHighPW = PiecewiseAgeRate([0.0] * len(vaccinationLowRiskAgesUniversal),
            vaccinationLowRiskAgesUniversal)
	
	self.proportionVaccinatedLow = self.proportionVaccinatedLowPW.full(self.ages)
	self.proportionVaccinatedHigh = self.proportionVaccinatedHighPW.full(self.ages)
	if len(vaccinationLowRiskAgesTypical) != len(vaccinationLowRiskAgesUniversal):
		raise ValueError, "The number of age group bins for low and high efficacy vaccine should be the same!"
	self.proportionVaccinatedLength = len(vaccinationLowRiskAgesTypical)
        

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
	sU0 = s0 * (1 - self.proportionVaccinatedLow - self.proportionVaccinatedHigh)
	sVL0 = s0 * self.proportionVaccinatedLow
	sVH0 = s0 * self.proportionVaccinatedHigh
	

        FU = self.transmissionScaling \
              * numpy.outer(self.susceptibility * sU0,
                            self.transmissibility) * self.contactMatrix
	FVL = self.transmissionScaling \
              * numpy.outer((1 - self.vaccineEfficacyVsInfectionLow)
                            * self.susceptibility * sVL0,
                            self.transmissibility) * self.contactMatrix

	FVH = self.transmissionScaling \
              * numpy.outer((1 - self.vaccineEfficacyVsInfectionHigh)
                            * self.susceptibility * sVH0,
                            self.transmissibility) * self.contactMatrix
        
        F = numpy.vstack((numpy.hstack((FU, FU, FU)),
                          numpy.hstack((FVL, FVL, FU)),
			  numpy.hstack((FVH, FVH, FU))))

	##death rate of low efficacy and high efficacy vaccine is same
        V = numpy.diag(numpy.hstack(
            (self.recoveryRate + self.deathRateU,
             self.recoveryRate + self.deathRateV,
	     self.recoveryRate + self.deathRateV)))


        G = numpy.dot(F, numpy.linalg.inv(V))

        (Lambda, Nu) = numpy.linalg.eig(G)

        i = numpy.argmax(Lambda.real)
        LambdaMax = numpy.real_if_close(Lambda[i])
        
        assert numpy.isreal(LambdaMax), \
               "Complex maximal eigenvalue %s!" % LambdaMax

        return LambdaMax

