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

	
	self.vaccineEfficacyVsInfectionTypical_H1 = self.passedParamValues["vacEfficacy"] [0] * self.relative_vaccineEfficacyVsInfection_H1
	self.vaccineEfficacyVsInfectionTypical_H3 = self.passedParamValues["vacEfficacy"] [0] * self.relative_vaccineEfficacyVsInfection_H3
	self.vaccineEfficacyVsInfectionTypical_B = self.passedParamValues["vacEfficacy"] [0] * self.relative_vaccineEfficacyVsInfection_B
	
	self.vaccineEfficacyVsInfectionUniversal_H1 = self.passedParamValues["vacEfficacy"] [1] * self.relative_vaccineEfficacyVsInfection_H1
	self.vaccineEfficacyVsInfectionUniversal_H3 = self.passedParamValues["vacEfficacy"] [1] * self.relative_vaccineEfficacyVsInfection_H3
	self.vaccineEfficacyVsInfectionUniversal_B = self.passedParamValues["vacEfficacy"] [1] * self.relative_vaccineEfficacyVsInfection_B

	# Death probability for low risk **vacccinated** individuals
        self.caseMortalityVL_H1 = self.lowRiskcaseMortality_H1 * (1 - self.vaccineEfficacyVsDeath_H1)
	self.caseMortalityVL_H3 = self.lowRiskcaseMortality_H3 * (1 - self.vaccineEfficacyVsDeath_H3)
	self.caseMortalityVL_B = self.lowRiskcaseMortality_B * (1 - self.vaccineEfficacyVsDeath_B)
	
	# Death probability for high risk **vacccinated** individuals
        self.caseMortalityVH_H1 = self.highRiskcaseMortality_H1 * (1 - self.highRiskvaccineEfficacyVsDeath_H1)
	self.caseMortalityVH_H3 = self.highRiskcaseMortality_H3 * (1 - self.highRiskvaccineEfficacyVsDeath_H3)
	self.caseMortalityVH_B = self.highRiskcaseMortality_B * (1 - self.highRiskvaccineEfficacyVsDeath_B)
	
	##equation S1.7. Death rate of low-risk unvaccinated individuals
        self.deathRateUL_H1 = self.recoveryRate  * self.lowRiskcaseMortality_H1/ (1 - self.lowRiskcaseMortality_H1)
	self.deathRateUL_H3 = self.recoveryRate  * self.lowRiskcaseMortality_H3/ (1 - self.lowRiskcaseMortality_H3)
	self.deathRateUL_B = self.recoveryRate  * self.lowRiskcaseMortality_B/ (1 - self.lowRiskcaseMortality_B)
	
	##Death rate of high-risk unvaccinated individuals
        self.deathRateUH_H1 = self.recoveryRate  * self.highRiskcaseMortality_H1 / (1 - self.highRiskcaseMortality_H1)
	self.deathRateUH_H3 = self.recoveryRate  * self.highRiskcaseMortality_H3 / (1 - self.highRiskcaseMortality_H3)
	self.deathRateUH_B = self.recoveryRate  * self.highRiskcaseMortality_B / (1 - self.highRiskcaseMortality_B)
	
	##Death rate of low-risk vaccinated individuals
        self.deathRateVL_H1 = self.recoveryRate  * self.caseMortalityVL_H1 / (1 - self.caseMortalityVL_H1)
	self.deathRateVL_H3 = self.recoveryRate  * self.caseMortalityVL_H3 / (1 - self.caseMortalityVL_H3)
	self.deathRateVL_B = self.recoveryRate  * self.caseMortalityVL_B / (1 - self.caseMortalityVL_B)
	

	
	##Death rate of high-risk vaccinated individuals
        self.deathRateVH_H1 = self.recoveryRate  * self.caseMortalityVH_H1  / (1 - self.caseMortalityVH_H1)
	self.deathRateVH_H3 = self.recoveryRate  * self.caseMortalityVH_H3  / (1 - self.caseMortalityVH_H3)
	self.deathRateVH_B = self.recoveryRate  * self.caseMortalityVH_B  / (1 - self.caseMortalityVH_B)
	
	
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
	self.proportionVaccinatedTypical = self.proportionVaccinatedTL + self.proportionVaccinatedTH
	self.proportionVaccinatedUniversal = self.proportionVaccinatedNL + self.proportionVaccinatedNH
	

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
	

        FUL_H1 = self.transmissionScaling  * numpy.outer(self.susceptibility_H1 * sUL0, self.transmissibility) * self.contactMatrix
	FUL_H3 = self.transmissionScaling  * numpy.outer(self.susceptibility_H3 * sUL0, self.transmissibility) * self.contactMatrix
	FUL_B = self.transmissionScaling  * numpy.outer(self.susceptibility_B * sUL0, self.transmissibility) * self.contactMatrix
	
	FUH_H1 = self.transmissionScaling  * numpy.outer(self.susceptibility_H1 * sUH0, self.transmissibility) * self.contactMatrix
	FUH_H3 = self.transmissionScaling  * numpy.outer(self.susceptibility_H3 * sUH0, self.transmissibility) * self.contactMatrix
	FUH_B = self.transmissionScaling  * numpy.outer(self.susceptibility_B * sUH0, self.transmissibility) * self.contactMatrix
	
	FTL_H1 = self.transmissionScaling  * numpy.outer((1 - self.vaccineEfficacyVsInfectionTypical_H1) * self.susceptibility_H1 * sTL0, self.transmissibility) * self.contactMatrix
	FTL_H3 = self.transmissionScaling  * numpy.outer((1 - self.vaccineEfficacyVsInfectionTypical_H3) * self.susceptibility_H3 * sTL0, self.transmissibility) * self.contactMatrix
	FTL_B = self.transmissionScaling  * numpy.outer((1 - self.vaccineEfficacyVsInfectionTypical_B) * self.susceptibility_B * sTL0, self.transmissibility) * self.contactMatrix
	
	FTH_H1 = self.transmissionScaling  * numpy.outer((1 - self.vaccineEfficacyVsInfectionTypical_H1) * self.susceptibility_H1* sTH0, self.transmissibility) * self.contactMatrix
	FTH_H3 = self.transmissionScaling  * numpy.outer((1 - self.vaccineEfficacyVsInfectionTypical_H3) * self.susceptibility_H3* sTH0, self.transmissibility) * self.contactMatrix
	FTH_B = self.transmissionScaling  * numpy.outer((1 - self.vaccineEfficacyVsInfectionTypical_B) * self.susceptibility_B* sTH0, self.transmissibility) * self.contactMatrix
	
	
	FNL_H1 = self.transmissionScaling  * numpy.outer((1 - self.vaccineEfficacyVsInfectionUniversal_H1) * self.susceptibility_H1 * sNL0, self.transmissibility) * self.contactMatrix
	FNL_H3 = self.transmissionScaling  * numpy.outer((1 - self.vaccineEfficacyVsInfectionUniversal_H3) * self.susceptibility_H3 * sNL0, self.transmissibility) * self.contactMatrix
	FNL_B = self.transmissionScaling  * numpy.outer((1 - self.vaccineEfficacyVsInfectionUniversal_B) * self.susceptibility_B * sNL0, self.transmissibility) * self.contactMatrix
	
	FNH_H1 = self.transmissionScaling  * numpy.outer((1 - self.vaccineEfficacyVsInfectionUniversal_H1)* self.susceptibility_H1 * sNH0, self.transmissibility) * self.contactMatrix
	FNH_H3 = self.transmissionScaling  * numpy.outer((1 - self.vaccineEfficacyVsInfectionUniversal_H3)* self.susceptibility_H3 * sNH0, self.transmissibility) * self.contactMatrix
	FNH_B = self.transmissionScaling  * numpy.outer((1 - self.vaccineEfficacyVsInfectionUniversal_B)* self.susceptibility_B * sNH0, self.transmissibility) * self.contactMatrix
        
	
	F = numpy.vstack((numpy.hstack((FUL_H1, FUL_H1, FUL_H1, FUL_H1, FUL_H1, FUL_H1,FUL_H1, FUL_H1, FUL_H1, FUL_H1, FUL_H1, FUL_H1,FUL_H1, FUL_H1, FUL_H1, FUL_H1, FUL_H1, FUL_H1)),
			  numpy.hstack((FUL_H3, FUL_H3, FUL_H3, FUL_H3, FUL_H3, FUL_H3,FUL_H3, FUL_H3, FUL_H3, FUL_H3, FUL_H3, FUL_H3,FUL_H3, FUL_H3, FUL_H3, FUL_H3, FUL_H3, FUL_H3)),
			  numpy.hstack((FUL_B, FUL_B, FUL_B, FUL_B, FUL_B, FUL_B, FUL_B, FUL_B, FUL_B, FUL_B, FUL_B, FUL_B, FUL_B, FUL_B, FUL_B, FUL_B, FUL_B, FUL_B)),
			  
			  numpy.hstack((FUH_H1, FUH_H1, FUH_H1, FUH_H1, FUH_H1, FUH_H1, FUH_H1, FUH_H1, FUH_H1, FUH_H1, FUH_H1, FUH_H1,FUH_H1, FUH_H1, FUH_H1, FUH_H1, FUH_H1, FUH_H1)),
			  numpy.hstack((FUH_H3, FUH_H3, FUH_H3, FUH_H3, FUH_H3, FUH_H3, FUH_H3, FUH_H3, FUH_H3, FUH_H3, FUH_H3, FUH_H3, FUH_H3, FUH_H3, FUH_H3, FUH_H3, FUH_H3, FUH_H3)),
			  numpy.hstack((FUH_B, FUH_B, FUH_B, FUH_B, FUH_B, FUH_B, FUH_B, FUH_B, FUH_B, FUH_B, FUH_B, FUH_B, FUH_B, FUH_B, FUH_B, FUH_B, FUH_B, FUH_B)),
			  
			  numpy.hstack((FTL_H1, FTL_H1, FTL_H1, FTL_H1, FTL_H1, FTL_H1, FTL_H1, FTL_H1, FTL_H1, FTL_H1, FTL_H1, FTL_H1, FTL_H1, FTL_H1, FTL_H1, FTL_H1, FTL_H1, FTL_H1)),
			  numpy.hstack((FTL_H3, FTL_H3, FTL_H3, FTL_H3, FTL_H3, FTL_H3, FTL_H3, FTL_H3, FTL_H3, FTL_H3, FTL_H3, FTL_H3, FTL_H3, FTL_H3, FTL_H3, FTL_H3, FTL_H3, FTL_H3)),
			  numpy.hstack((FTL_B, FTL_B, FTL_B, FTL_B, FTL_B, FTL_B, FTL_B, FTL_B, FTL_B, FTL_B, FTL_B, FTL_B, FTL_B, FTL_B, FTL_B, FTL_B, FTL_B, FTL_B)),
			  
			  numpy.hstack((FTH_H1, FTH_H1, FTH_H1, FTH_H1, FTH_H1, FTH_H1, FTH_H1, FTH_H1, FTH_H1, FTH_H1, FTH_H1, FTH_H1, FTH_H1, FTH_H1, FTH_H1, FTH_H1, FTH_H1, FTH_H1)),
			  numpy.hstack((FTH_H3, FTH_H3, FTH_H3, FTH_H3, FTH_H3, FTH_H3, FTH_H3, FTH_H3, FTH_H3, FTH_H3, FTH_H3, FTH_H3, FTH_H3, FTH_H3, FTH_H3, FTH_H3, FTH_H3, FTH_H3)),
			  numpy.hstack((FTH_B, FTH_B, FTH_B, FTH_B, FTH_B, FTH_B, FTH_B, FTH_B, FTH_B, FTH_B, FTH_B, FTH_B, FTH_B, FTH_B, FTH_B, FTH_B, FTH_B, FTH_B)),
			  
			  numpy.hstack((FNL_H1, FNL_H1, FNL_H1, FNL_H1, FNL_H1, FNL_H1, FNL_H1, FNL_H1, FNL_H1, FNL_H1, FNL_H1, FNL_H1, FNL_H1, FNL_H1, FNL_H1, FNL_H1, FNL_H1, FNL_H1)),
			  numpy.hstack((FNL_H3, FNL_H3, FNL_H3, FNL_H3, FNL_H3, FNL_H3, FNL_H3, FNL_H3, FNL_H3, FNL_H3, FNL_H3, FNL_H3, FNL_H3, FNL_H3, FNL_H3, FNL_H3, FNL_H3, FNL_H3)),
			  numpy.hstack((FNL_B, FNL_B, FNL_B, FNL_B, FNL_B, FNL_B, FNL_B, FNL_B, FNL_B, FNL_B, FNL_B, FNL_B, FNL_B, FNL_B, FNL_B, FNL_B, FNL_B, FNL_B)),
			  
			  numpy.hstack((FNH_H1, FNH_H1, FNH_H1, FNH_H1, FNH_H1, FNH_H1, FNH_H1, FNH_H1, FNH_H1, FNH_H1, FNH_H1, FNH_H1, FNH_H1, FNH_H1, FNH_H1, FNH_H1, FNH_H1, FNH_H1)),
			  numpy.hstack((FNH_H3, FNH_H3, FNH_H3, FNH_H3, FNH_H3, FNH_H3, FNH_H3, FNH_H3, FNH_H3, FNH_H3, FNH_H3, FNH_H3, FNH_H3, FNH_H3, FNH_H3, FNH_H3, FNH_H3, FNH_H3)),
			  numpy.hstack((FNH_B, FNH_B, FNH_B, FNH_B, FNH_B, FNH_B, FNH_B, FNH_B, FNH_B, FNH_B, FNH_B, FNH_B, FNH_B, FNH_B, FNH_B, FNH_B, FNH_B, FNH_B))))
                    

	##death rate of typical and universal vaccine is assumed to be the same.
        V = numpy.diag(numpy.hstack(
            (self.recoveryRate + self.deathRateUL_H1,
	     self.recoveryRate + self.deathRateUL_H3,
	     self.recoveryRate + self.deathRateUL_B,
	     
	     self.recoveryRate + self.deathRateUH_H1,
	     self.recoveryRate + self.deathRateUH_H3,
	     self.recoveryRate + self.deathRateUH_B,
	     
	     self.recoveryRate + self.deathRateVL_H1,
	     self.recoveryRate + self.deathRateVL_H3,
	     self.recoveryRate + self.deathRateVL_B,
	     
	     self.recoveryRate + self.deathRateVH_H1,
	     self.recoveryRate + self.deathRateVH_H3,
	     self.recoveryRate + self.deathRateVH_B,
	     
             self.recoveryRate + self.deathRateVL_H1,
	     self.recoveryRate + self.deathRateVL_H3,
	     self.recoveryRate + self.deathRateVL_B,
	     
	     self.recoveryRate + self.deathRateVH_H1,
	     self.recoveryRate + self.deathRateVH_H3,
	     self.recoveryRate + self.deathRateVH_B)))


        G = numpy.dot(F, numpy.linalg.inv(V))

        (Lambda, Nu) = numpy.linalg.eig(G)

        i = numpy.argmax(Lambda.real)
        LambdaMax = numpy.real_if_close(Lambda[i])
        
        assert numpy.isreal(LambdaMax), \
               "Complex maximal eigenvalue %s!" % LambdaMax

        return LambdaMax

