from PiecewiseAgeParameter import PiecewiseAgeParameter, PiecewiseAgeRate
from ages import ages, vaccinationAgesTypical, vaccinationAgesUniversal
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

    def __init__(self, index, calibration, **paramValues):

	self.passedParamValues = ParamDict(paramValues)	

	self.ages = numpy.array(ages)
	self.vaccinationAgesTypical = numpy.array(vaccinationAgesTypical)
	self.vaccinationAgesUniversal = numpy.array(vaccinationAgesUniversal)


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


	self.population_highrisk = [(a*b) for (a,b) in zip(self.population, self.proportionHighRisk)]
        self.population_lowrisk = self.population - self.population_highrisk
	
	self.TypicalvaccineEfficacyVsInfectionTypical_H1 = self.passedParamValues["vacEfficacy"] [0] * self.relative_TypicalvaccineEfficacyVsInfection_H1
	self.TypicalvaccineEfficacyVsInfectionTypical_H3 = self.passedParamValues["vacEfficacy"] [0] * self.relative_TypicalvaccineEfficacyVsInfection_H3
	self.TypicalvaccineEfficacyVsInfectionTypical_B = self.passedParamValues["vacEfficacy"] [0] * self.relative_TypicalvaccineEfficacyVsInfection_B
	
	self.UniversalvaccineEfficacyVsInfectionUniversal_H1 = self.passedParamValues["vacEfficacy"] [1] * self.relative_UniversalvaccineEfficacyVsInfection_H1
	self.UniversalvaccineEfficacyVsInfectionUniversal_H3 = self.passedParamValues["vacEfficacy"] [1] * self.relative_UniversalvaccineEfficacyVsInfection_H3
	self.UniversalvaccineEfficacyVsInfectionUniversal_B = self.passedParamValues["vacEfficacy"] [1] * self.relative_UniversalvaccineEfficacyVsInfection_B

        # Set up proportion vaccinated vectors
        self.proportionVaccinatedTypicalLPW = PiecewiseAgeRate([0.0] * len(vaccinationAgesTypical),
            vaccinationAgesTypical)
	self.proportionVaccinatedTypicalHPW = PiecewiseAgeRate([0.0] * len(vaccinationAgesTypical),
            vaccinationAgesTypical)
	self.proportionVaccinatedUniversalLPW = PiecewiseAgeRate([0.0] * len(vaccinationAgesUniversal),
            vaccinationAgesUniversal)
        self.proportionVaccinatedUniversalHPW = PiecewiseAgeRate([0.0] * len(vaccinationAgesUniversal),
            vaccinationAgesUniversal)
	
	if len(vaccinationAgesTypical) != len(vaccinationAgesUniversal):
		raise ValueError, "The number of age group bins for low and high efficacy vaccine should be the same!"
	self.proportionVaccinatedTLPW = PiecewiseAgeRate([0.0] * len(vaccinationAgesTypical),
            vaccinationAgesTypical)
	self.proportionVaccinatedNLPW = PiecewiseAgeRate([0.0] * len(vaccinationAgesUniversal),
            vaccinationAgesUniversal)
	self.proportionVaccinatedTHPW = PiecewiseAgeRate([0.0] * len(vaccinationAgesTypical),
            vaccinationAgesTypical)
	self.proportionVaccinatedNHPW = PiecewiseAgeRate([0.0] * len(vaccinationAgesUniversal),
            vaccinationAgesUniversal)
		
	
	
	
	
	self.proportionVaccinatedTL = self.proportionVaccinatedTypicalLPW.full(self.ages)
	self.proportionVaccinatedTH = self.proportionVaccinatedTypicalHPW.full(self.ages)
	self.proportionVaccinatedNL = self.proportionVaccinatedUniversalLPW.full(self.ages)
	self.proportionVaccinatedNH = self.proportionVaccinatedUniversalHPW.full(self.ages)
	self.proportionVaccinatedL = self.proportionVaccinatedTL + self.proportionVaccinatedNL
	self.proportionVaccinatedH = self.proportionVaccinatedTH + self.proportionVaccinatedNH
	self.proportionVaccinatedTypical = self.proportionVaccinatedTL + self.proportionVaccinatedTH
	self.proportionVaccinatedUniversal = self.proportionVaccinatedNL + self.proportionVaccinatedNH
	

	self.proportionVaccinatedLLength = len(vaccinationAgesTypical)
	self.proportionVaccinatedHLength = len(vaccinationAgesTypical)
        

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


        #self.transmissionScaling_H1 = 1.0
	#self.transmissionScaling_H1 *=  self.R0 / self.computeR0_H1()
	#self.transmissionScaling_H3 = 1.0
	#self.transmissionScaling_H3 *=  self.R0 / self.computeR0_H3()
	#self.transmissionScaling_B = 1.0
	#self.transmissionScaling_B *=  self.R0 / self.computeR0_B()
    
	if calibration:
	    if "betaList" in self.passedParamValues:
		self.transmissionScaling_H1 =  self.passedParamValues["betaList"][0]
		self.transmissionScaling_H3 =  self.passedParamValues["betaList"][1]
		self.transmissionScaling_B =   self.passedParamValues["betaList"][2]
	    
	    if "vac_eff_hospitalization" in self.passedParamValues: 
		self.vac_eff_hospitalization = self.passedParamValues["vac_eff_hospitalization"]
	    
	    if "vac_eff_mortality" in self.passedParamValues: 
		self.vac_eff_mortality = self.passedParamValues["vac_eff_mortality"]

	    
	
	#print ("scaling H1 =========="), self.transmissionScaling_H1
	#print ("scaling H3 =========="), self.transmissionScaling_H3
	#print ("scaling B =========="), self.transmissionScaling_B
	
    """    
    def computeR0_H1(self):
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
	

        FUL_H1 = self.transmissionScaling_H1  * numpy.outer(self.susceptibility_H1 * sUL0, self.transmissibility) * self.contactMatrix
	FUH_H1 = self.transmissionScaling_H1  * numpy.outer(self.susceptibility_H1 * sUH0, self.transmissibility) * self.contactMatrix
	FTL_H1 = self.transmissionScaling_H1  * numpy.outer((1 - self.vaccineEfficacyVsInfectionTypical_H1) * self.susceptibility_H1 * sTL0, self.transmissibility) * self.contactMatrix
	FTH_H1 = self.transmissionScaling_H1  * numpy.outer((1 - self.vaccineEfficacyVsInfectionTypical_H1) * self.susceptibility_H1* sTH0, self.transmissibility) * self.contactMatrix
	FNL_H1 = self.transmissionScaling_H1  * numpy.outer((1 - self.vaccineEfficacyVsInfectionUniversal_H1) * self.susceptibility_H1 * sNL0, self.transmissibility) * self.contactMatrix
	FNH_H1 = self.transmissionScaling_H1  * numpy.outer((1 - self.vaccineEfficacyVsInfectionUniversal_H1)* self.susceptibility_H1 * sNH0, self.transmissibility) * self.contactMatrix
        
	
	F = numpy.vstack((numpy.hstack((FUL_H1, FUL_H1, FUL_H1, FUL_H1, FUL_H1, FUL_H1)),
			  numpy.hstack((FUH_H1, FUH_H1, FUH_H1, FUH_H1, FUH_H1, FUH_H1)),
			  numpy.hstack((FTL_H1, FTL_H1, FTL_H1, FTL_H1, FTL_H1, FTL_H1)),
			  numpy.hstack((FTH_H1, FTH_H1, FTH_H1, FTH_H1, FTH_H1, FTH_H1)),
			  numpy.hstack((FNL_H1, FNL_H1, FNL_H1, FNL_H1, FNL_H1, FNL_H1)),
			  numpy.hstack((FNH_H1, FNH_H1, FNH_H1, FNH_H1, FNH_H1, FNH_H1))))
                    

	##death rate of typical and universal vaccine is assumed to be the same.
        V = numpy.diag(numpy.hstack(
            (self.recoveryRate + self.deathRateUL_H1,
	     self.recoveryRate + self.deathRateUH_H1,
	     self.recoveryRate + self.deathRateVL_H1,
	     self.recoveryRate + self.deathRateVH_H1,
             self.recoveryRate + self.deathRateVL_H1,
	     self.recoveryRate + self.deathRateVH_H1)))
	
        G = numpy.dot(F, numpy.linalg.inv(V))

        (Lambda, Nu) = numpy.linalg.eig(G)

        i = numpy.argmax(Lambda.real)
        LambdaMax = numpy.real_if_close(Lambda[i])
        
        assert numpy.isreal(LambdaMax), \
               "Complex maximal eigenvalue %s!" % LambdaMax

        return LambdaMax
    
    def computeR0_H3(self):
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
	

        FUL_H3 = self.transmissionScaling_H3  * numpy.outer(self.susceptibility_H3 * sUL0, self.transmissibility) * self.contactMatrix
	FUH_H3 = self.transmissionScaling_H3  * numpy.outer(self.susceptibility_H3 * sUH0, self.transmissibility) * self.contactMatrix
	FTL_H3 = self.transmissionScaling_H3  * numpy.outer((1 - self.vaccineEfficacyVsInfectionTypical_H3) * self.susceptibility_H3 * sTL0, self.transmissibility) * self.contactMatrix
	FTH_H3 = self.transmissionScaling_H3  * numpy.outer((1 - self.vaccineEfficacyVsInfectionTypical_H3) * self.susceptibility_H3* sTH0, self.transmissibility) * self.contactMatrix
	FNL_H3 = self.transmissionScaling_H3  * numpy.outer((1 - self.vaccineEfficacyVsInfectionUniversal_H3) * self.susceptibility_H3 * sNL0, self.transmissibility) * self.contactMatrix
	FNH_H3 = self.transmissionScaling_H3  * numpy.outer((1 - self.vaccineEfficacyVsInfectionUniversal_H3)* self.susceptibility_H3 * sNH0, self.transmissibility) * self.contactMatrix
        
	
	F = numpy.vstack((numpy.hstack((FUL_H3, FUL_H3, FUL_H3, FUL_H3, FUL_H3, FUL_H3)),
			  numpy.hstack((FUH_H3, FUH_H3, FUH_H3, FUH_H3, FUH_H3, FUH_H3)),
			  numpy.hstack((FTL_H3, FTL_H3, FTL_H3, FTL_H3, FTL_H3, FTL_H3)),
			  numpy.hstack((FTH_H3, FTH_H3, FTH_H3, FTH_H3, FTH_H3, FTH_H3)),
			  numpy.hstack((FNL_H3, FNL_H3, FNL_H3, FNL_H3, FNL_H3, FNL_H3)),
			  numpy.hstack((FNH_H3, FNH_H3, FNH_H3, FNH_H3, FNH_H3, FNH_H3))))
                    

	##death rate of typical and universal vaccine is assumed to be the same.
        V = numpy.diag(numpy.hstack(
            (self.recoveryRate + self.deathRateUL_H3,
	     self.recoveryRate + self.deathRateUH_H3,
	     self.recoveryRate + self.deathRateVL_H3,
	     self.recoveryRate + self.deathRateVH_H3,
             self.recoveryRate + self.deathRateVL_H3,
	     self.recoveryRate + self.deathRateVH_H3)))
	
        G = numpy.dot(F, numpy.linalg.inv(V))

        (Lambda, Nu) = numpy.linalg.eig(G)

        i = numpy.argmax(Lambda.real)
        LambdaMax = numpy.real_if_close(Lambda[i])
        
        assert numpy.isreal(LambdaMax), \
               "Complex maximal eigenvalue %s!" % LambdaMax

        return LambdaMax
    
    def computeR0_B(self):
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
	

        FUL_B = self.transmissionScaling_B  * numpy.outer(self.susceptibility_B * sUL0, self.transmissibility) * self.contactMatrix
	FUH_B = self.transmissionScaling_B  * numpy.outer(self.susceptibility_B * sUH0, self.transmissibility) * self.contactMatrix
	FTL_B = self.transmissionScaling_B  * numpy.outer((1 - self.vaccineEfficacyVsInfectionTypical_B) * self.susceptibility_B * sTL0, self.transmissibility) * self.contactMatrix
	FTH_B = self.transmissionScaling_B  * numpy.outer((1 - self.vaccineEfficacyVsInfectionTypical_B) * self.susceptibility_B* sTH0, self.transmissibility) * self.contactMatrix
	FNL_B = self.transmissionScaling_B  * numpy.outer((1 - self.vaccineEfficacyVsInfectionUniversal_B) * self.susceptibility_B * sNL0, self.transmissibility) * self.contactMatrix
	FNH_B = self.transmissionScaling_B  * numpy.outer((1 - self.vaccineEfficacyVsInfectionUniversal_B)* self.susceptibility_B * sNH0, self.transmissibility) * self.contactMatrix
        
	
	F = numpy.vstack((numpy.hstack((FUL_B, FUL_B, FUL_B, FUL_B, FUL_B, FUL_B)),
			  numpy.hstack((FUH_B, FUH_B, FUH_B, FUH_B, FUH_B, FUH_B)),
			  numpy.hstack((FTL_B, FTL_B, FTL_B, FTL_B, FTL_B, FTL_B)),
			  numpy.hstack((FTH_B, FTH_B, FTH_B, FTH_B, FTH_B, FTH_B)),
			  numpy.hstack((FNL_B, FNL_B, FNL_B, FNL_B, FNL_B, FNL_B)),
			  numpy.hstack((FNH_B, FNH_B, FNH_B, FNH_B, FNH_B, FNH_B))))
                    

	##death rate of typical and universal vaccine is assumed to be the same.
        V = numpy.diag(numpy.hstack(
            (self.recoveryRate + self.deathRateUL_B,
	     self.recoveryRate + self.deathRateUH_B,
	     self.recoveryRate + self.deathRateVL_B,
	     self.recoveryRate + self.deathRateVH_B,
             self.recoveryRate + self.deathRateVL_B,
	     self.recoveryRate + self.deathRateVH_B)))
	
        G = numpy.dot(F, numpy.linalg.inv(V))

        (Lambda, Nu) = numpy.linalg.eig(G)

        i = numpy.argmax(Lambda.real)
        LambdaMax = numpy.real_if_close(Lambda[i])
        
        assert numpy.isreal(LambdaMax), \
               "Complex maximal eigenvalue %s!" % LambdaMax

        return LambdaMax
    """
