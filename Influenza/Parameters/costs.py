# Parameter values for the costs
from PiecewiseAgeParameter import PiecewiseAgeNumber, PiecewiseAgeRate
import demography

import numpy

#: Placeholder value.  Set later by fixing R0.
transmissibilityPW = PiecewiseAgeRate(
    [1.],
    [0])


#Wielders et al., 2010
caseARDSfractionPW =  PiecewiseAgeRate(
    [0.023 * 0.0721, 0.023 * 0.0259, 0.023 * 0.0766, 0.023 * 0.1217, 0.023 * 0.2873, 0.023 * 0.4163],
    [0, 15, 30, 45, 60, 75])


#Meier 2000
casePneumoniafractionPW =  PiecewiseAgeRate(
    [0.002, 0.003, 0.003, 0.01],
    [0, 15, 50, 65])

#Meier 2000
caseOtitisfractionPW =  PiecewiseAgeRate(
    [0.04 * 0.994, 0.04* 0.994, 0.007* 0.994, 0.003* 0.994, 0.02* 0.994],
    [0, 5, 15, 50, 65])


caseDeafnessfractionPW =  PiecewiseAgeRate(
    [0.04 * 0.006, 0.04* 0.006, 0.007* 0.006, 0.003* 0.006, 0.02* 0.006],
    [0, 5, 15, 50, 65])

#: +--------------------------------------------+-----------------------+
#: | Disability weight                          | Reference             |
#: +============================================+=======================+
#: | 0.006 (Mild)                               | Global burden of      |
#: +--------------------------------------------+                       |
#: | 0.051  (Moderate)                          | disease 2016          |
#: +--------------------------------------------+-----------------------|
#: | 0.18  (ARDS)                               |  Plass 2014           |
#: +--------------------------------------------+-----------------------|
#: | 0.09  (otitis moderate)                    |  IHME            |
#: +--------------------------------------------+                       |
#: | 0.18 (Ages 0-4), 0.17 (Ages 5+)            |
#: |Otitis Deafness                             |     Murray 1996 GBD   |
#: +--------------------------------------------+-----------------------|
#: | 0.279                                      |   Niessen 2009        |
#: +--------------------------------------------+-----------------------+

disabilityWeightUncomplicated = 0.051
disabilityWeightHospitalizedUncomplicated = 0.133
disabilityWeightARDS = 0.18
disabilityWeightPneumonia = 0.279
disabilityWeightOtitis = 0.09
disabilityWeightDeafnessPW = PiecewiseAgeRate( 
	[0.18, 0.17],
        [0,5])


#Niessen 2009
durationPneumonia = 0.0383562 # 2 weeks
durationOtitis = 0.0191781  # 1 week

#: For years of life lost. 
#: From US Life Tables 2014
expectationOfLifePW = PiecewiseAgeNumber(
    [78.9,78.3,74.4,69.5,64.5,59.7,54.9,50.2,45.4,40.7,36.1,31.7, 27.4,23.3,19.4,15.7,12.3,9.2,6.7,4.6,3.2,2.3],
    [0,1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100])

