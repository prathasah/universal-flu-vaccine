# -*- coding: iso-8859-1 -*-
#
# Epidemiological parameter values
#

from PiecewiseAgeParameter import PiecewiseAgeRate
import pandas as pd
df  = pd.read_csv("/Users/prathasah/Dropbox (Bansal Lab)/Git-files/universal-flu-vaccine/Influenza/Parameters/sampled_parameter_set.csv")


def recoveryRatePW(index):
    return PiecewiseAgeRate(
    [df.at[index, "recovery_rate"]],
    [0])


def proportionHighRiskPW(index):
    return PiecewiseAgeRate(
    [df.at[index, "proportionHighRisk_0"],
     df.at[index, "proportionHighRisk_2"],
     df.at[index, "proportionHighRisk_5"],
     df.at[index, "proportionHighRisk_19"],
     df.at[index, "proportionHighRisk_25"],
     df.at[index, "proportionHighRisk_50"],
     df.at[index, "proportionHighRisk_65"]],
    [0, 2, 5, 19, 25, 50, 65])

def susceptibility_H1PW(index):
    return PiecewiseAgeRate(
    [df.at[index, "susceptibility_H1_0"],
     df.at[index, "susceptibility_H1_4"],
     df.at[index, "susceptibility_H1_18"],
     df.at[index, "susceptibility_H1_65"]],
    [0, 4, 18, 65])

def susceptibility_H3PW(index):
    return PiecewiseAgeRate(
    [df.at[index, "susceptibility_H3_0"],
     df.at[index, "susceptibility_H3_4"],
     df.at[index, "susceptibility_H3_18"],
     df.at[index, "susceptibility_H3_65"]],
    [0, 4, 18, 65])

def susceptibility_BPW(index):
    return PiecewiseAgeRate(
    [df.at[index, "susceptibility_B_0"],
     df.at[index, "susceptibility_B_4"],
     df.at[index, "susceptibility_B_18"],
     df.at[index, "susceptibility_B_65"]],
    [0, 4, 18, 65])

def relative_vaccineEfficacyVsInfection_H1PW(index):
    return PiecewiseAgeRate(
    [df.at[index, "relative_vaccineEfficacyVsInfection_H1_0"],
     df.at[index, "relative_vaccineEfficacyVsInfection_H1_0.5"],
      df.at[index, "relative_vaccineEfficacyVsInfection_H1_5"],
     df.at[index, "relative_vaccineEfficacyVsInfection_H1_18"],
     df.at[index, "relative_vaccineEfficacyVsInfection_H1_50"]],
    [0, 0.5,  5, 18, 50])

def relative_vaccineEfficacyVsInfection_H3PW(index):
    return PiecewiseAgeRate(
    [df.at[index, "relative_vaccineEfficacyVsInfection_H3_0"],
     df.at[index, "relative_vaccineEfficacyVsInfection_H3_0.5"],
      df.at[index, "relative_vaccineEfficacyVsInfection_H3_5"],
     df.at[index, "relative_vaccineEfficacyVsInfection_H3_18"],
     df.at[index, "relative_vaccineEfficacyVsInfection_H3_50"]],
    [0, 0.5,  5, 18, 50])

def relative_vaccineEfficacyVsInfection_BPW(index):
    return PiecewiseAgeRate(
    [df.at[index, "relative_vaccineEfficacyVsInfection_B_0"],
     df.at[index, "relative_vaccineEfficacyVsInfection_B_0.5"],
      df.at[index, "relative_vaccineEfficacyVsInfection_B_5"],
     df.at[index, "relative_vaccineEfficacyVsInfection_B_18"],
     df.at[index, "relative_vaccineEfficacyVsInfection_B_50"]],
    [0, 0.5,  5, 18, 50])


def vaccineEfficacyVsHospitalization_H1PW(index):
    return PiecewiseAgeRate(
    [df.at[index, "vaccineEfficacyVsHospitalization_H1_0"],
     df.at[index, "vaccineEfficacyVsHospitalization_H1_0.5"],
     df.at[index, "vaccineEfficacyVsHospitalization_H1_16"],
     df.at[index, "vaccineEfficacyVsHospitalization_H1_65"]],
    [0, 0.5, 16, 65])

def vaccineEfficacyVsHospitalization_H3PW(index):
    return PiecewiseAgeRate(
    [df.at[index, "vaccineEfficacyVsHospitalization_H3_0"],
     df.at[index, "vaccineEfficacyVsHospitalization_H3_0.5"],
     df.at[index, "vaccineEfficacyVsHospitalization_H3_16"],
     df.at[index, "vaccineEfficacyVsHospitalization_H3_65"]],
    [0, 0.5, 16, 65])

def vaccineEfficacyVsHospitalization_BPW(index):
    return PiecewiseAgeRate(
    [df.at[index, "vaccineEfficacyVsHospitalization_B_0"],
     df.at[index, "vaccineEfficacyVsHospitalization_B_0.5"],
     df.at[index, "vaccineEfficacyVsHospitalization_B_16"],
     df.at[index, "vaccineEfficacyVsHospitalization_B_65"]],
    [0, 0.5, 16, 65])

def vaccineEfficacyVsDeath_H1PW(index):
    return PiecewiseAgeRate(
    [df.at[index, "vaccineEfficacyVsDeath_H1_0"],
     df.at[index, "vaccineEfficacyVsDeath_H1_0.5"],
     df.at[index, "vaccineEfficacyVsDeath_H1_18"],
     df.at[index, "vaccineEfficacyVsDeath_H1_65"]],
    [0, 0.5, 18,65])

def vaccineEfficacyVsDeath_H3PW(index):
    return PiecewiseAgeRate(
    [df.at[index, "vaccineEfficacyVsDeath_H3_0"],
     df.at[index, "vaccineEfficacyVsDeath_H3_0.5"],
     df.at[index, "vaccineEfficacyVsDeath_H3_18"],
     df.at[index, "vaccineEfficacyVsDeath_H3_65"]],
    [0, 0.5, 18,65])

def vaccineEfficacyVsDeath_BPW(index):
    return PiecewiseAgeRate(
    [df.at[index, "vaccineEfficacyVsDeath_B_0"],
     df.at[index, "vaccineEfficacyVsDeath_B_0.5"],
     df.at[index, "vaccineEfficacyVsDeath_B_18"],
     df.at[index, "vaccineEfficacyVsDeath_B_65"]],
    [0, 0.5, 18,65])

def highRiskvaccineEfficacyVsDeath_H1PW(index):
    return PiecewiseAgeRate(
    [df.at[index, "highRiskvaccineEfficacyVsDeath_H1_0"],
     df.at[index, "highRiskvaccineEfficacyVsDeath_H1_0.5"],
     df.at[index, "highRiskvaccineEfficacyVsDeath_H1_18"],
     df.at[index, "highRiskvaccineEfficacyVsDeath_H1_65"]],
    [0, 0.5, 18,65])

def highRiskvaccineEfficacyVsDeath_H3PW(index):
    return PiecewiseAgeRate(
    [df.at[index, "highRiskvaccineEfficacyVsDeath_H3_0"],
     df.at[index, "highRiskvaccineEfficacyVsDeath_H3_0.5"],
     df.at[index, "highRiskvaccineEfficacyVsDeath_H3_18"],
     df.at[index, "highRiskvaccineEfficacyVsDeath_H3_65"]],
    [0, 0.5, 18,65])

def highRiskvaccineEfficacyVsDeath_BPW(index):
    return PiecewiseAgeRate(
    [df.at[index, "highRiskvaccineEfficacyVsDeath_B_0"],
     df.at[index, "highRiskvaccineEfficacyVsDeath_B_0.5"],
     df.at[index, "highRiskvaccineEfficacyVsDeath_B_18"],
     df.at[index, "highRiskvaccineEfficacyVsDeath_B_65"]],
    [0, 0.5, 18,65])


def lowRiskcaseMortality_H1PW(index):
    return PiecewiseAgeRate(
    [df.at[index, "lowRiskcaseMortality_H1_0"],
     df.at[index, "lowRiskcaseMortality_H1_5"],
     df.at[index, "lowRiskcaseMortality_H1_18"],
     df.at[index, "lowRiskcaseMortality_H1_50"],
    df.at[index, "lowRiskcaseMortality_H1_65"],
    df.at[index, "lowRiskcaseMortality_H1_75"]],
    [0, 5, 18, 50, 65, 75])

def lowRiskcaseMortality_H3PW(index):
    return PiecewiseAgeRate(
    [df.at[index, "lowRiskcaseMortality_H3_0"],
     df.at[index, "lowRiskcaseMortality_H3_5"],
     df.at[index, "lowRiskcaseMortality_H3_18"],
     df.at[index, "lowRiskcaseMortality_H3_50"],
    df.at[index, "lowRiskcaseMortality_H3_65"],
    df.at[index, "lowRiskcaseMortality_H3_75"]],
    [0, 5, 18, 50, 65, 75])

def lowRiskcaseMortality_BPW(index):
    return PiecewiseAgeRate(
    [df.at[index, "lowRiskcaseMortality_B_0"],
     df.at[index, "lowRiskcaseMortality_B_5"],
     df.at[index, "lowRiskcaseMortality_B_18"],
     df.at[index, "lowRiskcaseMortality_B_50"],
    df.at[index, "lowRiskcaseMortality_B_65"],
    df.at[index, "lowRiskcaseMortality_B_75"]],
    [0, 5, 18, 50, 65, 75])

def highRiskcaseMortality_H1PW(index):
    return PiecewiseAgeRate(
    [df.at[index, "highRiskcaseMortality_H1_0"],
     df.at[index, "highRiskcaseMortality_H1_5"],
     df.at[index, "highRiskcaseMortality_H1_18"],
     df.at[index, "highRiskcaseMortality_H1_50"],
    df.at[index, "highRiskcaseMortality_H1_65"],
    df.at[index, "highRiskcaseMortality_H1_75"]],
    [0, 5, 18, 50, 65, 75])

def highRiskcaseMortality_H3PW(index):
    return PiecewiseAgeRate(
    [df.at[index, "highRiskcaseMortality_H3_0"],
     df.at[index, "highRiskcaseMortality_H3_5"],
     df.at[index, "highRiskcaseMortality_H3_18"],
     df.at[index, "highRiskcaseMortality_H3_50"],
    df.at[index, "highRiskcaseMortality_H3_65"],
    df.at[index, "highRiskcaseMortality_H3_75"]],
    [0, 5, 18, 50, 65, 75])

def highRiskcaseMortality_BPW(index):
    return PiecewiseAgeRate(
    [df.at[index, "highRiskcaseMortality_B_0"],
     df.at[index, "highRiskcaseMortality_B_5"],
     df.at[index, "highRiskcaseMortality_B_18"],
     df.at[index, "highRiskcaseMortality_B_50"],
    df.at[index, "highRiskcaseMortality_B_65"],
    df.at[index, "highRiskcaseMortality_B_75"]],
    [0, 5, 18, 50, 65, 75])


def lowRiskcaseHospitalization_H1PW(index):
    return PiecewiseAgeRate(
    [df.at[index, "lowRiskcaseHospitalization_H1_0"],
     df.at[index, "lowRiskcaseHospitalization_H1_5"],
     df.at[index, "lowRiskcaseHospitalization_H1_18"],
     df.at[index, "lowRiskcaseHospitalization_H1_50"],
    df.at[index, "lowRiskcaseHospitalization_H1_65"],
    df.at[index, "lowRiskcaseHospitalization_H1_75"]],
    [0, 5, 18, 50, 65, 75])

def lowRiskcaseHospitalization_H3PW(index):
    return PiecewiseAgeRate(
    [df.at[index, "lowRiskcaseHospitalization_H3_0"],
     df.at[index, "lowRiskcaseHospitalization_H3_5"],
     df.at[index, "lowRiskcaseHospitalization_H3_18"],
     df.at[index, "lowRiskcaseHospitalization_H3_50"],
    df.at[index, "lowRiskcaseHospitalization_H3_65"],
    df.at[index, "lowRiskcaseHospitalization_H3_75"]],
    [0, 5, 18, 50, 65, 75])

def lowRiskcaseHospitalization_BPW(index):
    return PiecewiseAgeRate(
    [df.at[index, "lowRiskcaseHospitalization_B_0"],
     df.at[index, "lowRiskcaseHospitalization_B_5"],
     df.at[index, "lowRiskcaseHospitalization_B_18"],
     df.at[index, "lowRiskcaseHospitalization_B_50"],
    df.at[index, "lowRiskcaseHospitalization_B_65"],
    df.at[index, "lowRiskcaseHospitalization_B_75"]],
    [0, 5, 18, 50, 65, 75])

def highRiskcaseHospitalization_H1PW(index):
    return PiecewiseAgeRate(
    [df.at[index, "highRiskcaseHospitalization_H1_0"],
     df.at[index, "highRiskcaseHospitalization_H1_5"],
     df.at[index, "highRiskcaseHospitalization_H1_18"],
     df.at[index, "highRiskcaseHospitalization_H1_50"],
    df.at[index, "highRiskcaseHospitalization_H1_65"],
    df.at[index, "highRiskcaseHospitalization_H1_75"]],
    [0, 5, 18, 50, 65, 75])

def highRiskcaseHospitalization_H3PW(index):
    return PiecewiseAgeRate(
    [df.at[index, "highRiskcaseHospitalization_H3_0"],
     df.at[index, "highRiskcaseHospitalization_H3_5"],
     df.at[index, "highRiskcaseHospitalization_H3_18"],
     df.at[index, "highRiskcaseHospitalization_H3_50"],
    df.at[index, "highRiskcaseHospitalization_H3_65"],
    df.at[index, "highRiskcaseHospitalization_H3_75"]],
    [0, 5, 18, 50, 65, 75])

def highRiskcaseHospitalization_BPW(index):
    return PiecewiseAgeRate(
    [df.at[index, "highRiskcaseHospitalization_B_0"],
     df.at[index, "highRiskcaseHospitalization_B_5"],
     df.at[index, "highRiskcaseHospitalization_B_18"],
     df.at[index, "highRiskcaseHospitalization_B_50"],
    df.at[index, "highRiskcaseHospitalization_B_65"],
    df.at[index, "highRiskcaseHospitalization_B_75"]],
    [0, 5, 18, 50, 65, 75])

def R0(index):
    return df.at[index, "R0"]


