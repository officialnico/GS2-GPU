#! /usr/bin/env python3

import netCDF4 as nc
import numpy as np
import os

VERBOSE_LEVEL = 2 # 0: no extra prints, 1: shape prints, 2: data print


test_time_step_index = -2 # must match the one used in "compute_kinetic_energy_transfer_golden_answers.py"
filename = os.path.join(os.getcwd(), "test_kinetic_energy_transfer.out.nc")
data = nc.Dataset(filename)

all_tests_passed = True

###################################################################################
###################################################################################
##############################  test Tv_theta #####################################
###################################################################################
###################################################################################

if VERBOSE_LEVEL > 0:
    print("testing kinetic_energy_transfer_theta")

# get result from GS2 diagnostic
transfer = np.array(data["kinetic_energy_transfer_theta"][:])
A = transfer[test_time_step_index,:] # select timestep where we test

# get golden answer
B = np.loadtxt("kinetic_energy_transfer_theta_golden_answer.txt")

if VERBOSE_LEVEL > 0:
    print("     A.shape: ", A.shape)
    print("     B.shape: ", B.shape)
if VERBOSE_LEVEL > 1:
    print("         A: ", A)
    print("         B: ", B)

resultsAgree = np.allclose(A,B)

if VERBOSE_LEVEL > 0:
    if resultsAgree:
        print("#######################################")
        print("#######################################")
        print("############# TEST PASSED #############")
        print("#######################################")
        print("#######################################")
    else:
        print("#######################################")
        print("#######################################")
        print("############# TEST FAILED #############")
        print("#######################################")
        print("#######################################")

if not resultsAgree:
    all_tests_passed = False


###################################################################################
###################################################################################
########################## test Tv_theta_kxsource #################################
###################################################################################
###################################################################################

if VERBOSE_LEVEL > 0:
    print("testing kinetic_energy_transfer_theta_kxsource")

# get result from GS2 diagnostic
transfer = np.array(data["kinetic_energy_transfer_theta_kxsource"][:])
A = transfer[test_time_step_index,:] # select timestep where we test

# get golden answer
B = np.loadtxt("kinetic_energy_transfer_theta_kxsource_golden_answer.txt")

if VERBOSE_LEVEL > 0:
    print("     A.shape: ", A.shape)
    print("     B.shape: ", B.shape)
if VERBOSE_LEVEL > 1:
    print("         A: ", A)
    print("         B: ", B)

resultsAgree = np.allclose(A,B)

if VERBOSE_LEVEL > 0:
    if resultsAgree:
        print("#######################################")
        print("#######################################")
        print("############# TEST PASSED #############")
        print("#######################################")
        print("#######################################")
    else:
        print("#######################################")
        print("#######################################")
        print("############# TEST FAILED #############")
        print("#######################################")
        print("#######################################")

if not resultsAgree:
    all_tests_passed = False

###################################################################################
###################################################################################
########################## test Tv_theta_kxtarget #################################
###################################################################################
###################################################################################

if VERBOSE_LEVEL > 0:
    print("testing kinetic_energy_transfer_theta_kxtarget")

# get result from GS2 diagnostic
transfer = np.array(data["kinetic_energy_transfer_theta_kxtarget"][:])
A = transfer[test_time_step_index,:] # select timestep where we test

# get golden answer
B = np.loadtxt("kinetic_energy_transfer_theta_kxtarget_golden_answer.txt")

if VERBOSE_LEVEL > 0:
    print("     A.shape: ", A.shape)
    print("     B.shape: ", B.shape)
if VERBOSE_LEVEL > 1:
    print("         A: ", A)
    print("         B: ", B)

resultsAgree = np.allclose(A,B)

if VERBOSE_LEVEL > 0:
    if resultsAgree:
        print("#######################################")
        print("#######################################")
        print("############# TEST PASSED #############")
        print("#######################################")
        print("#######################################")
    else:
        print("#######################################")
        print("#######################################")
        print("############# TEST FAILED #############")
        print("#######################################")
        print("#######################################")

if not resultsAgree:
    all_tests_passed = False

###################################################################################
###################################################################################
########################## test Tv_theta_kysource #################################
###################################################################################
###################################################################################

if VERBOSE_LEVEL > 0:
    print("testing kinetic_energy_transfer_theta_kysource")

# get result from GS2 diagnostic
transfer = np.array(data["kinetic_energy_transfer_theta_kysource"][:])
A = transfer[-2,:] # select second to last timestep of gs2 output

# get golden answer
B = np.loadtxt("kinetic_energy_transfer_theta_kysource_golden_answer.txt")

if VERBOSE_LEVEL > 0:
    print("     A.shape: ", A.shape)
    print("     B.shape: ", B.shape)
if VERBOSE_LEVEL > 1:
    print("         A: ", A)
    print("         B: ", B)

resultsAgree = np.allclose(A,B)

if VERBOSE_LEVEL > 0:
    if resultsAgree:
        print("#######################################")
        print("#######################################")
        print("############# TEST PASSED #############")
        print("#######################################")
        print("#######################################")
        print()
        print()
    else:
        print("#######################################")
        print("#######################################")
        print("############# TEST FAILED #############")
        print("#######################################")
        print("#######################################")
        print()
        print()

if not resultsAgree:
    all_tests_passed = False

###################################################################################
###################################################################################
################################# total result ####################################
###################################################################################
###################################################################################


if resultsAgree:
    print("#######################################")
    print("#######################################")
    print("######## ALL TESTS PASSED #############")
    print("#######################################")
    print("#######################################")
else:
    print("#######################################")
    print("#######################################")
    print("###### AT LEAST 1 TEST FAILED #########")
    print("#######################################")
    print("#######################################")