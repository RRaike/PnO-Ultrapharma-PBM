# PnO-Ultrapharma-PBM
#
# This README.md will explain the function of each file in the PnO-Ultrapharma-PBM repository. The most important files will be discussed first.
#
################################################################################################################################################
# Most important files
#
# CSD/FinalFile.m : This file is the automation of finding a model that fits on the experimental data. It does this for all experimental data.
# CSD/BruteForceFunction.m: This is the function with which the FinalFile.m tries to find a minimum. In this case a brute force function.
# ComparisonTesting.m: Visualization of the different fits found in FinalFile.m and their respective experimental data. Also has some handpicked values.
#		These values are the values for every experimental dataset for which the RSME is minimal (for all models tried on it).
#
#################################################################################################################################################
# Underlying/ lesser files
#
#
# DataCharacteristics.m: Calculates for a given set all their respective D10, D50, D90 and span. Also gives number of elements in distribution.
#
#
# CSD/Model/main.m: Script that computes the PBM. Has all parameters in it.
# CSD/Model/CSD_model.m: Function version of Model/main.m , so it can be used in a higher script to call for a model.
# CSD/Model/Solving_PBE.m: Defines boundary conditions for the ode15s solver and calls it. 
#			Splits the ouput in their respective pieces: L (Size data), C (concentration data), T (Temperature Data). All as function of z.
# CSD/Model/Stelsel.m: Inside the ode15s solver. This is the discretized model as differential equations.
#
#
# CSD/Experimental/CSD_Experimental.m: Function that returns the experimental data. Also is capable of filtering unwanted datapoints, which must be
#				given to it. Gets data from the Excelfiles and not directly from the data gotten from Fiji.
# CSD/Experimental/Test_CSD.m: Script with same functionality as CSD_Experimental.m.
#
#
# CSD/HelperFunctions/CalculateRSME.m: Function that returns the Root Mean Square Error between the two input vectors.
# CSD/HelperFunctions/MakeUsableCSDDataForRSME.m: Function that makes a CSD model for input kinetic parameters. 
#						Then reduce the size by 1 (by using the AverageModelVector.m), so the RSME can be calculated.
#
#
###################################################################################################################################################
# Helper files
#
#
# DigFiles.m: Helper function that can go into folders to search for certain files.
# FindAllInFolder.m: Find everything in a certain folder. Can be more folders or files.
#
#
# CSD/AverageModelVector.m: Take the average value of every two succesive values in a vector. 
#			The length of the returned vector is one element shorter than the input vector.
# CSD/LoadExperimentalData.m: Find the Excel files that contain the experimental data.
#
#
# CSD/Experimental/EquilizeBins.m: Returns the bins of Fiji, but with equal bins. This is not automatic.
# CSD/Experimental/FilterExtendedSummary.m: Function that takes unwanted data out the set. It does not change the input data.
#					    The data used is the Excel files and not data directly from Fiji.
# CSD/Experimental/InputExcel.m: Function that reads the Excel files into MATLAB.
# CSD/Experimental/getSampleName.m: Function that makes the name of a sample. Input is the date, sample number and number in sample. This is useful to quickly
#			  	    find a certain picture once the pictures are named properly.
#
#
# CSD/HelperFunctions/LeftNRightValue.m: Given a value and the step size. Take the value left and right of that given value. Assumes equal steps always.
#
#
#
#


