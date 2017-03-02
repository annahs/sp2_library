#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from pprint import pprint
import math
import mysql.connector
from datetime import datetime
import calendar
import scipy.optimize
import SP2_utilities

"""
This module contains the methods needed for handling SP2 calibration data
"""


def getCalibrationPoints(calibration_ID,cnx,cursor):
	"""
	retrieve calibration points from database
	"""
	cursor.execute('''
		SELECT 
			mobility_diameter,
			incand_pk_ht
		FROM
			sp2_calibration_points  
		WHERE
			calibration_ID = %s
			AND id > %s
		''',
		(calibration_ID,0))

	calibration_pts = cursor.fetchall()

	return calibration_pts


def getPeakHeightAndMassLists(calibration_pts):
	rBC_masses = []
	incand_pk_hts = []
	for point in calibration_pts:
		mobility_diameter = point[0]
		incand_pk_ht = point[1]
		rBC_mass = SP2_utilities.calcGyselMass(mobility_diameter)
		
		incand_pk_hts.append(incand_pk_ht)
		rBC_masses.append(rBC_mass)

	return incand_pk_hts, rBC_masses



def makeVariableDict(variable_list):
	variable_dict = {}
	i = 0
	for value in variable_list:
		variable_dict[i] = float(value)
		i+=1
	
	if 2 not in variable_dict:
		variable_dict[2] = None

	return variable_dict



def writeFitToDatabase(calibration_ID,popt,perr,cnx,cursor):

	popt_dict = makeVariableDict(popt)
	pcov_dict = makeVariableDict(perr)

	if len(popt) == 2:
		calibration_function = 'linear'
	if len(popt) == 3:
		calibration_function = 'quadratic'


	cursor.execute('''
		UPDATE 
			sp2_calibrations
		SET
			0_term = %s,
			1_term = %s,
			2_term = %s,
			0_term_err = %s,
			1_term_err = %s,
			2_term_err = %s,
			calibration_function = %s
		WHERE
			id = %s
		''',
		(popt_dict[0],popt_dict[1],popt_dict[2],pcov_dict[0],pcov_dict[1],pcov_dict[2],calibration_function,calibration_ID))
	cnx.commit()

def calcRSquared(xdata,ydata,popt,fit_function):
	combined_data = zip(xdata,ydata)
	mean_y = np.mean(ydata)
	residuals2 = []
	diffs = []
	for value in combined_data:
		residual = value[1] - fit_function(value[0], popt[0],popt[1])
		residuals2.append(residual**2)
		diffs.append((value[1]-mean_y)**2)
	ss_res = np.sum(residuals2)
	ss_tot = np.sum(diffs)
	r_squared = 1 - (ss_res / ss_tot)
	return r_squared


