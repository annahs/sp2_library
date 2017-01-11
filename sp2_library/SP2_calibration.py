#module for dealing with SP2 calibrations
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


def quadratic(x_vals, c0, c1, c2):
	return c0 + c1*x_vals + c2*x_vals*x_vals


def linear(x_vals, c0, c1):
	return c0 + c1*x_vals


def getCalibrationPoints(database_name, calibration_ID):

	#database connection
	cnx = mysql.connector.connect(user='root', password='', host='localhost', database=database_name)
	cursor = cnx.cursor()

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

	cnx.close()

	return calibration_pts


def calcrBCMass(mobility_diameter):
	Gysel_mass_parameters = [-9.60e-02,4.22e-03,-4.82e-05,6.26e-07,-5.65e-10,-1.99e-13,3.80e-16]  #from Gysel et al. Atmos. Meas. Tech., 4, 2851–2858, 2011  DOI:10.5194/amt-4-2851-2011
	rBC_mass = Gysel_mass_parameters[0] + Gysel_mass_parameters[1]*mobility_diameter + Gysel_mass_parameters[2]*mobility_diameter**2 + Gysel_mass_parameters[3]*mobility_diameter**3 + Gysel_mass_parameters[4]*mobility_diameter**4 + Gysel_mass_parameters[5]*mobility_diameter**5 + Gysel_mass_parameters[6]*mobility_diameter**6

	return rBC_mass


def getLists(calibration_pts):
	rBC_masses = []
	incand_pk_hts = []
	for point in calibration_pts:
		mobility_diameter = point[0]
		incand_pk_ht = point[1]
		rBC_mass = calcrBCMass(mobility_diameter)
		
		incand_pk_hts.append(incand_pk_ht)
		rBC_masses.append(rBC_mass)

	return incand_pk_hts, rBC_masses


def fitLinear(calibration_pts):
	incand_pk_hts, rBC_masses = getLists(calibration_pts)

	xdata = np.array(incand_pk_hts)
	ydata = np.array(rBC_masses)

	popt, pcov = scipy.optimize.curve_fit(linear, xdata, ydata)				
	perr = np.sqrt(np.diag(pcov))

	residuals = ydata- linear(xdata, popt[0],popt[1])
	ss_res = np.sum(residuals**2)
	ss_tot = np.sum((ydata-np.mean(ydata))**2)
	r_squared = 1 - (ss_res / ss_tot)

	return popt, perr, r_squared


def fitQuadratic(calibration_pts):
	incand_pk_hts, rBC_masses = getLists(calibration_pts)

	xdata = np.array(incand_pk_hts)
	ydata = np.array(rBC_masses)

	popt, pcov = scipy.optimize.curve_fit(quadratic, xdata, ydata)				
	perr = np.sqrt(np.diag(pcov))
	r_squared = np.nan

	return popt, perr, r_squared


def plotCalib(calibration_pts, popt):

	incand_pk_hts, rBC_masses = getLists(calibration_pts)

	fit_rBC_masses = []
	if len(popt) == 2:
		for pk_ht in incand_pk_hts:
			fit_rBC_mass = linear(pk_ht,popt[0],popt[1])
			fit_rBC_masses.append(fit_rBC_mass)
	if len(popt) == 3:
		for pk_ht in incand_pk_hts:
			fit_rBC_mass = quadratic(pk_ht,popt[0],popt[1],popt[2])
			fit_rBC_masses.append(fit_rBC_mass)

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.scatter(incand_pk_hts,rBC_masses,color='r')
	ax.plot(incand_pk_hts,fit_rBC_masses, '-r')
	plt.xlabel('Incandescent peak height (a.u.)')
	plt.ylabel('rBC mass (fg)')

	ax.grid(True)
	plt.show()
	

def makeVariableDict(variable_list):
	variable_dict = {}
	i = 0
	for value in variable_list:
		variable_dict[i] = float(value)
		i+=1
	
	if 2 not in variable_dict:
		variable_dict[2] = None

	return variable_dict



def writeFitToDatabase(database_name,calibration_ID,popt,perr):

	popt_dict = makeVariableDict(popt)
	pcov_dict = makeVariableDict(perr)

	if len(popt) == 2:
		calibration_function = 'linear'
	if len(popt) == 3:
		calibration_function = 'quadratic'


	#database connection
	cnx = mysql.connector.connect(user='root', password='', host='localhost', database=database_name)
	cursor = cnx.cursor()

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

	cnx.close()	



