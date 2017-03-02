#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
from pprint import pprint
from datetime import datetime
from datetime import timedelta
import mysql.connector
import math
import calendar
import os.path
from scipy.optimize import curve_fit
import inspect
from dateutil import parser


"""
this module is a collection of widely used generic methods 
"""

def lognorm(x_vals, A, mu, sig):
	return A/(np.sqrt(2*math.pi)*mu*x_vals)*np.exp(-(np.log(x_vals/sig))**2/(2*mu**2))

def linear(x_vals, c0, c1):
	return c0 + c1*x_vals

def quadratic(x_vals, c0, c1, c2):
	return c0 + c1*x_vals + c2*x_vals*x_vals

def polynomial3(x_vals, c0, c1, c2, c3):
	return c0 + c1*x_vals + c2*x_vals*x_vals + c3*x_vals*x_vals*x_vals

def polynomial4(x_vals, c0, c1, c2, c3, c4):
	return c0 + c1*x_vals + c2*x_vals*x_vals + c3*x_vals*x_vals*x_vals+ c4*x_vals*x_vals*x_vals*x_vals

def polynomial5(x_vals, c0, c1, c2, c3, c4, c5):
	return c0 + c1*x_vals + c2*x_vals*x_vals + c3*x_vals*x_vals*x_vals+ c4*x_vals*x_vals*x_vals*x_vals+ c5*x_vals*x_vals*x_vals*x_vals*x_vals

def expDecay(x_vals, c0, c1,c2):
	return c0 + c1*np.exp(-c2*x_vals)

def log(x_vals, c0, c1, c2):
	return c0/(1+c1*np.exp(-c2*x_vals))


def generateInitialValues(arg_number,initial_values):
	if initial_values == {}:
		p0 = ()
		i=0
		while i < arg_number:
			p0 = p0 + (1,)
			i+=1	
	else:
		p0 = initial_values['p0']

	return p0



def fitFunction(function,x_vals,y_vals,**initial_values):
	arg_specs = inspect.getargspec(function)
	arg_number =  len(arg_specs[0]) - 1

	#format the initial values for the fit, if provided
	p0_val = generateInitialValues(arg_number,initial_values)

	try:
		popt, pcov = curve_fit(function, np.array(x_vals), np.array(y_vals),p0=p0_val)	
	except Exception,e: 
		popt = []
		pcov = []
		i=0
		while i < arg_number:
			popt.append(np.nan)
			pcov.append(np.nan)
			i+=1	

	perr = np.sqrt(np.diag(pcov))

	return popt,perr


def calcGyselMass(mobility_diameter):
	Gysel_mass_parameters = [-9.60e-02,4.22e-03,-4.82e-05,6.26e-07,-5.65e-10,-1.99e-13,3.80e-16]  #from Gysel et al. Atmos. Meas. Tech., 4, 2851–2858, 2011 DOI is 10.5194/amt-4-2851-2011
	rBC_mass = Gysel_mass_parameters[0] + Gysel_mass_parameters[1]*mobility_diameter + Gysel_mass_parameters[2]*mobility_diameter**2 + Gysel_mass_parameters[3]*mobility_diameter**3 + Gysel_mass_parameters[4]*mobility_diameter**4 + Gysel_mass_parameters[5]*mobility_diameter**5 + Gysel_mass_parameters[6]*mobility_diameter**6

	return rBC_mass


def calculateVED(rBC_density,rBC_mass):  #fg to nm
	VED = (((rBC_mass/(10**15*rBC_density))*6/math.pi)**(1/3.0))*10**7
	return VED

def calculateMass(rBC_density,rBC_VED):  #nm to fg
	mass = ((rBC_VED/(10.**7))**3)*(math.pi/6.)*rBC_density*(10.**15)
	return mass

def valid_date(s):
	return parser.parse(s)


def getInstrID(instr_owner,instr_number,database_name):

	cnx = mysql.connector.connect(user='root', password='', host='localhost', database=database_name)
	cursor = cnx.cursor()

	cursor.execute('''
		SELECT 
			id
		FROM
			sp2_instrument_info 
		WHERE
			instr_owner = %s
			AND instr_number = %s
		''',
		(instr_owner,instr_number))

	instr_info = cursor.fetchall()
	instr_id   = instr_info[0][0]
	
	cnx.close()

	return instr_id


def getInstrInfo(instr_owner,instr_number,cursor):
	cursor.execute('''
		SELECT 
			id,
			number_of_channels,
			acquisition_rate,
			bytes_per_record,
			min_detectable_signal
		FROM
			sp2_instrument_info 
		WHERE
			instr_owner = %s
			AND id = %s
		''',
		(instr_owner,instr_number))

	instr_info = cursor.fetchall()
	instr_id 				= instr_info[0][0]
	number_of_channels 		= instr_info[0][1]
	acquisition_rate		= instr_info[0][2]
	bytes_per_record 		= int(instr_info[0][3])
	min_detectable_signal 	= instr_info[0][4]

	return instr_id,number_of_channels,acquisition_rate,bytes_per_record,min_detectable_signal



