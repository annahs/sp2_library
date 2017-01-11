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
import mysql.connector



def lognorm(x_vals, A, w, xc):
	return A/(np.sqrt(2*math.pi)*w*x_vals)*np.exp(-(np.log(x_vals/xc))**2/(2*w**2))

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


def getInstrInfo(instr_owner,instr_number,database_name):

	#database connection
	cnx = mysql.connector.connect(user='root', password='', host='localhost', database=database_name)
	cursor = cnx.cursor()


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
			AND instr_number = %s
		''',
		(instr_owner,instr_number))

	instr_info = cursor.fetchall()
	instr_id 				= instr_info[0][0]
	number_of_channels 		= instr_info[0][1]
	acquisition_rate		= instr_info[0][2]
	bytes_per_record 		= int(instr_info[0][3])
	min_detectable_signal 	= instr_info[0][4]
	
	cnx.close()
	return instr_id,number_of_channels,acquisition_rate,bytes_per_record,min_detectable_signal