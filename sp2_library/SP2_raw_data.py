#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from pprint import pprint
from struct import *
import math
import mysql.connector
from datetime import datetime
import calendar
from SP2_particle_record import ParticleRecord

"""
This module contains methods for dealing with raw .sp2b files
"""


def make_plot(particle_record):
	"""
	plot raw data signals
	"""
	x_vals_all 			= particle_record.getAcqPoints()
	y_vals_scat_HG 		= particle_record.scatData	
	y_vals_scat_LG 		= particle_record.lowGainScatData
	y_vals_incand_HG 	= particle_record.wideBandIncandData
	y_vals_incand_LG 	= particle_record.lowGainWideBandIncandData
	y_vals_split_HG 	= particle_record.splitData

	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax1.plot(x_vals_all, y_vals_scat_HG,'o', markerfacecolor='None', label = 'HG scattering signal')  
	ax1.plot(x_vals_all, y_vals_incand_HG, color ='red',marker = 'o', linestyle = 'None', label = 'HG incandescent signal')
	
	if y_vals_scat_LG != []:
		ax1.plot(x_vals_all, y_vals_scat_LG,'s', markerfacecolor='None', label = 'LG scattering signal')  
	if y_vals_incand_LG != []:
		ax1.plot(x_vals_all, y_vals_incand_LG, color ='red',marker = 's', linestyle = 'None', label = 'LG incandescent signal')
	
	ax1.set_xlabel('data point #')
	ax1.set_ylabel('amplitude (a.u.)')
	#ax1.plot(x_vals_all, y_vals_split, 'o', color ='green')

	plt.legend()
	plt.show()


	
def defineIncandInsertStatement(number_of_channels,table_name):
	"""
	define the datbase insert statement for incandescent particles
	"""
	if number_of_channels == 8:
		insert_statement = ('''INSERT INTO  ''' + table_name + '''							  
				  (instr_ID,
				  instr_location_ID,
				  sp2b_file, 
				  file_index, 
				  UNIX_UTC_ts_int_start,
				  UNIX_UTC_ts_int_end,
				  BB_incand_HG_pkht,
				  BB_incand_HG_pkpos,
				  BB_incand_LG_pkht,
				  BB_incand_LG_pkpos,
				  NB_incand_HG_pkht,
				  NB_incand_LG_pkht)
				  VALUES (
				  %(instr_ID)s,
				  %(instr_location_ID)s,
				  %(sp2b_file)s, 
				  %(file_index)s, 
				  %(UNIX_UTC_ts_int_start)s,
				  %(UNIX_UTC_ts_int_end)s,
				  %(BB_incand_HG_pkht)s,
				  %(BB_incand_HG_pkpos)s,
				  %(BB_incand_LG_pkht)s,
				  %(BB_incand_LG_pkpos)s,
				  %(NB_incand_HG_pkht)s,
				  %(NB_incand_LG_pkht)s
				  )''')

	if number_of_channels == 4:
		insert_statement = ('''INSERT INTO  ''' + table_name + '''							  
				  (instr_ID,
				  instr_location_ID,
				  sp2b_file, 
				  file_index, 
				  UNIX_UTC_ts_int_start,
				  UNIX_UTC_ts_int_end,
				  BB_incand_HG_pkht,
				  BB_incand_HG_pkpos,
				  NB_incand_HG_pkht)
				  VALUES (
				  %(instr_ID)s,
				  %(instr_location_ID)s,
				  %(sp2b_file)s, 
				  %(file_index)s, 
				  %(UNIX_UTC_ts_int_start)s,
				  %(UNIX_UTC_ts_int_end)s,
				  %(BB_incand_HG_pkht)s,
				  %(BB_incand_HG_pkpos)s,
				  %(NB_incand_HG_pkht)s
				  )''')

	return insert_statement

def defineNonincandInsertStatement(number_of_channels,table_name):
	"""
	define the database insert statement for non-incandescent particles
	"""
	if number_of_channels == 8:
		insert_statement = ('''INSERT INTO  ''' + table_name + '''							  
				  (instr_ID,
				  instr_location_ID,
				  sp2b_file, 
				  file_index, 
				  UNIX_UTC_ts_int_start,
				  UNIX_UTC_ts_int_end,
				  BB_scat_HG_pkht,
				  BB_scat_HG_pkpos,
				  BB_scat_LG_pkht,
				  BB_scat_LG_pkpos)
				  VALUES (
				  %(instr_ID)s,
				  %(instr_location_ID)s,
				  %(sp2b_file)s, 
				  %(file_index)s, 
				  %(UNIX_UTC_ts_int_start)s,
				  %(UNIX_UTC_ts_int_end)s,
				  %(BB_scat_HG_pkht)s,
				  %(BB_scat_HG_pkpos)s,
				  %(BB_scat_LG_pkht)s,
				  %(BB_scat_LG_pkpos)s
				  )''')

	if number_of_channels == 4:
		insert_statement = ('''INSERT INTO  ''' + table_name + '''							  
				  (instr_ID,
				  instr_location_ID,
				  sp2b_file, 
				  file_index, 
				  UNIX_UTC_ts_int_start,
				  UNIX_UTC_ts_int_end,
				  BB_scat_HG_pkht,
				  BB_scat_HG_pkpos)
				  VALUES (
				  %(instr_ID)s,
				  %(instr_location_ID)s,
				  %(sp2b_file)s, 
				  %(file_index)s, 
				  %(UNIX_UTC_ts_int_start)s,
				  %(UNIX_UTC_ts_int_end)s,
				  %(BB_scat_HG_pkht)s,
				  %(BB_scat_HG_pkpos)s
				  )''')

	return insert_statement


def writeIncandParticleData(sp2b_file,parameters,prev_particle_ts,count,insert_statement,cnx,cursor):
	"""
	Parse the raw data records and write information from incandescent particles to the database
	"""
	record_index = 0      
	multiple_records = []
	i=1
	while record_index < parameters['number_of_records']:
		#read the binary for a particle
		record = sp2b_file.read(parameters['bytes_per_record'])	
		particle_record = ParticleRecord(record, parameters['acq_rate'])

		#run the broadband high gain incandPeakInfo method to retrieve the BBHG incandescence peak attributes, this will be used to determine if this is an incandescent particle	
		particle_record.incandPeakInfo() 			
		bbhg_incand_pk_amp = float(particle_record.incandMax)
		bbhg_incand_pk_pos = float(particle_record.incandMaxPos)

		#if this is an incandescent particle that we can detect with the HG channel then continue 
		if bbhg_incand_pk_amp >= parameters['min_detectable_signal']:
			event_time = particle_record.timestamp  #UTC
			
			#get the remaining high gain signal information
			particle_record.narrowIncandPeakInfo()
			nbhg_incand_pk_amp = float(particle_record.narrowIncandMax)
			
			#get low gain signals if this is an 8-channel instrument
			if parameters['number_of_channels'] == 8:
				particle_record.incandPeakInfoLG()
				bblg_incand_pk_amp = float(particle_record.incandMax_LG)
				bblg_incand_pk_pos = float(particle_record.incandMaxPos_LG)
				particle_record.narrowIncandPeakInfoLG()
				nblg_incand_pk_amp = float(particle_record.narrowIncandMax_LG)

				single_record ={
				'instr_ID':parameters['instr_id'],
				'instr_location_ID':parameters['instr_locn_ID'],
				'sp2b_file':parameters['file_name'], 
				'file_index':record_index, 
				'UNIX_UTC_ts_int_start':prev_particle_ts,
				'UNIX_UTC_ts_int_end':event_time,
				'BB_incand_HG_pkht':bbhg_incand_pk_amp,
				'BB_incand_HG_pkpos':bbhg_incand_pk_pos,
				'BB_incand_LG_pkht':bblg_incand_pk_amp,
				'BB_incand_LG_pkpos':bblg_incand_pk_pos,
				'NB_incand_HG_pkht':nbhg_incand_pk_amp,
				'NB_incand_LG_pkht':nblg_incand_pk_amp,	
				}

			
			else:
				single_record ={
				'instr_ID':parameters['instr_id'],
				'instr_location_ID':parameters['instr_locn_ID'],
				'sp2b_file':parameters['file_name'], 
				'file_index':record_index, 
				'UNIX_UTC_ts_int_start':prev_particle_ts,
				'UNIX_UTC_ts_int_end':event_time,
				'BB_incand_HG_pkht':bbhg_incand_pk_amp,
				'BB_incand_HG_pkpos':bbhg_incand_pk_pos,
				'NB_incand_HG_pkht':nbhg_incand_pk_amp,
				}

			if 'mob_dia' in parameters:
				single_record['mob_dia'] = parameters['mob_dia']

			multiple_records.append((single_record))
			#bulk insert to db table
			if i%2000 == 0:
				cursor.executemany(insert_statement, multiple_records)
				cnx.commit()
				multiple_records = []


			prev_particle_ts = event_time

			#increment count of detectable incandescent particles
			i+= 1
			
		record_index+=1   

	#bulk insert of remaining records to db
	if multiple_records != []:
		cursor.executemany(insert_statement, multiple_records)
		cnx.commit()
		multiple_records = []
	count+=i

	return prev_particle_ts,count


def writeNonincandParticleData(sp2b_file,parameters,prev_particle_ts,count,insert_statement,cnx,cursor):
	"""
	Parse the raw data records and write information from non-incandescent particles to the database
	"""	
	record_index = 0      
	multiple_records = []
	i=1
	while record_index < parameters['number_of_records']:
		#read the binary for a particle
		record = sp2b_file.read(parameters['bytes_per_record'])	
		particle_record = ParticleRecord(record, parameters['acq_rate'])

		#run the broadband high gain incandPeakInfo method to retrieve the BBHG incandescence peak attributes, this will be used to determine if this is an incandescent particle	
		particle_record.incandPeakInfo() 
		particle_record.scatteringPeakInfo()			
		bbhg_incand_pk_amp = float(particle_record.incandMax)
		bbhg_scat_pk_amp = float(particle_record.scatteringMax)
		bbhg_scat_pk_pos = float(particle_record.scatteringMaxPos)

		#if this is an incandescent particle that we can detect with the HG channel then continue 
		if bbhg_incand_pk_amp < parameters['min_detectable_incand_signal'] and bbhg_scat_pk_amp > parameters['min_detectable_scat_signal']:
			event_time = particle_record.timestamp  #UTC

			#get low gain signals if this is an 8-channel instrument
			if parameters['number_of_channels'] == 8:
				particle_record.scatteringPeakInfoLG()
				bblg_scat_pk_amp = float(particle_record.scatteringMax_LG)
				bblg_scat_pk_pos = float(particle_record.scatteringMaxPos_LG)

				single_record ={
				'instr_ID':parameters['instr_id'],
				'instr_location_ID':parameters['instr_locn_ID'],
				'sp2b_file':parameters['file_name'], 
				'file_index':record_index, 
				'UNIX_UTC_ts_int_start':prev_particle_ts,
				'UNIX_UTC_ts_int_end':event_time,
				'BB_scat_HG_pkht':bbhg_scat_pk_amp,
				'BB_scat_HG_pkpos':bbhg_scat_pk_pos,
				'BB_scat_LG_pkht':bblg_scat_pk_amp,
				'BB_scat_LG_pkpos':bblg_scat_pk_pos,
				}

			
			else:
				single_record ={
				'instr_ID':parameters['instr_id'],
				'instr_location_ID':parameters['instr_locn_ID'],
				'sp2b_file':parameters['file_name'], 
				'file_index':record_index, 
				'UNIX_UTC_ts_int_start':prev_particle_ts,
				'UNIX_UTC_ts_int_end':event_time,
				'BB_scat_HG_pkht':bbhg_scat_pk_amp,
				'BB_scat_HG_pkpos':bbhg_scat_pk_pos,
				}


			multiple_records.append((single_record))
			#bulk insert to db table
			if i%2000 == 0:
				cursor.executemany(insert_statement, multiple_records)
				cnx.commit()
				multiple_records = []


			prev_particle_ts = event_time

			#increment count of detectable incandescent particles
			i+= 1
			
		record_index+=1   

	#bulk insert of remaining records to db
	if multiple_records != []:
		cursor.executemany(insert_statement, multiple_records)
		cnx.commit()
		multiple_records = []
	count+=i

	return prev_particle_ts,count	
