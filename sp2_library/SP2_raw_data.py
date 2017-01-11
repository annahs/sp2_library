#module for dealing with SP2 raw data
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
from SP2_particle_record_UTC import ParticleRecord

def make_plot(record):
	x_vals_all = record.getAcqPoints()
	y_vals_all = record.getScatteringSignal()	
	y_vals_split = record.getSplitDetectorSignal()
	y_vals_incand = record.getWidebandIncandSignal()

	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax1.plot(x_vals_all,y_vals_all,'o', markerfacecolor='None', label = 'HG scattering signal')  
	ax1.plot(x_vals_all, y_vals_incand, color ='red',marker = 'o', linestyle = 'None', label = 'HG incandescent signal')
	ax1.set_xlabel('data point #')
	ax1.set_ylabel('amplitude (a.u.)')
	#ax1.plot(x_vals_all, y_vals_split, 'o', color ='green')

	plt.legend()
	plt.show(block=False)
	raw_input("Hit Enter see next")
	plt.close()

	
def defineIncandInsertStatement(number_of_channels,instrument_locn):
	if number_of_channels == 8:
		insert_statement = ('''INSERT INTO  sp2_single_particle_data_locn''' + str(instrument_locn) + '''							  
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
		insert_statement = ('''INSERT INTO  sp2_single_particle_data_locn''' + str(instrument_locn) + '''							  
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

def defineNonincandInsertStatement(number_of_channels,instrument_locn):
	if number_of_channels == 8:
		insert_statement = ('''INSERT INTO  sp2_single_nonincand_particle_data_locn''' + str(instrument_locn) + '''							  
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
		insert_statement = ('''INSERT INTO  sp2_single_nonincand_particle_data_locn''' + str(instrument_locn) + '''							  
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

		if parameters['show_plot'] == True:
			print parameters['file_name'], record_index
			make_plot(particle_record)

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