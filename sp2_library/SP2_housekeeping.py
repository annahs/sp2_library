#!/usr/bin/env python
# -*- coding: UTF-8 -*-
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
import SP2_utilities

"""
This module has methods for dealing with SP2 housekeeping files
"""

def defineHKInsertStatement(hk_table):
	#create insert statement variable for database updating
	add_interval = ('''INSERT INTO  ''' + hk_table + '''
              (instr_ID,
              instr_location_ID,
              UNIX_UTC_ts_int_start,
              UNIX_UTC_ts_int_end,
              sample_flow,
              yag_power,
              sheath_flow,
              yag_xtal_temp,
              chamber_temp,
              chamber_pressure)
              VALUES (
              %(instr_ID)s,
              %(instr_location_ID)s,
              %(UNIX_UTC_ts_int_start)s,
              %(UNIX_UTC_ts_int_end)s,
              %(sample_flow)s,
              %(yag_power)s,
              %(sheath_flow)s,
              %(yag_xtal_temp)s,
              %(chamber_temp)s,
              %(chamber_pressure)s)''')

	return add_interval


def HKfileToDatabase(hk_file,add_interval,parameters,last_ts,cnx,cursor):
	
	multiple_records = []
	prev_UNIX_time_stamp_UTC = last_ts
	i=1

	hk_file.readline()
	for line in hk_file: 
		newline = line.split()
		seconds_past_midnight = float(newline[parameters['seconds_past_midnight_col']])
		UNIX_date = calendar.timegm(parameters['file_date'].utctimetuple())
		UNIX_time_stamp_UTC_end = UNIX_date + seconds_past_midnight - parameters['timezone']*3600
		UNIX_time_stamp_UTC_start = prev_UNIX_time_stamp_UTC
		prev_UNIX_time_stamp_UTC = UNIX_time_stamp_UTC_end

		single_record = {
		'instr_ID':				parameters['instr_ID'],
		'instr_location_ID':	parameters['instr_location_ID'],
		'UNIX_UTC_ts_int_start':UNIX_time_stamp_UTC_start,
		'UNIX_UTC_ts_int_end': 	UNIX_time_stamp_UTC_end,
		'sample_flow':  		float(newline[parameters['sample_flow_col']  ]),
		'yag_power':    		float(newline[parameters['yag_power_col']    ]),
		'sheath_flow':  		float(newline[parameters['sheath_flow_col']  ]),
		'yag_xtal_temp':		float(newline[parameters['yag_xtal_temp_col']]),
		'chamber_temp':			float(newline[parameters['chamber_temp_col']]),
		'chamber_pressure':		float(newline[parameters['chamber_pressure_col']])*100  #mBar to Pa ,
		}

		#convert any NaNs to None for db
		for key in single_record:
			if math.isnan(single_record[key]):
				single_record[key] = None

		#add single records to multiple records list
		multiple_records.append((single_record))
		i+= 1

		#bulk insert to db table
		if i%1000 == 0:
			cursor.executemany(add_interval, multiple_records)
			cnx.commit()
			multiple_records = []

	try:
		last_ts = UNIX_time_stamp_UTC_end
	except:
		last_ts = np.nan

	return multiple_records, last_ts



def addHKKeysToRawDataTable(parameters,cnx,cursor):

	#select hk data from the hk table for a given time block
	cursor.execute('''
		SELECT 
			id, 
			UNIX_UTC_ts_int_start,
			UNIX_UTC_ts_int_end 
		FROM ''' + parameters['hk_table'] + ''' 
		WHERE 
			UNIX_UTC_ts_int_start >= %s 
			AND UNIX_UTC_ts_int_start < %s
			AND instr_ID = %s''',
		(parameters['UNIX_start'],parameters['UNIX_end'],parameters['instr_ID']))
	hk_data = cursor.fetchall()


	#loop through the hk data and add keys to the raw data table
	print len(hk_data), 'starting loop'
	i=0
	for data_point in hk_data:
		id = data_point[0]
		hk_start_time = data_point[1]
		hk_end_time = data_point[2] 
		
		cursor.execute(('''
			UPDATE ''' + parameters['raw_data_table']+ ''' 
			SET 
				hk_id = %s 
			WHERE 
				UNIX_UTC_ts_int_end >= %s 
				AND UNIX_UTC_ts_int_end < %s
				AND instr_ID = %s'''),
			(id,hk_start_time,hk_end_time,parameters['instr_ID']))	
		cnx.commit()
		
		i+=1
		if (i % 5000) == 0:
			print 'fraction done: ', round(i*1./len(hk_data),2)
			
	






	
		

