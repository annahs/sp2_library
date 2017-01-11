#module for dealing with SP2 configuration files
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
import re
from dateutil import parser



def getConfigData(config_file):
	
	with open(config_file, 'r') as f:

		#regular expressions for lines we want to parse
		time_stamp = r"Last\sDate\sUpdated\s?=\s?(.*)"
		particle_sample_rate = r"1\sof\sEvery\s?=\s?(\d+)"	 #this is write 1 out of every x particle records to file - this requires correction late on 
		time_sample_rate = r"Out\sof\s?=\s?(\d+)"          	 #this is write data 1 out of every x minutes - this we can filter by ignoring long intervals
		
		for line in f:
			match_time_stamp = re.search(time_stamp, line, flags=re.IGNORECASE)
			if match_time_stamp:
				ini_date = parser.parse(match_time_stamp.group(1))
				

			match_particle_sample_rate = re.search(particle_sample_rate, line, flags=re.IGNORECASE)
			if match_particle_sample_rate:
				sample_factor_particle = float(match_particle_sample_rate.group(1))
			

			match_time_sample_rate = re.search(time_sample_rate, line, flags=re.IGNORECASE)
			if match_time_sample_rate:
				sample_factor_time = float(match_time_sample_rate.group(1))


		return ini_date,sample_factor_particle,sample_factor_time



def _createConfigInsertStatement(config_table):
	#create insert statement variable for database updating
	add_interval = ('''INSERT INTO  ''' + config_table + '''
              (instr_ID,
              instr_locn_ID,
              UNIX_UTC_ts_int_start,
              UNIX_UTC_ts_int_end,
              1_in_x_particles,
              1_in_x_minutes)
              VALUES (
              %(instr_ID)s,
              %(instr_locn_ID)s,
              %(UNIX_UTC_ts_int_start)s,
              %(UNIX_UTC_ts_int_end)s,
              %(1_in_x_particles)s,
              %(1_in_x_minutes)s)''')

	return add_interval



def writeConfigData(parameters,UNIX_time_stamp_UTC_start,UNIX_time_stamp_UTC_end,sample_factor_particle,sample_factor_time):
	
	#open database connection
	cnx = mysql.connector.connect(user='root', password='', host='localhost', database=parameters['database'])
	cursor = cnx.cursor()

	add_interval = _createConfigInsertStatement(parameters['config_table'])

	single_record = {
					'instr_ID':parameters['instr_ID'],
              		'instr_locn_ID':parameters['instrument_locn'],
					'UNIX_UTC_ts_int_start': float(UNIX_time_stamp_UTC_start),
					'UNIX_UTC_ts_int_end': float(UNIX_time_stamp_UTC_end),
					'1_in_x_particles':  sample_factor_particle, 
					'1_in_x_minutes':    sample_factor_time
					}
					
	cursor.execute(add_interval, single_record)
	cnx.commit()

	cnx.close()


