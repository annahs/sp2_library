#SP2 size distribution class
#timestamps are in UTC

import sys
import os
import numpy as np
from pprint import pprint
from datetime import datetime
from datetime import timedelta
import mysql.connector
import math
import matplotlib.pyplot as plt
import calendar
from scipy.optimize import curve_fit
from mysql_db_connection import dbConnection

class SizeDistribution(object,dbConnection):
		
	def __init__(self,database_name,instr_location_ID, instr_ID, interval_start, interval_end,minimum_bin,maximum_bin,bin_width):
		dbConnection.__init__(self, database_name)
		

		self.instr_location_ID 		= instr_location_ID 
		self.instr_ID     			= instr_ID
		self.interval_start 		= interval_start
		self.interval_end   		= interval_end
		self.binning_min    		= minimum_bin
		self.binning_max    		= maximum_bin
		self.binning_increment    	= bin_width
		self.temperature    		= 273.15  	#default is STP
		self.pressure       		= 101325  	#default is STP
		self.altitude				= 1 		#default is 1m above amsl
		self.rBC_density 			= 1.8 		#g/mol - Bond and Bergstrom 2006

		self.calibration_ID = None
		self.assembled_interval_data = None
		self.binned_data = None
		self.interval_db_id = None

		self._retrieveHousekeepingLimits()		#set values for discrimination based on housekeeping parameters
		self._retrieveInstrInfo()				#get basic instrument info
	
	def _retrieveHousekeepingLimits(self):
		self.db_cur.execute('''
		SELECT 
			yag_min,
			yag_max,
			sample_flow_min,
			sample_flow_max
		FROM
			sp2_hk_limits 
		WHERE
			instr_ID = %s
			AND location_id = %s
		LIMIT 1
		''',
		(self.instr_ID, self.instr_location_ID))
		hk_limits = self.db_cur.fetchall()

		self.yag_min = hk_limits[0][0]
		self.yag_max = hk_limits[0][1]
		self.sample_flow_min = hk_limits[0][2]
		self.sample_flow_max = hk_limits[0][3]


	def _retrieveInstrInfo(self):
		self.db_cur.execute('''
		SELECT 
			number_of_channels
		FROM
			sp2_instrument_info 
		WHERE
			id = %s
			AND id > %s
		LIMIT 1
		''',
		(self.instr_ID,0))
		
		instr_info = self.db_cur.fetchall()

		self.number_of_channels = instr_info[0][0]


	def retrieveSingleParticleData(self):
		self.db_cur.execute('''
		SELECT 
			sp.UNIX_UTC_ts_int_start,
			sp.UNIX_UTC_ts_int_end,
			sp.BB_incand_HG_pkht,
			sp.BB_incand_LG_pkht,
			hk.sample_flow
		FROM
			sp2_single_particle_data sp
				JOIN
			sp2_hk_data hk ON sp.HK_id = hk.id
		WHERE
			sp.UNIX_UTC_ts_int_end BETWEEN %s AND %s
			AND hk.yag_power BETWEEN %s AND %s
			AND hk.sample_flow BETWEEN %s AND %s
		''',
		(self.interval_start,self.interval_end,self.yag_min,self.yag_max,self.sample_flow_min,self.sample_flow_max))
		
		self.single_particle_data = self.db_cur.fetchall()
		
		

	def addIntervalDataToDatabase(self):
		
		#step 1: get the rBC core diameters and volume of air sampled in the interval
		self._assembleIntervalData()
		#step 2: bin the rBC core diameters
		self._binAssembledData()
		#step 3: write the interval data to the database
		self._writeIntervalDataToDatabase()
		#step 4: write the binned rBC core data to the database
		self._writeBinnedDataToDatabase()
		
	
	def _assembleIntervalData(self):
		interval_data_dict = {}
		interval_sampled_volume = 0
		interval_mass = 0
		interval_mass_uncertainty = 0
		ved_list = []
		i =0
		for row in self.single_particle_data:
			ind_start_time = row[0] #UNIX UTC timestamp
			ind_end_time = row[1]	#UNIX UTC timestamp
			BB_incand_HG = row[2]  	#in arbitrary units
			BB_incand_LG = row[3]  	#in arbitrary units
			sample_flow = row[4]  	#in vccm
			
			if sample_flow == None: #ignore particles if we can't calculate a volume
				continue
			if (ind_end_time-ind_start_time) > 500:  #ignore particles with a huge sample interval (this arises when the SP2 was set to sample only from 1 of every x minutes)
				continue
			
			
			STP_correction_factor = (self.pressure/101325)*(273.15/self.temperature)
			sample_factor = self._getSampleFactor(ind_end_time)
			particle_sample_vol =  sample_flow*(ind_end_time-ind_start_time)*STP_correction_factor/(60*sample_factor)   #factor of 60 needed because flow is in sccm and time is in seconds
			interval_sampled_volume += particle_sample_vol
			
			VED,rBC_mass,rBC_mass_uncertainty = self._getVEDAndMass(BB_incand_HG,BB_incand_LG,ind_end_time)
			if np.isnan(VED) == False:   #VED and mass are nan if the partcile is outside of the min and max size limits.  In this case we still need to account for the volume of air sampled, but we want to ignore the particle
				interval_mass += rBC_mass
				interval_mass_uncertainty += rBC_mass_uncertainty
				ved_list.append(VED)
			if i % 1000 == 0:
				print i
			i += 1	
		
		interval_data_dict['VED list'] = ved_list
		interval_data_dict['total mass'] = interval_mass
		interval_data_dict['total number'] = len(ved_list)
		interval_data_dict['total mass uncertainty'] = interval_mass_uncertainty
		interval_data_dict['sampled volume'] = interval_sampled_volume
		
		self.assembled_interval_data = interval_data_dict
		

	
	def _binAssembledData(self):
		raw_dia_list = self.assembled_interval_data['VED list']
		total_vol_sccm = self.assembled_interval_data['sampled volume']

		bin_dict = self._makeBinDict()

		for dia in raw_dia_list:
			for point in bin_dict:
				LL_bin = bin_dict[point][0]
				UL_bin = bin_dict[point][1]

				if (LL_bin <= dia < UL_bin):
					mass = ((dia/(10.**7))**3)*(math.pi/6.)*1.8*(10.**15)	
					bin_dict[point][2] += mass
					bin_dict[point][3] += 1
	
		self.binned_data = bin_dict
	
	
	def _writeIntervalDataToDatabase(self):
		
		interval_record ={
		'UNIX_UTC_ts_int_start':self.interval_start,
		'UNIX_UTC_ts_int_end':self.interval_end,
		'instr_ID':self.instr_ID,
		'instr_location_ID':self.instr_location_ID,
		'calibration_ID':self.calibration_ID,
		'total_interval_mass':self.assembled_interval_data['total mass'],
		'total_interval_mass_uncertainty':self.assembled_interval_data['total mass uncertainty'],
		'total_interval_number':self.assembled_interval_data['total number'],
		'total_interval_volume':self.assembled_interval_data['sampled volume'],
		}
		
		add_data = ('''INSERT INTO sp2_interval_data							  
		  (UNIX_UTC_ts_int_start,
		  UNIX_UTC_ts_int_end,
		  instr_ID,
		  instr_location_ID,
		  calibration_ID,
		  total_interval_mass,
		  total_interval_mass_uncertainty,
		  total_interval_number,
		  total_interval_volume)
		  VALUES (
		  %(UNIX_UTC_ts_int_start)s,
		  %(UNIX_UTC_ts_int_end)s,
		  %(instr_ID)s,
		  %(instr_location_ID)s,
		  %(calibration_ID)s,
		  %(total_interval_mass)s,
		  %(total_interval_mass_uncertainty)s,
		  %(total_interval_number)s,
		  %(total_interval_volume)s)''')
		  
	  	#do not duplicate intervals or their associated data!
		
		self.db_cur.execute('''DELETE FROM sp2_interval_data 
			WHERE UNIX_UTC_ts_int_start = %s
			AND UNIX_UTC_ts_int_end = %s
			AND instr_ID = %s
			AND instr_location_ID = %s
			AND calibration_ID = %s
			''',
			(self.interval_start,self.interval_end,self.instr_ID,self.instr_location_ID,self.calibration_ID))

		#add new data
		self.db_cur.execute(add_data, interval_record)
		self.db_connection.commit()
		
		interval_id = self.db_cur.lastrowid
		self.interval_db_id = interval_id
	

	def _writeBinnedDataToDatabase(self):
		
		add_data = ('''INSERT INTO sp2_distribution_data							  
		  (bin_ll,
		  bin_ul,
		  bin_mass,
		  bin_number,
		  total_volume,
		  interval_id)
		  VALUES (
		  %(bin_ll)s,
		  %(bin_ul)s,
		  %(bin_mass)s,
		  %(bin_number)s,
		  %(total_volume)s,
		  %(interval_id)s)''')

		multiple_records = []
		for point in self.binned_data:
				
			single_record ={
			'bin_ll':point[0],
			'bin_ul':point[1],
			'mass':point[2],
			'number':point[3],
			'total_volume':self.assembled_interval_data['sampled volume'],
			'interval_ID':self.interval_db_id
			}
			
			multiple_records.append((single_record))
  
		self.db_cur.executemany(add_data, multiple_records)
		self.db_connection.commit()


	def fit_distr(bin_midpoints,bin_values):
		try:
			popt, pcov = curve_fit(lognorm, np.array(bin_midpoints), np.array(bin_values), p0=(2000,0.6,150))	

		except Exception,e: 
			popt = [np.nan,np.nan,np.nan]
			pcov = [np.nan,np.nan,np.nan]
				
		perr = np.sqrt(np.diag(pcov))
		sigma = math.exp(popt[1])

		return popt,perr,sigma

	
	def _getVEDAndMass(self,BB_incand_HG,BB_incand_LG,ind_end_time):
		
		#get calibration coefficients and detection limits for this calibration (may be detector limits or limits of calibration)
		#HG channel is present in all instrument types	
		HG_lower_dl, HG_upper_dl, HG_calib_0, HG_calib_1, HG_calib_2, HG_calib_0_err, HG_calib_1_err, HG_calib_2_err = self._getCalibrationCoefficients('BBHG_incand',ind_end_time)
		
		#check if this particle has a signal outside the detection limits. If it does, return NaN for the mass and VED
		if self.number_of_channels == 4:
			if (BB_incand_HG <= HG_lower_dl) or (BB_incand_HG >= HG_upper_dl): #check if this particle has a signal outside the detection limits. If it does, return NaN for the mass and VED
				return np.nan, np.nan, np.nan

		#if this is an 8 channel instrument, get the low gain channel detection limits and check if this particle has a signal outside the detection limits. If it does, return NaN for the mass and VED
		if self.number_of_channels == 8:
			LG_lower_dl, LG_upper_dl, LG_calib_0, LG_calib_1, LG_calib_2, LG_calib_0_err, LG_calib_1_err, LG_calib_2_err = self._getCalibrationCoefficients('BBLG_incand',ind_end_time)
			if (BB_incand_HG <= HG_lower_dl) or (BB_incand_LG >= LG_upper_dl):
				return np.nan, np.nan, np.nan

		#if the signal is within the detection limits, continue to calculate mass and VED.  
		#based on the signal and instrument parameters, choose to use low or high gain channel.  The default is to use the high gain channel if possible.
		if BB_incand_HG < HG_upper_dl:
			signal = BB_incand_HG
			rBC_mass = HG_calib_0 + HG_calib_1*signal + HG_calib_2*signal*signal
			rBC_mass_max = (HG_calib_0+HG_calib_0_err)  + (HG_calib_1+HG_calib_1_err)*signal + (HG_calib_2+HG_calib_2_err)*signal*signal
			rBC_mass_uncertainty = (rBC_mass_max - rBC_mass)
		else:
			signal = BB_incand_LG
			rBC_mass = LG_calib_0 + LG_calib_1*signal + LG_calib_2*signal*signal
			rBC_mass_max = (LG_calib_0+LG_calib_0_err)  + (LG_calib_1+LG_calib_1_err)*signal + (LG_calib_2+LG_calib_2_err)*signal*signal
			rBC_mass_uncertainty = (rBC_mass_max - rBC_mass)

		VED = self._calculateVED(rBC_mass)

		return VED,rBC_mass,rBC_mass_uncertainty


	def _getCalibrationCoefficients(self,calibrated_channel,ind_end_time):

		self.db_cur.execute('''
		SELECT 
			0_term,
			1_term,
			2_term,
			0_term_err,
			1_term_err,
			2_term_err,
			lower_dl,
			upper_dl,
			id	
		FROM
			sp2_calibrations
		WHERE
			instr_ID = %s
			AND instr_location_ID = %s
			AND calibrated_channel = %s
			AND calibration_date <= %s
			ORDER BY calibration_date DESC LIMIT 1
			
		''',
		(self.instr_ID,self.instr_location_ID,calibrated_channel,ind_end_time))

		calib_coeffs = self.db_cur.fetchall()

		for row in calib_coeffs:
			calib_0 = row[0]
			calib_1 = row[1]
			calib_2 = row[2]
			calib_0_err = row[3]
			calib_1_err = row[4]
			calib_2_err = row[5]
			lower_dl = row[6]
			upper_dl = row[7]
			calib_ID = row[8]


		self.calibration_ID = calib_ID
		return lower_dl, upper_dl, calib_0, calib_1, calib_2, calib_0_err, calib_1_err, calib_2_err


	def _getSampleFactor(self,ind_end_time):
		
		sample_factor = 1  #default value is 1
		
		#check if a sample factor table exists
		self.db_cur.execute('''SHOW TABLES LIKE 'sp2_sample_factors'
			''')
		table = self.db_cur.fetchall()
		
		#if a sample factor table exists check if there is a sample factor for this interval
		if table != []:
			self.db_cur.execute('''
			SELECT 
				sample_factor
			FROM
				sp2_sample_factors
			WHERE
				UNIX_UTC_int_start <= %s
				AND UNIX_UTC_int_end > %s
				AND instr_locn_ID = %s
				AND instr_ID = %s
			''',
			(ind_end_time,ind_end_time,self.instr_location_ID,self.instr_ID))
			sample_factor_result = self.db_cur.fetchall()
			
			#if there is an entry for this interval, get the sample factor
			if sample_factor_result != []:
				sample_factor = sample_factor_result[0]
			
		return sample_factor
				

	def _makeBinDict(self):
		new_dict = {}
		for bin in range(self.min_rBC_VED,self.max_rBC_VED,self.binning_increment):
			new_dict[bin] = [bin,(bin+self.binning_increment),0,0]
		return new_dict


	def _calculateVED(self,rBC_mass):
		VED = (((rBC_mass/(10**15*self.rBC_density))*6/math.pi)**(1/3.0))*10**7
		return VED
		
	def lognorm(x_vals, A, w, xc):
		return A/(np.sqrt(2*math.pi)*w*x_vals)*np.exp(-(np.log(x_vals/xc))**2/(2*w**2))
		
