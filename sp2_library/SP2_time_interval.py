# -*- coding: UTF-8 -*-
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
import SP2_utilities
from mysql_db_connection import dbConnection

class TimeInterval(object,dbConnection):

	"""

	This class represents a time interval in which data was collected by an SP2 (Droplet Measurement Technolgies Inc).
	It has methods for importing, filtering, and analyzing data collected within the interval.

	"""
	
	def __init__(self,database_name,instr_location_ID, instr_ID, interval_start, interval_end):
		dbConnection.__init__(self, database_name)
		

		self.instr_location_ID 			= instr_location_ID 
		self.instr_ID     				= instr_ID
		self.interval_start 			= interval_start
		self.interval_end   			= interval_end
		self.temperature    			= 273.15  	#default is STP
		self.pressure       			= 101325  	#default is STP
		self.altitude					= 1 		#default is 1m above amsl
		self.rBC_density 				= 1.8 		#g/mol - Bond and Bergstrom 2006
		self.interval_max				= 500.   	#maximum time between particles (in sec) 
		self.extrapolate_calibration 	= False
		
		self.calibration_ID = None
		self.assembled_interval_data = None
		self.binned_data = None
		self.interval_db_id = None

		self.retrieveInstrInfo()				#get basic instrument info
		self.retrieveSampleFactors()			#get all sample factors for this interval
		self.retrieveCalibrationData()			#get HG and LG calibration info for this interval
		self.retrieveHousekeepingLimits()		#get values for QC based on housekeeping parameters



	def retrieveHousekeepingLimits(self):
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
			AND instr_location_ID = %s
			AND UNIX_UTC_ts_int_start <= %s
			AND UNIX_UTC_ts_int_end > %s
		LIMIT 1
		''',
		(self.instr_ID, self.instr_location_ID,self.interval_start,self.interval_end))
		hk_limits = self.db_cur.fetchall()

		self.yag_min 		 = hk_limits[0][0]
		self.yag_max 		 = hk_limits[0][1]
		self.sample_flow_min = hk_limits[0][2]
		self.sample_flow_max = hk_limits[0][3]


	def retrieveInstrInfo(self):
		self.db_cur.execute('''
		SELECT 
			number_of_channels,
			min_detectable_signal,
			saturation_limit
		FROM
			sp2_instrument_info 
		WHERE
			id = %s
			AND id > %s
		LIMIT 1
		''',
		(self.instr_ID,0))
		
		instr_info = self.db_cur.fetchall()

		self.number_of_channels 	= instr_info[0][0]
		self.min_detectable_signal 	= instr_info[0][1]
		self.saturation_limit 		= instr_info[0][2]


	def retrieveSampleFactors(self):
		
		#set default
		sample_factors = [(self.interval_start,self.interval_end,1)]

		self.db_cur.execute('''
		SELECT 
			UNIX_UTC_ts_int_start,
			UNIX_UTC_ts_int_end,
			1_in_x_particles
		FROM
			sp2_config_parameters
		WHERE
			UNIX_UTC_ts_int_start <= %s
			AND UNIX_UTC_ts_int_end > %s
			AND instr_location_ID = %s
			AND instr_ID = %s
		''',
		(self.interval_start,self.interval_end,self.instr_location_ID,self.instr_ID))
		sample_factor_results = self.db_cur.fetchall()
		
		#if no results, then use the default value
		if sample_factor_results != []:
			sample_factors = sample_factor_results

		self.sample_factors = sample_factors

	
	def retrieveCalibrationData(self):
		
		calibration_data = {}
		for channel in ['BBHG_incand','BBLG_incand']:
			self.db_cur.execute('''
			SELECT 
				0_term,
				1_term,
				2_term,
				0_term_err,
				1_term_err,
				2_term_err,
				calibration_material,
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
			(self.instr_ID,self.instr_location_ID,channel,self.interval_start))

			calib_coeffs = self.db_cur.fetchall()
			if calib_coeffs == []:
				calib_coeffs_np = [[np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,'nan',np.nan]]
			else:
				calib_coeffs_np = np.array(calib_coeffs, dtype=[('term0', 'f4'),('term1', 'f4'),('term2', 'f4'),('term0err', 'f4'),('term1err', 'f4'),('term2err', 'f4'),('mat', 'S7'),('ID', 'f4'),])  #converts Nones to nans for calculations

			#Aqudag correction
			for row in calib_coeffs_np:
				calib_material 	= row[6]
				calib_ID 		= row[7]
				calib_0 		= row[0]
				calib_0_err		= row[3]
				if calib_material == 'Aquadag':
					calib_1 	= row[1]/0.7
					calib_1_err = row[4]/0.7
					calib_2 	= row[2]/0.7
					calib_2_err = row[5]/0.7
					
			#set calibration ids		
			if channel == 'BBHG_incand':
				self.HG_calibration_ID = float(calib_ID)
			if channel == 'BBLG_incand':
				self.LG_calibration_ID = float(calib_ID)

			#get the signal limits for calculating mass
			if self.extrapolate_calibration == False:
				pkht_ll, pkht_ul = self._retrieveCalibrationLimits(calib_ID)
			else:
				pkht_ll = self.min_detectable_signal
				pkht_ul = self.saturation_limit

			calibration_data[channel] = [pkht_ll, pkht_ul, calib_0, calib_1, calib_2, calib_0_err, calib_1_err, calib_2_err]

		self.calibration_info = calibration_data


	def _retrieveCalibrationLimits(self,calib_ID):
		
		if np.isnan(calib_ID):
			return np.nan, np.nan

		self.db_cur.execute('''
		SELECT 
			min(incand_pk_ht),
			max(incand_pk_ht)
		FROM
			sp2_calibration_points
		WHERE
			calibration_ID = %s
			and id >= %s
		LIMIT 1	
		''',
		(float(calib_ID),0))
		calib_points = self.db_cur.fetchall()

		pkht_ll = calib_points[0][0]
		pkht_ul = calib_points[0][1]
		
		return pkht_ll, pkht_ul		

	
	def retrieveSingleParticleData(self):
		"""
		Get the single particle data for this interval. Use housekeeping information to exclude periods of poor instrument performance .
		"""
		self.db_cur.execute('''
		SELECT 
			sp.UNIX_UTC_ts_int_start,
			sp.UNIX_UTC_ts_int_end,
			sp.BB_incand_HG_pkht,
			sp.BB_incand_LG_pkht,
			hk.sample_flow,
			sp.NB_incand_HG_pkht,
			hk.chamber_temp,
			hk.chamber_pressure
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
		
	
	def getParticleSampleFactor(self,ind_end_time):
		"""
		Get the sample factor for this particle (i.e the fraction of particles recorded.  To save space the SP2 can be configured to record only 1/X particles detected)
		
		Parameters
		----------
		ind_end_time : float 
			Time at which the particle record was written to file
		"""
		if len(self.sample_factors) == 1:
			sample_factor = self.sample_factors[0][2]
		else:
			for row in self.sample_factors:
				sf_start = row[0]
				sf_end   = row[1]
				if sf_start <= ind_end_time < sf_end:
					sample_factor = row[2]
					break

		return sample_factor



	#Interval methods	
	def assembleIntervalData(self):
		"""
		Assemble the interval data.  
		This is stored as a dictionary with: the total interval rBC mass, the uncertainty in total rBC mass, rBC particle number, total sampled volume, and a list of the diameters of detected particles.
		"""
		interval_data_dict = {}
		interval_sampled_volume = 0
		interval_mass = 0
		interval_mass_uncertainty = 0
		ved_list = []

		for row in self.single_particle_data:
			ind_start_time 	= row[0] 	#UNIX UTC timestamp
			ind_end_time 	= row[1]	#UNIX UTC timestamp
			BB_incand_HG 	= row[2]  	#in arbitrary units
			BB_incand_LG 	= row[3]  	#in arbitrary units
			sample_flow 	= row[4]  	#in vccm
			chamber_temp 	= row[6]+273.15 	#in deg C -> K
			chamber_pressure= row[7]  			#in Pa
			
			if sample_flow == None: #ignore particles if we can't calculate a volume
				continue
			if (ind_end_time-ind_start_time) > self.interval_max  or (ind_end_time-ind_start_time) < 0:  #ignore particles with a huge sample interval (this arises when the SP2 was set to sample only from 1 of every x minutes)
				continue
			
			#get appropriate sample factor
			sample_factor = self.getParticleSampleFactor(ind_end_time)
			STP_correction_factor = (chamber_pressure/101325)*(273.15/chamber_temp) 

			particle_sample_vol =  sample_flow*(ind_end_time-ind_start_time)*STP_correction_factor/(60*sample_factor)   #factor of 60 needed because flow is in sccm and time is in seconds
			interval_sampled_volume += particle_sample_vol

			rBC_mass,rBC_mass_uncertainty = self.calculateMass(BB_incand_HG,BB_incand_LG,ind_end_time)
			VED = SP2_utilities.calculateVED(self.rBC_density,rBC_mass)

			if self.min_VED <= VED <= self.max_VED:  #we limit mass and number concentrations to within the set size limits
				interval_mass += rBC_mass
				interval_mass_uncertainty += rBC_mass_uncertainty
				ved_list.append(VED)
				
		interval_data_dict['VED list'] = ved_list
		interval_data_dict['total mass'] = interval_mass
		interval_data_dict['total number'] = len(ved_list)
		interval_data_dict['total mass uncertainty'] = interval_mass_uncertainty
		interval_data_dict['sampled volume'] = interval_sampled_volume

		self.assembled_interval_data = interval_data_dict


	#Binned data methods
	def binAssembledData(self,binning_increment):
		"""
		Bin the assembled interval data.
		
		Parameters
		----------
		binning_increment : int 
			bin width
		"""
		raw_dia_list = self.assembled_interval_data['VED list']
		total_vol_sccm = self.assembled_interval_data['sampled volume']
		self.binning_increment = binning_increment
		
		bin_dict = self.makeBinDict()

		for dia in raw_dia_list:
			for point in bin_dict:
				LL_bin = bin_dict[point][0]
				UL_bin = bin_dict[point][1]

				if (LL_bin <= dia < UL_bin):
					mass = SP2_utilities.calculateMass(self.rBC_density,dia)	
					bin_dict[point][2] += mass
					bin_dict[point][3] += 1
	
		self.binned_data = bin_dict
	

	def lognormFit(self,bin_midpoints,bin_values):
		"""
		Fit a single lognormal function to the binned data.
		
		Parameters
		----------
		bin_midpoints : list of floats 
			bin midpoints
		bin_values : list of floats
			bin values
		"""
		try:
			popt, pcov = curve_fit(SP2_utilities.lognorm, np.array(bin_midpoints), np.array(bin_values), p0=(2000,0.6,150))	

		except Exception,e: 
			popt = [np.nan,np.nan,np.nan]
			pcov = [np.nan,np.nan,np.nan]
				
		perr = np.sqrt(np.diag(pcov))
		sigma = math.exp(popt[1])

		return popt,perr,sigma
	

	def setBinningLimits(self):
		"""
		Set the binning limits based on the calibration range
		"""
		
		HG_pkht_ll, HG_pkht_ul, HG_calib_0, HG_calib_1, HG_calib_2, HG_calib_0_err, HG_calib_1_err, HG_calib_2_err = self.calibration_info['BBHG_incand']
		min_detectible_mass = np.nansum([HG_calib_0,(HG_calib_1*HG_pkht_ll),(HG_calib_2*HG_pkht_ll*HG_pkht_ll)])
		max_detectible_mass = np.nansum([HG_calib_0,(HG_calib_1*HG_pkht_ul),(HG_calib_2*HG_pkht_ul*HG_pkht_ul)])

		if self.number_of_channels == 8:
			LG_pkht_ll, LG_pkht_ul, LG_calib_0, LG_calib_1, LG_calib_2, LG_calib_0_err, LG_calib_1_err, LG_calib_2_err = self.calibration_info['BBLG_incand']
			max_detectible_mass = np.nansum([LG_calib_0,(LG_calib_1*LG_pkht_ul),(LG_calib_2*LG_pkht_ul*LG_pkht_ul)])

		#min and max VEDs for binning
		self.min_VED = int(math.floor(SP2_utilities.calculateVED(self.rBC_density,min_detectible_mass))) #floor - gives nearest integer below value, binning range is slightly wider than data
		self.max_VED = int(math.ceil(SP2_utilities.calculateVED(self.rBC_density,max_detectible_mass))) #ceil - gives nearest integer above value, as above



	def makeBinDict(self):
		"""
		Make the binning dictionary
		"""
		new_dict = {}
		for bin in range(self.min_VED,(self.max_VED+self.binning_increment),self.binning_increment):
			new_dict[bin] = [bin,(bin+self.binning_increment),0,0]
		return new_dict


	#auxiliary methods
	def calculateMass(self,BB_incand_HG,BB_incand_LG,ind_end_time):
		"""
		Calculate the mass and uncertainty for rBC in a single particle
		
		Parameters
		----------
		BB_incand_HG : float
			Broadband high-gain incandescence channel signal height
		BB_incand_LG : float
			Broadband low-gain incandescence channel signal height
		ind_end_time : float 
			Time at which the particle record was written to file
		"""
		
		#get calibration coefficients and detection limits for this calibration (may be detector limits or limits of calibration)
		#HG channel is present in all instrument types	
		HG_pkht_ll, HG_pkht_ul, HG_calib_0, HG_calib_1, HG_calib_2, HG_calib_0_err, HG_calib_1_err, HG_calib_2_err = self.calibration_info['BBHG_incand']
		#check if this particle has a signal outside the detection limits. If it does, return NaN for the mass and VED
		if self.number_of_channels == 4:			
			if (BB_incand_HG <= HG_pkht_ll) or (BB_incand_HG >= HG_pkht_ul): #check if this particle has a signal outside the calibration limits. If it does, return NaN for the mass and VED
				return np.nan, np.nan

		#if this is an 8 channel instrument, get the low gain channel detection limits and check if this particle has a signal outside the detection limits. If it does, return NaN for the mass and VED
		if self.number_of_channels == 8:
			LG_pkht_ll, LG_pkht_ul, LG_calib_0, LG_calib_1, LG_calib_2, LG_calib_0_err, LG_calib_1_err, LG_calib_2_err = self.calibration_info['BBLG_incand']
			if (BB_incand_HG <= HG_pkht_ll) or (BB_incand_LG >= LG_pkht_ul):
				return np.nan, np.nan

		#if the signal is within the detection limits, continue to calculate mass 
		#based on the signal and instrument parameters, choose to use low or high gain channel.  The default is to use the high gain channel if possible.
		if BB_incand_HG < HG_pkht_ul:
			signal = BB_incand_HG
			rBC_mass = np.nansum([HG_calib_0,(HG_calib_1*signal),(HG_calib_2*signal*signal)])
			rBC_mass_max = np.nansum([(HG_calib_0+HG_calib_0_err),(HG_calib_1+HG_calib_1_err)*signal,(HG_calib_2+HG_calib_2_err)*signal*signal])
			rBC_mass_uncertainty = (rBC_mass_max - rBC_mass)*1.

		else:
			signal = BB_incand_LG
			rBC_mass = np.nansum([LG_calib_0,(LG_calib_1*signal),(LG_calib_2*signal*signal)])
			rBC_mass_max = np.nansum([(LG_calib_0+LG_calib_0_err),(LG_calib_1+LG_calib_1_err)*signal,(LG_calib_2+LG_calib_2_err)*signal*signal])
			rBC_mass_uncertainty = (rBC_mass_max - rBC_mass)*1.

		return rBC_mass,rBC_mass_uncertainty




