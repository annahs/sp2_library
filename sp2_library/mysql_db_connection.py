#!/usr/bin/env python
# -*- coding: utf-8 -*-
import mysql.connector

class dbConnection():

	"""

	Creates and manages a mysql database connection and cursor

	"""

	def __init__(self, schema):
	    self.db_connection = mysql.connector.connect(user='root', password='', host='localhost', database=schema) 
	    self.db_cur = self.db_connection.cursor()

	def __del__(self):
	    self.db_connection.close()


