# -*- coding: utf-8 -*-
# @Author: Charles Starr
# @Date:   2017-05-23 11:47:10
# @Last Modified by:   Charles Starr
# @Last Modified time: 2017-05-23 14:58:38
import pymzml

class MassSpectrum(object):
	#this class will hold the data for a mass spec experiment

	def __init__(self, ms_datafile):

		self.ms_datafile = ms_datafile
		self.ms_spectrum = self.parse_mzml()
		self.centroided_peaks = self.ms_spectrum.next().centroidedPeaks

	def parse_mzml(self):
		# Instructions for parsing an mzml file.

		msrun = pymzml.run.Reader(self.ms_datafile)

		return msrun



