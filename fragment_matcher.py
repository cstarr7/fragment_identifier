# -*- coding: utf-8 -*-
# @Author: Charles Starr
# @Date:   2017-05-23 13:33:09
# @Last Modified by:   Charles Starr
# @Last Modified time: 2017-05-23 14:21:03

from spectrum_handler import *
from peptide_fragmenter import *
import bisect
import pandas as pd

class FragmentMatcher(object):
	#this class will take a mass spectrum and theoretical peptide fragments
	#and align them
	def __init__(self, in_file, mass_spectrum, peptide_fragments, accuracy_tol, intensity_tol):

		self.in_file = in_file
		self.mass_spectrum = mass_spectrum
		self.masses_only = [peak[0] for peak in mass_spectrum]
		self.peptide_fragments = peptide fragments
		self.accuracy_tol = accuracy_tol
		self.intensity_tol = intensity_tol
		self.potential_matches = self.match_fragments()

	def match_fragments(self):

		potential_matches = []

		for fragment in self.peptide_fragments:
			fragment_mass = fragment[1] - 1.007825 #subtract proton from flight
			
			index = bisect.bisect(self.masses_only, fragment_mass)
			down_index = index - 1
			up_index = index + 1

			if abs(self.mass_spectrum[index][0] - fragment_mass) < self.accuracy_tol and 
				self.mass_spectrum[index][1] > self.intensity_tol:
				
				potential_matches.append(
					(self.mass_spectrum[index][0],
						self.mass_spectrum[index][1], fragment[0], fragment[1]))

			while abs(self.mass_spectrum[down_index][0] - fragment_mass) < self.accuracy_tol:
				if self.mass_spectrum[down_index][1] > self.intensity_tol:
					potential_matches.append(
					(self.mass_spectrum[down_index][0],
						self.mass_spectrum[down_index][1], fragment[0], fragment[1]))

				down_index -= 1

			while abs(self.mass_spectrum[up_index][0] - fragment_mass) < self.accuracy_tol:
				if self.mass_spectrum[up_index][1] > self.intensity_tol:
					potential_matches.append(
					(self.mass_spectrum[up_index][0], 
						self.mass_spectrum[up_index][1], fragment[0], fragment[1]))

				up_index -= 1

		return potential_matches

	def export_results(self):

		frame = pd.DataFrame(self.potential_matches, columns=['Observed Mass', 'Intensity', 'Sequence', 'Exact Mass'])
		outfile = pd.ExcelFile(self.in_file + '_matches')
		frame.to_excel(outfile)
		outfile.save()







