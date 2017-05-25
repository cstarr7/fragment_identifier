# -*- coding: utf-8 -*-
# @Author: Charles Starr
# @Date:   2017-05-23 13:33:09
# @Last Modified by:   Charles Starr
# @Last Modified time: 2017-05-25 13:09:00

from spectrum_handler import *
from peptide_fragmenter import *
import bisect
import pandas as pd

class FragmentMatcher(object):
	#this class will take a mass spectrum and theoretical peptide fragments
	#and align them
	def __init__(self, in_file, mass_spectrum, peptide_fragments, accuracy_tol, intensity_tol):

		self.in_file = in_file
		self.mass_spectrum = sorted(mass_spectrum.centroided_peaks, key=lambda tup: tup[0])
		self.masses_only = [peak[0] for peak in self.mass_spectrum]
		self.peptide_fragments = peptide_fragments.fragment_list
		self.accuracy_tol = accuracy_tol
		self.intensity_tol = intensity_tol
		self.potential_proton_matches = self.match_proton_flight()
		self.potential_sodium_matches = self.match_sodium_flight()
		self.potential_potassium_matches = self.match_potassium_flight()

	def match_proton_flight(self):
		#matches fragments with a proton added, the most common flight

		potential_matches = []

		for fragment in self.peptide_fragments:
			adjusted_mass = fragment[1] + 1.007825 #add proton for flight

			frag_matches = self.match_fragments(fragment, adjusted_mass)

			potential_matches += frag_matches

		return potential_matches

	def match_sodium_flight(self):
		#matches fragments with sodium added, another common flight
		potential_matches = []

		for fragment in self.peptide_fragments:
			adjusted_mass = fragment[1] + 22.989770 #add sodium for flight

			frag_matches = self.match_fragments(fragment, adjusted_mass)

			potential_matches += frag_matches

		return potential_matches

	def match_potassium_flight(self):
		#matches fragments with potassium added, uncommon flight
		potential_matches = []

		for fragment in self.peptide_fragments:
			adjusted_mass = fragment[1] + 38.963708 #add potassium for flight

			frag_matches = self.match_fragments(fragment, adjusted_mass)

			potential_matches += frag_matches

		return potential_matches

	def match_fragments(self, fragment, adjusted_mass):

		potential_matches = []
			
		index = bisect.bisect(self.masses_only, adjusted_mass)
		down_index = index - 1
		up_index = index + 1

		if (abs(self.mass_spectrum[index][0] - adjusted_mass) < self.accuracy_tol and 
			self.mass_spectrum[index][1] > self.intensity_tol):
			
			potential_matches.append(
				(self.mass_spectrum[index][0],
					self.mass_spectrum[index][1], fragment[0], fragment[1]))

		while abs(self.mass_spectrum[down_index][0] - adjusted_mass) < self.accuracy_tol:
			if self.mass_spectrum[down_index][1] > self.intensity_tol:
				potential_matches.append(
				(self.mass_spectrum[down_index][0],
					self.mass_spectrum[down_index][1], fragment[0], fragment[1]))

			down_index -= 1

		while abs(self.mass_spectrum[up_index][0] - adjusted_mass) < self.accuracy_tol:
			if self.mass_spectrum[up_index][1] > self.intensity_tol:
				potential_matches.append(
				(self.mass_spectrum[up_index][0], 
					self.mass_spectrum[up_index][1], fragment[0], fragment[1]))

			up_index += 1

		return potential_matches

	def export_results(self):

		proton_frame = pd.DataFrame(self.potential_proton_matches, columns=['Observed Mass', 'Intensity', 'Sequence', 'Exact Mass'])
		sodium_frame = pd.DataFrame(self.potential_sodium_matches, columns=['Observed Mass', 'Intensity', 'Sequence', 'Exact Mass'])
		potassium_frame = pd.DataFrame(self.potential_potassium_matches, columns=['Observed Mass', 'Intensity', 'Sequence', 'Exact Mass'])
		
		outfile = pd.ExcelWriter(self.in_file[:-5] + '_matchestest.xlsx')
		proton_frame.to_excel(outfile, sheet_name = 'Proton')
		sodium_frame.to_excel(outfile, sheet_name = 'Sodium')
		potassium_frame.to_excel(outfile, sheet_name = 'Potassium')
		outfile.save()


'''
		test_frame = pd.DataFrame(self.masses_only)
		test_outfile = pd.ExcelWriter('test_out.xlsx')
		test_frame.to_excel(test_outfile)
		test_outfile.save()

		test_frametwo = pd.DataFrame(self.mass_spectrum)
		test_outfiletwo = pd.ExcelWriter('test_outtwo.xlsx')
		test_frametwo.to_excel(test_outfiletwo)
		test_outfile.save()
'''



