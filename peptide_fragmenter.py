# -*- coding: utf-8 -*-
# @Author: Charles Starr
# @Date:   2017-05-23 12:52:51
# @Last Modified by:   Charles Starr
# @Last Modified time: 2017-05-23 23:10:18

from atomic_mass_data import *
import pandas as pd


class PeptideFragmenter(object):
	#this class takes a peptide and breaks it into all possible pieces
	def __init__(self, sequence):

		self.peptide_sequence = sequence
		self.fragment_list = self.create_fragments()
		
		
		writer = pd.ExcelWriter('testthree.xlsx')
		frame = pd.DataFrame(self.fragment_list)
		frame.to_excel(writer)
		writer.save()
		

	def create_fragments(self):
	#do the peptide fragmenting sequentially
		frags = []
		for i in range(len(self.peptide_sequence)):
			#this higher order loop chops the n-terminus and calculates mass for that peptide
			new_frag = self.peptide_sequence[i:]
			mass = self.exact_mass_calculator(new_frag)
			mass += termini_map['c_amide'] # add amidated C-terminus
			frags.append((new_frag, mass))

			for j in range(1, len(self.peptide_sequence[i:])):
				#this lower order loop chops sequentially from the c-terminus
				new_frag = self.peptide_sequence[i:-j]
				mass = self.exact_mass_calculator(new_frag)
				mass += termini_map['c_carboxyl'] # add normal C-terminus
				frags.append((new_frag, mass))
		return frags

	def residue_mass(self, sequence):
        # Handles consistent masses for each sequence. First calculates
        # backbone mass based on length of sequence, then appends each
        # residue mass.
		
		exact_mass = (len(sequence) - 1) * termini_map['bond']
        
		for amino_acid in sequence:
		    exact_mass += amino_acid_map[amino_acid]

		return exact_mass

	def exact_mass_calculator(self, sequence):
        # Calculates the mw weight for a peptide represented as a 
        # string of amino acids
		
		exact_mass = self.residue_mass(sequence)
		exact_mass += termini_map['normal_n'] # add N-terminus mass

		return exact_mass




PeptideFragmenter('RRGWALRLVLAY')
