# -*- coding: utf-8 -*-
# @Author: Charles Starr
# @Date:   2017-05-23 14:22:52
# @Last Modified by:   Charles Starr
# @Last Modified time: 2017-05-23 14:28:51

from spectrum_handler import *
from peptide_fragmenter import *
from fragment_matcher import *


def main(mass_data_file, peptide_sequence, accuracy, intensity):

	spectrum = MassSpectrum(mass_data_file)
	fragments = PeptideFragmenter(peptide_sequence)
	matcher = FragmentMatcher(mass_data_file, spectrum, fragments, accuracy, intensity)
	matcher.export_results()

main('ARVA_7_CYTOSOL%2033%25.mzML', 'RRGWALRLVLAY', 0.5, 100.0)