# -*- coding: utf-8 -*-
# @Author: Charles Starr
# @Date:   2017-05-23 14:22:52
# @Last Modified by:   Charles Starr
# @Last Modified time: 2017-05-24 15:13:19

from spectrum_handler import *
from peptide_fragmenter import *
from fragment_matcher import *


def main(mass_data_file, peptide_sequence, accuracy, intensity):

	spectrum = MassSpectrum(mass_data_file)
	fragments = PeptideFragmenter(peptide_sequence)
	matcher = FragmentMatcher(mass_data_file, spectrum, fragments, accuracy, intensity)
	matcher.export_results()

main('indol_8.mzML', 'ILPWKWPWWPWRR', 0.5, 100.0)