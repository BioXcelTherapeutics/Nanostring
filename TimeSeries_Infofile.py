#!/usr/bin/env python

"""

	Modules to use NIH/NLM E-Direct to mine PubMed

"""

__author__ =        "Nirmal Keshava"
__copyright__ =     "Copyright 2018, BioXcel Therapeutics"
__credits__ =       []
__license__ =       "GPL"
__version__ =       "1.0"
__maintainer__ =    "Nirmal Keshava"
__email__ =         "nkeshava@bioxceltherapeutics.com"
__status__ =        "R & D"

###########################################################################









###########################################################################
#
# Experiment 1
#
# #########################################################################

inputNormalizedCsvFile = '/data/sandbox/Nanostring/18_11_09_PkPd_TimeCourse/2018 10 10 Bioxcel Normalized.csv'
inputRawCsvFile = '/data/sandbox/Nanostring/18_11_09_PkPd_TimeCourse/2018 10 10 Bioxcel Raw.csv'

outputFolder = '/home/nirmal/Output/Exp_1/'

# keywords = [ 'adiporon' ] #, 'florida', 'EGFR', 'lung', 'smoking', 'metastasis' ]
keywords = 'diabetes[title]'

sortOrder = 'First Author'
dateType = 'PDAT'
minDate = '0'
maxDate = '2018'
retMax = '10000'

drugRepoFile = '/data/sandbox/BroadRepurposingHub/repurposing_drugs_20180907.xlsx'

