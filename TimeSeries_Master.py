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

from TimeSeries_Infofile import *
# import biopython

def ReadXlsFile(filename, masterSheetNum = 2, headerRowNum = 1, topRowBlank = True, firstColBlank = True, numLines = 0):
	"""

		Method to read in a .xls file

	:param filename:
	:param opFilename:
	:return:
	"""

	import xlrd

	# Open up xls file
	# book = xlrd.open_workbook(filename, encoding_override="utf_8")
	book = xlrd.open_workbook(filename)
	master_sheet = book.sheet_by_index(masterSheetNum)

	# Get number of rows and columns
	Nrows = master_sheet.nrows
	Ncols = master_sheet.ncols

	# Handle top row
	sheet = list()
	if topRowBlank == True:
		startRow = 1
	else:
		startRow = 0

	if numLines != 0:
		endRow = startRow + numLines
	else:
		endRow = Nrows

	# Handle first column
	if firstColBlank == True:
		startCol = 1
	else:
		startCol = 0

	# Iterate over master_sheet
	sheet = list()
	for ii in range(startRow, endRow):
		if ii == headerRowNum:
			# Grab header
			header = list()
			for jj in range(startCol, Ncols):
				curVal = master_sheet.cell_value(ii,jj)
				header.append(curVal)
				# header.append(curVal.encode('utf8'))
		else:
			sheet.append([])
			for jj in range(startCol, Ncols):
				# print ii,jj, master_sheet.cell_value(ii,jj)
				curVal = master_sheet.cell_value(ii,jj)
				# print ii, jj, curVal
				# sheet[-1].append(curVal)
				if isinstance(curVal, int):
					# Looking for floating points
					sheet[ii].append(int(curVal))
				elif isinstance(curVal, unicode):
					# Looking for unicode strings
					sheet[-1].append(curVal.encode('utf8'))
				# elif isinstance(curVal, float) and curVal > 40000:
				#     # Looking for dates
				#     sheet[-1].append(xlrd.xldate.xldate_as_datetime(curVal,book.datemode))
				elif isinstance(curVal, float):
					# Looking for floating points
					sheet[-1].append(curVal)
				# else:
				#     # There is nothing in that column
				#     sheet[-1].append('')

	return header, sheet

def ReadCsvFile(csvFilename):
	"""

		Method to read .CSV file in the format in which the Nanostring time-series
		gene expression data is delivered.

		Note, row0 is the header, row 1 and 2 are not useful, and col0 contains
		the gene names

	:param csvFilename:     <string> containig csv file name
	:return: header:        <list> of strings, for each column heading
			 csvData:       <list> of lists containing .csv data for each gene
	"""

	import csv

	# Create lists for header and data
	header = list()
	csvData = list()

	# Open file
	with open(csvFilename, 'r') as csvfile:
		reader = csv.reader(csvfile, delimiter=',')
		rowCtr = 0
		for row in reader:
			if rowCtr == 0:
				# Digest header row
				for kk in range(0,len(row)):
					header.append(row[kk])
				rowCtr += 1
			elif rowCtr == 1 or rowCtr == 2:
				rowCtr += 1
			else:
				# Capture rows as list
				curRow = list()
				for kk in range(0, len(row)):
					curRow.append(row[kk])

				# Current row, reformatted
				curRowR = list()
				for kk in range(len(curRow)):
					if kk <= 6:
						curRowR.append(curRow[kk])
					else:
						curRowR.append(float(curRow[kk]))

				# Append reformatted row
				csvData.append(curRowR)

	return header, csvData[1:]

def GetMouseTimeCourseData(csvFilename):
	"""

		Method to format CSV data into format for digesting time-course data, including
		mouse replicates

	:param ccsvFilename:            <string> containig csv file name
	:return: geneList:              <list> of strings containing gene names
			 geneTimeCourseData     <list> of sub-lists where each list entry is for a
			                        different gene, and each sub-list entry contains the
			                        gene expression values for each time-point
	"""

	# Read in CSV file
	csvHeader, csvData = ReadCsvFile(csvFilename)

	# Extract list of gene names
	geneList = [ xx[0] for ii,xx in enumerate(csvData) ]
	Ngenes = len(geneList)

	# Extract time-course gene expression values for each gene
	geneTimeCourseData = list()
	for ii in range(Ngenes):
		curGeneData = csvData[ii]
		curGeneTimeCourseData = list()
		for jj in range(7):
			curGeneVals = [ curGeneData[7+jj*3], curGeneData[7+jj*3 + 1], curGeneData[7+jj*3 + 2] ]
			curGeneTimeCourseData.append(curGeneVals)
		geneTimeCourseData.append(curGeneTimeCourseData)

	return geneList, geneTimeCourseData

def PlotGeneTimeCourseProgressions(inputCsvFile, figsFolder):
	"""

		Method to plot time-course data for each gene

	:param inputNormalizedCsvFile:  <string> containig csv file name
	:param figsFolder:              <string> to folder where figures are to be written
	:return:                        NULL
	"""

	from matplotlib import pyplot as plt
	from numpy import log

	geneList, geneTimeCourseData = GetMouseTimeCourseData(inputCsvFile)

	Ngenes = len(geneList)

	for ii in range(Ngenes):
		fig = plt.figure()

		ax1 = fig.add_subplot(2, 1, 1)
		plt.boxplot(geneTimeCourseData[ii])
		ax1.set_title(geneList[ii])
		ax1.set_ylabel('Normalized Counts')
		ax1.set_xlabel('Hours Since Dosing')
		ax1.set_xticks([ 1, 2, 3, 4, 5, 6, 7])
		ax1.set_xticklabels(['0', '1', '2', '4', '8', '16', '24'])
		plt.grid()

		logTimeCourseData = list()
		for jj in range(len(geneTimeCourseData[ii])):
			logTimeCourseData.append([])
			for kk in range(len(geneTimeCourseData[ii][jj])):
				logTimeCourseData[jj].append(log(geneTimeCourseData[ii][jj][kk]))

		ax2 = fig.add_subplot(2, 1, 2)
		plt.boxplot(logTimeCourseData)
		ax2.set_title(geneList[ii])
		ax2.set_ylabel('log(Normalized Counts)')
		ax2.set_xlabel('Hours Since Dosing')
		ax2.set_xticks([ 1, 2, 3, 4, 5, 6, 7])
		ax2.set_xticklabels(['0', '1', '2', '4', '8', '16', '24'])

		fig.tight_layout()
		plt.grid()

		plt.savefig(figsFolder + 'TimeCourse_' + geneList[ii] + '.png')
		plt.close()
		# plt.show()

	return

def main():

	# Read in PkPd file
	# header, csvData = ReadCsvFile(inputNormalizedCsvFile)

	PlotGeneTimeCourseProgressions(inputNormalizedCsvFile, outputFolder + 'Figs/')

	# geneList, geneTimeCourseData = GetMouseTimeCourseData(inputNormalizedCsvFile)


	return

main()
