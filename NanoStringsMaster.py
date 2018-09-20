
#!/usr/bin/env python

# This script was written to post-process measurements performed by a Nanostring system which
# counts the number of matching RNA base-pair patterns (ie, fingerprints) in a sample
#
# This script processes the results for several distinct drug pair comparisons and for each gene
# performs a t-test between the assay msmts for each drug pair.  If the p-value < 0.05, the gene
# is flagged as being differentially expressed for that drug pair

import xlrd
import numpy as np
from scipy.stats import ttest_ind
import csv


def main():

    # Input folder and file
    dataFolder = 'C:/Users/Nirmal Keshava/OneDrive - BioXcel Therapeutics/Data/Immuno-Oncology/'
    dataFile = 'Canopy normalized data.xlsx'
    dataPath = dataFolder + dataFile

    # Output file
    opFileGenes = 'C:/Users/Nirmal Keshava/OneDrive - BioXcel Therapeutics/Output/701/Nanostring/Nanostring_Differential_Genes.csv'
    opFilePvals = 'C:/Users/Nirmal Keshava/OneDrive - BioXcel Therapeutics/Output/701/Nanostring/Nanostring_Differential_Pvals.csv'

    # Open .xls workbook
    g = xlrd.open_workbook(dataFolder + dataFile)
    sheet1 = g.sheet_by_index(0)

    # Iterate thru rows in sheet1 and extract data, row by row
    data = list()
    for ii in range(sheet1.nrows):

        if ii == 0:
            header = sheet1.row_values(ii)
        else:
            curData = sheet1.row_values(ii)
            curDataEnc = list()
            for jj in range(len(curData)):
                if type(curData[jj]) ==unicode:
                    curDataEnc.append(curData[jj].encode())
                else:
                    curDataEnc.append(curData[jj])

            data.append(curDataEnc)

    #
    # Plot mean values for each gene, for vehicle, 701, anti-pd-1, and combo
    #

    Ngenes = len(data)
    Ndrugs = 4
    NdrugComps = 6
    geneList = list()
    avgVals = np.zeros((Ngenes,4))
    extremeValInds = np.zeros((Ngenes,2))
    avgVehicleInd = 30
    avg701Ind = 31
    avgAntiPd1Ind = 32
    avgComboInd = 33

    # Iterate thru each gene and get names of genes, average value of measurements and extreme values
    for ii in range(Ngenes):
        # Get gene name
        geneList.append(data[ii][0])
        # Get average values
        avgVals[ii,0] = data[ii][avgVehicleInd]
        avgVals[ii,1] = data[ii][avg701Ind]
        avgVals[ii,2] = data[ii][avgAntiPd1Ind]
        avgVals[ii,3] = data[ii][avgComboInd]
        # Tally which drug had highest and lowest value
        extremeValInds[ii,0] = np.argmin(avgVals[ii,:])
        extremeValInds[ii,1] = np.argmax(avgVals[ii,:])

    # Assess genes in which there are significant deviations in expression levels between the drugs in each drug pair
    alpha = 0.05
    drugInds = [ [17,18,19], [20,21,22], [23,24,25], [26,27,28] ]
    geneComps = np.zeros((Ngenes,NdrugComps))
    geneCompsRawPval = np.zeros((Ngenes,NdrugComps))
    geneCompsRawTstat = np.zeros((Ngenes,NdrugComps))
    geneCompStrs = [ 'Vehicle vs. 701', 'Vehicle vs. AntiPD1', 'Vehicle vs. Combo', '701 vs. AntiPD1', '701 vs. Combo', 'AntiPD1 vs. Combo' ]

    # Iterate thru each gene
    for ii in range(Ngenes):

        # Get msmts for current gene
        curGeneData = data[ii]

        # For each pair, evaluate significance test to determine which genes are expressed differently
        drugCompCtr = 0
        for pairInd0 in range(Ndrugs-1):
            curPairInds0 = drugInds[pairInd0]
            # Extract msmts for current gene and drug pair member 0
            curPairData0 = [ xx for jj,xx in enumerate(curGeneData) if jj in curPairInds0 ]

            for pairInd1 in range(pairInd0+1,Ndrugs):
                curPairInds1 = drugInds[pairInd1]
                # Extract msmts for current gene and drug pair member 1
                curPairData1 = [ xx for jj,xx in enumerate(curGeneData) if jj in curPairInds1 ]

                # Perform significance test to assess whether curPairData0 and curPairData1 are different
                tstat, pval = ttest_ind(curPairData1, curPairData0)

                # Record Pval
                geneCompsRawPval[ii,drugCompCtr] = pval
                geneCompsRawTstat[ii,drugCompCtr] = tstat

                # Record whether there was a significant difference, and if so, the direction of the
                # difference
                if pval < alpha:
                    if tstat < 0:
                        # First member of drug pair has larger value than second:  second drug is *underexpressed*
                        geneComps[ii,drugCompCtr] = 1
                    else:
                        # Second member of drug pair has larger value than second:  second drug is *overexpressed*
                        geneComps[ii,drugCompCtr] = -1

                # Advance counter for drug pair comparison
                drugCompCtr += 1

    # Iterate thru each drug pair and record list of genes that were marked as expressing different amounts
    alpha = 0.05
    geneCompSig = list()
    geneCompSigRawPval = list()
    for ii in range(NdrugComps):
        curGeneComps = geneComps[:,ii]
        geneCompSig.append(sorted([ geneList[jj] for jj,xx in enumerate(curGeneComps) if xx != 0 ]))
        geneCompSigRawPval.append(sorted([ geneCompsRawPval[jj,ii] for jj,xx in enumerate(curGeneComps) if xx != 0 ]))

    ###
    ### START OF FDR ANALYSIS
    ###

    # Perform statistical correction of p-values to ascertain genes with actual significance for each drug comparison
    geneCompsFdr = np.zeros((Ngenes,NdrugComps))
    geneCompsPvalsFdr = np.zeros((Ngenes,NdrugComps))
    for ii in range(NdrugComps):
        # Get pvalues from current drug comp
        curPvals = geneCompsRawPval[:,ii]

        # Get FDR adjusted list of significant genes in each drug comparison
        curFdrSigs, adjPvals = PerformFdrThresh(curPvals, alpha)

        # Assign polarity to reported significant fdr-adjusted genes, by using polarity
        # of tstat value
        for jj in range(Ngenes):
            if curFdrSigs[jj] != 0:
                if geneCompsRawTstat[jj,ii] < 0:
                    geneCompsFdr[jj, ii] = 1
                else:
                    geneCompsFdr[jj, ii] = -1

        # Recod adjusted pvalues using FDR
        geneCompsPvalsFdr[:,ii] = adjPvals

    geneCompFdrSig = list()
    for ii in range(NdrugComps):
        curGeneComps = geneCompsFdr[:,ii]
        geneCompFdrSig.append(sorted([ geneList[jj] for jj,xx in enumerate(curGeneComps) if xx != 0 ]))

    ###
    ### END OF FDR ANALYSIS
    ###

    # Create new output file in which the rows are different drug pair comparisons and the columns
    # are genes, where the columns are ALL the genes that are found significant across all comparisons
    # This permits easier comparison of up/down regulated genes across different comparisons;
    #
    # Since the header row for this file is the list of total genes found up/down regulated,
    # the entries in the cells are either UP or DOWN, indicating whether the SECOND drug in the drug
    # pair is up or down regulated; if a gene is neither over/under expressed, the entry for that cell is empty

    # Incorporate over/under expression info and organize genes for easy comparison
    geneCompSigAll = sorted(list(set([ yy for xx in geneCompSig for yy in xx ])))

    # Create variable to hold string output to be written to file
    geneCompSigAllOutput = list()
    geneCompSigAllOutput.append([''] + geneCompSigAll)

    # Create matching variable to geneCompSigAllOutput to hold p-values
    geneCompSigPvalsAllOutput = list()
    geneCompSigPvalsAllOutput.append([''] + geneCompSigAll)

    # Iterate thru each drug comparison
    for ii in range(NdrugComps):
        # Make new list for outputs for current drug comp
        geneCompSigAllOutput.append([])
        geneCompSigPvalsAllOutput.append([])
        # Append name of drug comp
        geneCompSigAllOutput[ii+1].append(geneCompStrs[ii])
        geneCompSigPvalsAllOutput[ii+1].append(geneCompStrs[ii])
        # Iterate thru each gene in geneCompSigAll and append appropriate string
        for jj in range(len(geneCompSigAll)):
            # print ii,jj
            # Get name of current gene
            curGene = geneCompSigAll[jj]
            # Check if curGene is in list of significant genes for current drug pair comparison
            if curGene in geneCompSig[ii]:
                # Get index of curGene in geneList
                curGeneInd = geneList.index(curGene)
                if geneComps[curGeneInd,ii] == 1:
                    geneCompSigAllOutput[ii+1].append('DOWN')
                elif geneComps[curGeneInd,ii] == -1:
                    geneCompSigAllOutput[ii+1].append('UP')
                curGeneInd = [ kk for kk,xx in enumerate(geneCompSig[ii]) if curGene == xx ][0]
                geneCompSigPvalsAllOutput[ii+1].append(str(geneCompSigRawPval[ii][curGeneInd]))
            else:
                geneCompSigAllOutput[ii+1].append('')
                geneCompSigPvalsAllOutput[ii+1].append('')

    # Repeat analysis for FDR adjusted tests
    # Create variable to hold string output to be written to file
    geneCompSigAllFdrOutput = list()
    geneCompSigAllFdrOutput.append([''] + geneCompSigAll)

    # Create matching variable to geneCompSigAllOutput to hold p-values
    geneCompSigPvalsAllFdrOutput = list()
    geneCompSigPvalsAllFdrOutput.append([''] + geneCompSigAll)

    # Iterate thru each drug comparison
    for ii in range(NdrugComps):
        # Make new list for outputs for current drug comp
        geneCompSigAllFdrOutput.append([])
        geneCompSigPvalsAllFdrOutput.append([])
        # Append name of drug comp
        geneCompSigAllFdrOutput[ii+1].append(geneCompStrs[ii])
        # Iterate thru each gene in geneCompSigAll and append appropriate string
        for jj in range(len(geneCompSigAll)):
            # print ii,jj
            # Get name of current gene
            curGene = geneCompSigAll[jj]
            # Check if curGene is in list of significant genes for current drug pair comparison
            if curGene in geneCompFdrSig[ii]:
                # Get index of curGene in geneList
                curGeneInd = geneList.index(curGene)
                # curGeneInd = [ kk for kk,xx in enumerate(geneCompSig[ii]) if curGene == xx ][0]
                if geneCompsFdr[curGeneInd,ii] == 1:
                    # print ii,jj, curGene, 'down'
                    geneCompSigAllFdrOutput[ii+1].append('DOWN')
                elif geneCompsFdr[curGeneInd,ii] == -1:
                    # print ii,jj, curGeneComps, 'up'
                    geneCompSigAllFdrOutput[ii+1].append('UP')
                curGeneInd = [ kk for kk,xx in enumerate(geneCompFdrSig[ii]) if curGene == xx ][0]
                geneCompSigPvalsAllFdrOutput[ii+1].append(str(geneCompsPvalsFdr[curGeneInd][ii]))
            else:
                geneCompSigAllFdrOutput[ii+1].append('')
                geneCompSigPvalsAllFdrOutput[ii+1].append('')

    # Write out genes with differential expression to .csv file
    with open(opFileGenes,'wb') as csvfile:
        # Open csv file for holding differentially expressed genes
        csvWriter = csv.writer(csvfile)

        # Write out info about results
        csvWriter.writerow(['Rows are comparisons of two drugs and columns are genes'])
        csvWriter.writerow(['Only genes that were significant in at least one drug comparison are shown'])
        csvWriter.writerow(['UP/DOWN designation refers to direction of expression of the SECOND drug with respect to the first'])

        # Write out description of basic analysis
        csvWriter.writerow([''])
        csvWriter.writerow(['Traditional significance testing using p = 0.05 threshold and no Family-wise correction'])

        # Iterate thru each drug comparison and write out output for each, as formatted above
        for ii in range(NdrugComps):
            csvWriter.writerow(geneCompSigAllOutput[ii])

        # Write out results for identical analysis above except using B-H corrected significance testing
        csvWriter.writerow('')
        csvWriter.writerow(['The same genes identified for the traditional significance testing (above) are used here as column headings'])
        csvWriter.writerow(['False discovery rate (FDR) correction in which the Benjamini-Hochberg correction technique is employed'])
        csvWriter.writerow('')

        # Iterate thru each drug comparison
        for ii in range(NdrugComps):
            csvWriter.writerow(geneCompSigAllFdrOutput[ii])

        # Close .csv file
        csvfile.close()

    # Write out identical analysis as above, except write out p values
    with open(opFilePvals,'wb') as csvfile:
        # Open csv file for holding differentially expressed genes
        csvWriter = csv.writer(csvfile)

        # Write out info about results
        csvWriter.writerow(['Rows are comparisons of two drugs and columns are genes'])
        csvWriter.writerow(['Only genes that were significant in at least one drug comparison are shown'])
        csvWriter.writerow(['UP/DOWN designation refers to direction of expression of the SECOND drug with respect to the first'])

        # Write out description of basic analysis
        csvWriter.writerow([''])
        csvWriter.writerow(['Traditional significance testing using p = 0.05 threshold and no Family-wise correction'])

        # Iterate thru each drug comparison and write out output for each, as formatted above
        for ii in range(NdrugComps):
            csvWriter.writerow(geneCompSigPvalsAllOutput[ii])

        # Write out results for identical analysis above except using B-H corrected significance testing
        csvWriter.writerow('')
        csvWriter.writerow(['The same genes identified for the traditional significance testing (above) are used here as column headings'])
        csvWriter.writerow(['False discovery rate (FDR) correction in which the Benjamini-Hochberg correction technique is employed'])
        csvWriter.writerow('')

        # Iterate thru each drug comparison
        for ii in range(NdrugComps):
            csvWriter.writerow(geneCompSigPvalsAllFdrOutput[ii])

        # Close .csv file
        csvfile.close()



    return

def PerformFdrThresh(pvals,fdrRate = 0.05):
    """

        Method to perform B-H corrections on a set of p-values, with a specified false discovery rate (fdrRate)

    :param pvals:   list of p values from independent experiments
    :param alpha:
    :return:
    """

    import numpy as np

    # Sort pvalues in x
    sortedPvals = sorted(pvals)
    indsSortedPvals = np.argsort(pvals)

    # Get adjusted threshold
    adjThresh = fdrRate/float(len(pvals))

    # Get list of adjusted significances
    adjPvals = [ zz/float(ii+1) for ii,zz in enumerate(sortedPvals) ]

    # Get significant test indices in x
    sigInds = [ ii for ii,zz in enumerate(adjPvals) if zz <= adjThresh ]

    # Take largest value of sigInds and that becomes the index at which and
    # below, the tests are significant

    sigList = np.zeros(len(pvals))
    if sigInds != []:
        maxSigInds = sigInds[-1]
        for jj in range(maxSigInds):
            sigList[indsSortedPvals[jj]] = 1

    # Reorder adjPvals to match order in pvals
    unsortedAdjPvals = [ adjPvals[indsSortedPvals[ii]] for ii,xx in enumerate(range(len(pvals)))  ]

    return sigList, unsortedAdjPvals

    x=5









main()








#