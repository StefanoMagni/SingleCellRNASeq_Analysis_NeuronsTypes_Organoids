#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Created on Mon 21 Jan 2019
# @authors: stefano.magni

from Script_data import *
from Script_information import *
from genes_1 import *
from scipy.stats.stats import pearsonr
import Script_data as sd
import numpy
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
import matplotlib.ticker as plticker
import seaborn as sns
import copy
import matplotlib as mpl
import math


print("Step 0")
# Tags to decide what to run
# Pre Processing
GetDataFromFiles = True                                                             # PRE-PROCESSING
ComputeKneePlots = True                                                             # PRE-PROCESSING
PlotKneePlots = True                                                                # PRE-PROCESSING
CutCellNumber = True                                                                # PRE-PROCESSING          
Normalize = True                                                                    # PRE-PROCESSING
ComparePreliminaries = True                                                         # PRE-PROCESSING

ProduceListOfDifferentiallyExpressedGenes = True
DifferentiallyExpressedGenesIntersectionWithListCellTypesFunctions = True

FoldChangeInterestingGenes = True
FoldChangeInterestingGenesCellTypeProcesses = True

FoldChangeGeneListsNeuronalSubtypes = True
FoldChangeGeneListsGeneralCellTypesProcesses = True

DoComparisonDR = True
DoComparisonDR_ColoredByScores = True
DoComparisonDR_PCA_before_tSNE = True

DoGeneGeneCorrelation_GenesFrom4Lists_cells = True
DoGeneGeneCorrelation_GenesFromCellTypesProcessesLists = True
DoGeneGeneCorrelation_GenesFromCellTypesProcessesLists_OnlyNeuroStem = True

N_cells = 500

L_Gluta_Neuro = sd.get_from_txt('./ListsOfGenes/ListsOfGenesLisaNeuroSubtypes/L_Glutamatergic_Neurons.txt')[0]
L_GABA_Neuro = sd.get_from_txt('./ListsOfGenes/ListsOfGenesLisaNeuroSubtypes/L_GABAergic_Neurons.txt')[0]
L_Seroton_Neuro = sd.get_from_txt('./ListsOfGenes/ListsOfGenesLisaNeuroSubtypes/L_Serotonergic_Neurons.txt')[0]
L_DA_Neuro = sd.get_from_txt('./ListsOfGenes/ListsOfGenesLisaNeuroSubtypes/L_DA_Neurons.txt')[0]

#LisOfListsOfGenesNeuronalSubtypes = [L_Gluta_Neuro, L_GABA_Neuro, L_Seroton_Neuro, L_DA_Neuro]
#
#ListOfListsNames_NeuronalSubtypes = ['List of\nGenes for\nGlutamatergic\nNeurons',
#                        'List of\nGenes for\nGABAergic\nNeurons',
#                        'List of\nGenes for\nSerotonergic\nNeurons',
#                        'List of\nGenes for\nDopaminergic\nNeurons']
#
#ListOfListsNames_NeuronalSubtypes_Short = ['GlutamaNeuro',
#                        'GABANeuro',
#                        'SerotoNeuro',
#                        'DopamNeuro']

LisOfListsOfGenesNeuronalSubtypes = [L_DA_Neuro, L_Gluta_Neuro, L_GABA_Neuro, L_Seroton_Neuro, ]

ListOfListsNames_NeuronalSubtypes = ['List of\nGenes for\nDopaminergic\nNeurons',
                                     'List of\nGenes for\nGlutamatergic\nNeurons',
                                     'List of\nGenes for\nGABAergic\nNeurons',
                                     'List of\nGenes for\nSerotonergic\nNeurons']

ListOfListsNames_NeuronalSubtypes_Short = ['DopamNeuro',
                                           'GlutamaNeuro',
                                           'GABANeuro',
                                           'SerotoNeuro']

L_stem = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_stem.txt')[0] 
L_CellCycle = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_CellCycle.txt')[0]
L_ProApoptotic = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_Pro-Apoptotic.txt')[0]
L_AntiApoptotic = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_Anti-Apoptotic.txt')[0]
L_Caspases = sd.get_from_txt('./ListsOfGenes/ListsOfGenesJonas/list_Caspases.txt')[0]
L_Astro = sd.get_from_txt('./ListsOfGenes/ListsGenesLisa/Lastro2-1_GFAP_S100b_Added.txt')[0]
L_Neurons_Lisa = sd.get_from_txt('./ListsOfGenes/ListsGenesLisa/EnrichedNeuronList.txt')[0]

ListsOfListsOfGenesCellTypesProcesses = [L_stem, L_Neurons_Lisa, L_Astro, L_CellCycle, L_DA_Neuro, L_ProApoptotic, L_AntiApoptotic, L_Caspases]

ListOfListsNames_CellTypesProcesses = ['List of\nGenes for\nStemness',
                        'List of\nGenes for\nNeurons',
                        'List of\nGenes for\nAstrocytes',
                        'List of\nGenes for\nCell Cycle',
                        'List of\nGenes for\nDopaminergic\nNeurons',
                        'List of\nGenes\nPro Apoptotic',
                        'List of\nGenes\nAnti Apoptotic',
                        'List of\nGenes\nCaspases']

ListOfListsNames_CellTypesProcesses_Short = ['Stemness',
                        'Neurons',
                        'Astrocytes',
                        'CellCycle',
                        'DANeurons',
                        'ProApoptotic',
                        'AntiApoptotic',
                        'Caspases']


##### FONTS #####
font05 = {'family': 'serif',
    'color':  'darkred',
    'weight': 'normal',
    'size': 10,
    }
                
font01 = {'family': 'serif',
    'color':  'orange',
    'weight': 'normal',
    'size': 10,
    }
                
font001 = {'family': 'serif',
    'color':  'green',
    'weight': 'normal',
    'size': 10,
    }

font0001 = {'family': 'serif',
    'color':  'darkgreen',
    'weight': 'normal',
    'size': 10,
    }

font2 = {'family': 'serif',
'color':  'blue',
'weight': 'normal',
'size': 15,
}

font2s = {'family': 'serif',
'color':  'blue',
'weight': 'normal',
'size': 8,
}

font3 = {'family': 'serif',
    'color':  'black',
    'weight': 'normal',
    'size': 20,
    }

font3s = {'family': 'serif',
    'color':  'black',
    'weight': 'normal',
    'size': 10,
    }

font3ss = {'family': 'serif',
    'color':  'black',
    'weight': 'normal',
    'size': 6,
    }


if GetDataFromFiles == True:
    print("Step 1: GetDataFromFiles")
    # B = Wild Tipe
    B_35_genes, B_35_expr, B_35_expr_np = sd.get_from_txt('./Data_wt/org3wt35_S1_DGE.txt')
    B_70_genes, B_70_expr, B_70_expr_np = sd.get_from_txt('./Data_wt/org1wt70II_S1_DGE.txt')
    genesBlist = [B_35_genes,B_70_genes]
    
    
if ComputeKneePlots == True:
    print("Step 2: ComputeKneePlots")
    wt_a35,wt_b35 = sd.get_knee_plot(B_35_expr_np)
    wt_a70,wt_b70 = sd.get_knee_plot(B_70_expr_np)


if PlotKneePlots == True:
    print("Step 3: PlotKneePlots")
    plt.clf()    
    plt.xlim([0,3000]) 
    plt.xlabel('Cell Number') 
    plt.ylabel('Cumulative fraction of transcripts')
    plt.plot(wt_a35, label='Wt Day 35')  
    plt.plot(wt_a70, label='Wt Day 70')
    plt.legend(loc='center right')
    plt.title('Knee Plot Cumulative')
    plt.savefig('./KneePlotCumulative.pdf')
    
    plt.clf() 
    plt.xlim([-10,3000]) 
    plt.xlabel('Cell Number') 
    plt.ylabel('Number of transcripts')   
    plt.plot(wt_b35, label='Wt Day 35') 
    plt.plot(wt_b70, label='Wt Day 70') 
    plt.legend(loc='center right')
    plt.title('Knee Plot Absolute')
    plt.savefig('./KneePlotAbsolute.pdf')    


if CutCellNumber == True:
    print("Step 4: CutCellNumber")
    B_35 = sd.get_top_cells_matrix(B_35_expr_np,N_cells)[0]
    B_70 = sd.get_top_cells_matrix(B_70_expr_np,N_cells)[0]


if Normalize == True:    
    print("Step 5: Normalize")
    B_35_norm, B_35_log, B_35_sca = sd.get_normalization(B_35)
    B_70_norm, B_70_log, B_70_sca = sd.get_normalization(B_70)
    auxBlist = [B_35_norm,B_70_norm]
    scaBlist = [B_35_sca,B_70_sca]


if ComparePreliminaries == True:
    print("Step 6: ComparePreliminaries")
    C_35_70 = Compare((B_35_norm,B_35_genes),(B_70_norm,B_70_genes)) 
    print("Comp 1")
    aux_35_70 = C_35_70.merge()[0]
    print("Comp 2")        
    C_35_70_norm, C_35_70_log, C_35_70_sca = sd.get_normalization(aux_35_70)


if ProduceListOfDifferentiallyExpressedGenes == True:
    print("Step 7: ProduceListOfDifferentiallyExpressedGenes")

    #Saves list of genes that are different between Day 35 and Day 70
    aux = C_35_70.avg() 
    a = np.array(aux) #some transformations
    a = np.transpose(a) #some transformations
    N_Bonferroni_Correction = len(C_35_70_norm)
    Corrected_Pvalue = 0.05/N_Bonferroni_Correction
    L_min, L_min_stats = C_35_70.avg_list(Corrected_Pvalue, aux[3], a)
    sd.save_array(L_min_stats, './DiffExprGenes_005_Bonf.csv')

    ListOfDifferentiallyExpressed_Gluta_Genes = []
    ListOfDifferentiallyExpressed_GABA_Genes = []
    ListOfDifferentiallyExpressed_Seroton_Genes = []
    ListOfDifferentiallyExpressed_DAneuro_Genes = []

    for GeneName in L_min:
        if GeneName in L_Gluta_Neuro:
            ListOfDifferentiallyExpressed_Gluta_Genes.append(GeneName)
        if GeneName in L_GABA_Neuro:
            ListOfDifferentiallyExpressed_GABA_Genes.append(GeneName)
        if GeneName in L_Seroton_Neuro:
            ListOfDifferentiallyExpressed_Seroton_Genes.append(GeneName)
        if GeneName in L_DA_Neuro:
            ListOfDifferentiallyExpressed_DAneuro_Genes.append(GeneName)

    sd.save_from_list(ListOfDifferentiallyExpressed_Gluta_Genes, './DiffExprGenes_Gluta.txt')
    sd.save_from_list(ListOfDifferentiallyExpressed_GABA_Genes, './DiffExprGenes_GABA.txt')
    sd.save_from_list(ListOfDifferentiallyExpressed_Seroton_Genes, './DiffExprGenes_Seroton.txt')
    sd.save_from_list(ListOfDifferentiallyExpressed_DAneuro_Genes, './DiffExprGenes_DAneuro.txt')

    ListOfDEGs_WithinLists = [ListOfDifferentiallyExpressed_Gluta_Genes, 
                              ListOfDifferentiallyExpressed_GABA_Genes,
                              ListOfDifferentiallyExpressed_Seroton_Genes,
                              ListOfDifferentiallyExpressed_DAneuro_Genes]

    ListOfListsOfLisaGenesActuallyPresentInData = []
    ListGenesC = C_35_70.comm_genes() 
    i = 0
    LisOfListsOfGenesNeuronalSubtypes_NoAbsentGenes = []
    for CurrentList in LisOfListsOfGenesNeuronalSubtypes:
        ListGenesPresentInDataset = []
        for GeneInMyList in CurrentList:
            if GeneInMyList in ListGenesC:
                ListGenesPresentInDataset.append(GeneInMyList)
        LisOfListsOfGenesNeuronalSubtypes_NoAbsentGenes.append(ListGenesPresentInDataset)
        # Here a fastes way to intersect lists, 
        # BUT WITH THE DRAWBACK OF NOT PRESERVING ORDER (due to set()) -.-' So not used.
        # ListGenesPresentInDataset = list(set(ListGenesC).intersection(CurrentList))
        print(' ')
        print(len(CurrentList))
        print(len(ListGenesPresentInDataset))
        print("The following genes are absent from the common data:")
        print([item for item in CurrentList if item not in ListGenesC])
        ListOfListsOfLisaGenesActuallyPresentInData.append(ListGenesPresentInDataset)
        print("The following genes are DEGs within the current list of genes:")
        print(ListOfDEGs_WithinLists[i])
        print(len(ListOfDEGs_WithinLists[i]))
        print(len(ListOfDEGs_WithinLists[i]) / len(ListOfListsOfLisaGenesActuallyPresentInData[i]))
        i = i+1

    print(" ")
    print(len(ListGenesC))
    print(len(L_min))
    print(len(L_min)/len(ListGenesC))
    print(" ")
    print(" ")


if DifferentiallyExpressedGenesIntersectionWithListCellTypesFunctions == True:
    print("Step 7b: DifferentiallyExpressedGenesIntersectionWithListCellTypesFunctions")

    ListOfListsOfGenes = ListsOfListsOfGenesCellTypesProcesses

    ListOfDEGs_WithinLists = []
    i = 0
    for CurrentList in ListOfListsOfGenes:
        ListOfDEGsWithinCurrentList = []
        for GeneName in L_min:
            if GeneName in CurrentList:
                ListOfDEGsWithinCurrentList.append(GeneName)
                GeneListName = ListOfListsNames_CellTypesProcesses_Short[i]
        sd.save_from_list(ListOfDEGsWithinCurrentList, './DEGs' + str(GeneListName) + '.txt')
        ListOfDEGs_WithinLists.append(ListOfDEGsWithinCurrentList)
        i = i + 1

    ListOfListsOfGenesActuallyPresentInData = []
    ListGenesC = C_35_70.comm_genes() 
    i = 0
    for CurrentList in ListOfListsOfGenes:
        ListGenesPresentInDataset = []
        for GeneInMyList in CurrentList:
            if GeneInMyList in ListGenesC:
                ListGenesPresentInDataset.append(GeneInMyList)
        ListOfListsOfGenesActuallyPresentInData.append(ListGenesPresentInDataset)
        print(' ')
        print(len(CurrentList))
        print(len(ListGenesPresentInDataset))
        print("The following genes are absent from the common data:")
        print([item for item in CurrentList if item not in ListGenesC])
        print("The following genes are DEGs within the current list of genes:")
        print(ListOfDEGs_WithinLists[i])
        print(len(ListOfDEGs_WithinLists[i]))
        print(len(ListOfDEGs_WithinLists[i]) / len(ListOfListsOfGenesActuallyPresentInData[i]))
        i = i+1
    ListsOfListsOfGenesCellTypesProcessesActuallyinData = ListOfListsOfGenesActuallyPresentInData

def ComputeFoldChange(GenesExpressionA,GenesExpressionB):
            
            MeanA = sum(GenesExpressionA)/len(GenesExpressionA)
            MeanB = sum(GenesExpressionB)/len(GenesExpressionB)
            
            StdDevA = numpy.std(GenesExpressionA)
            StdDevB = numpy.std(GenesExpressionB)
            
            StandardErrorOfTheMeanA = StdDevA / math.sqrt(len(GenesExpressionA))
            StandardErrorOfTheMeanB = StdDevB / math.sqrt(len(GenesExpressionB))
            
            RatioMutOnCtrl = MeanA/MeanB
            MyBase = 2
            if RatioMutOnCtrl != 0:
                FoldChangeCurrentGene = math.log(RatioMutOnCtrl, MyBase) 
            else:
                FoldChangeCurrentGene = 0
            
            Temp1 = StandardErrorOfTheMeanA ** 2 / (MeanA ** 2 * math.log(2) ** 2)
            Temp2 = StandardErrorOfTheMeanB ** 2 / (MeanB ** 2 * math.log(2) ** 2)
            
            STDDEVofFoldChangeCurrentGene = math.sqrt(Temp1 + Temp2)
            
            return FoldChangeCurrentGene, STDDEVofFoldChangeCurrentGene
     
        
import statsmodels.stats.weightstats as wstats

def ApplyZtestOnMeansOfDistributions(DataPopulation1, DataPopulation2):
    ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
    # of 2 populations are statistically significantly different ######
    DataPopulation1_Object = wstats.DescrStatsW(DataPopulation1)
    DataPopulation2_Object = wstats.DescrStatsW(DataPopulation2)
    MyComparison = wstats.CompareMeans(DataPopulation1_Object,DataPopulation2_Object)
    TestStatistics, pvalue = MyComparison.ztest_ind(alternative='two-sided', usevar='unequal', value=0)

    return pvalue   


def autolabel_forOneComparison(ax,rects, values, ListOfFoldChangesThisDay, GeneStd, LabelUpBy, LabelDownBy):
    # attach pvalue text labels
    i = 0
    for rect in rects:
        height = ListOfFoldChangesThisDay
        Displacement = 0
        if ListOfFoldChangesThisDay > 0:
            Displacement = GeneStd + LabelUpBy
        elif ListOfFoldChangesThisDay < 0:
            Displacement = - GeneStd - LabelDownBy
        ax.text(rect.get_x()+rect.get_width()/2., height + Displacement, values,
                ha='center', va='bottom')
        i=i+1

 
def PlotFoldChangesManyGenesOneComparison(CurrentNamesList, 
                                              ListOfFoldChanges, 
                                              ListOfSTDEVforFoldChanges, 
                                              ListOfPvalues, 
                                              NBonferroni, 
                                              My_ylabel,
                                              LabelUpBy,
                                              LabelDownBy):
    
    N = 1 # Number of days
    ind = np.arange(N)  # the x locations for the groups
    width = 0.15       # the width of the bars
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ListOfRects = []
    i = 0
    ColorsList = ['b','g','r','c','m','y','k','b','g','r','c','m','y','k','b','g','r','c','m','y','k'
                  'b','g','r','c','m','y','k''b','g','r','c','m','y','k''b','g','r','c','m','y','k'
                  'b','g','r','c','m','y','k''b','g','r','c','m','y','k''b','g','r','c','m','y','k'
                  'b','g','r','c','m','y','k''b','g','r','c','m','y','k''b','g','r','c','m','y','k'
                  'b','g','r','c','m','y','k''b','g','r','c','m','y','k''b','g','r','c','m','y','k']
    for IndexGenes in range(len(CurrentNamesList)):
        Gene0FoldChanges = ListOfFoldChanges[IndexGenes]
        Gene0Std = ListOfSTDEVforFoldChanges[IndexGenes]
        ListOfRects.append(ax.bar(ind+i*width, Gene0FoldChanges, width, color=ColorsList[i], yerr=Gene0Std))
        i = i+1
    
    ax.set_ylabel(My_ylabel)
    # ax.set_title('Fold changes for CRC core genes')
    ListOfSteps = []
    for step in range(len(CurrentNamesList)):
        ListOfSteps.append(ind+step*width)
    ax.set_xticks(ListOfSteps)
    ax.set_xticklabels(CurrentNamesList, rotation=90)
    # ax.legend( , ('NR2F1', 'NR2F2', 'POU3F2', 'POU3F3', 'SOX2') )

    plt.show()
            
    ListOfAsterisks = []
    for current_pvalue in ListOfPvalues:
        if current_pvalue < 0.0001 / NBonferroni:
            Asterisks = '****'
        elif current_pvalue < 0.001 / NBonferroni:
            Asterisks = '***'
        elif current_pvalue < 0.01 / NBonferroni:
            Asterisks = '**'
        elif current_pvalue < 0.05 / NBonferroni:
            Asterisks = '*'
        else:
            Asterisks = ' '
        ListOfAsterisks.append(Asterisks)
            
    i = 0
    for i in range(len(CurrentNamesList)):
        print(i)
        autolabel_forOneComparison(ax,ListOfRects[i],ListOfAsterisks[i], ListOfFoldChanges[i], ListOfSTDEVforFoldChanges[i], LabelUpBy[i], LabelDownBy[i])
        i= i + 1

               
if FoldChangeInterestingGenes == 1:
    print("Step 8: FoldChangeInterestingGenes")

    NBonferroni = len(ListOfListsOfLisaGenesActuallyPresentInData[0]) 
    + len(ListOfListsOfLisaGenesActuallyPresentInData[1]) 
    + len(ListOfListsOfLisaGenesActuallyPresentInData[2]) 
    + len(ListOfListsOfLisaGenesActuallyPresentInData[3])
    
    for CurrentList in ListOfListsOfLisaGenesActuallyPresentInData:
        
        ListOfFoldChangesThisDay = []
        ListOfSTDEVforFoldChangesThisDay = []
        ListOfPvalues = []
        for i in range(len(CurrentList)): 
            GeneName = CurrentList[i]

            GenesExpression35 = C_35_70_norm[C_35_70.c_names.index(GeneName),0:N_cells] # !!! THIS IS WT DAY 35 !!!
            GenesExpression70 = C_35_70_norm[C_35_70.c_names.index(GeneName),N_cells:N_cells+N_cells] # !!! THIS IS WT DAY 70 !!!
            
            FoldChangeCurrentGene, STDDEVofFoldChangeCurrentGene = ComputeFoldChange(
                    GenesExpression35,GenesExpression70)
            
            ListOfFoldChangesThisDay.append(FoldChangeCurrentGene)
            ListOfSTDEVforFoldChangesThisDay.append(STDDEVofFoldChangeCurrentGene)
            
            ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
            # of 2 populations are statistically significantly different ######
            pvalue = ApplyZtestOnMeansOfDistributions(GenesExpression35, GenesExpression70)

            # Remember Bonferroni correction for multiple testing, 
            # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
            # by the number of repetitions of the test, 4 days x 7 lists = 28
            # The 7 lists are Stem, DAneur, Mito, Cell Cycle, Pro-apoptosis, Anti-apoptosis, Caspases
            print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
            print(' ')
            ListOfPvalues.append(pvalue)

  
        ########## NOW LET'S PLOT THIS FOLD-CHANGES ##########
        My_ylabel = r'Fold change ($log_2(\langle Day 35 \rangle / \langle Day 70 \rangle)$)'
#        LabelUpBy = [-0.2, 0, -0.2, -0.2]
#        LabelDownBy = [0.5, 0.2, 0.5, 0.5]
        LabelUpBy = -0.2 * np.ones(100)
        LabelDownBy = 0.5 * np.ones(100)
        PlotFoldChangesManyGenesOneComparison(CurrentList, 
                                                  ListOfFoldChangesThisDay, 
                                                  ListOfSTDEVforFoldChangesThisDay, 
                                                  ListOfPvalues, 
                                                  NBonferroni, 
                                                  My_ylabel,
                                                  LabelUpBy,
                                                  LabelDownBy)
        
        
if FoldChangeInterestingGenesCellTypeProcesses == 1:
    print("Step 8b: FoldChangeInterestingGenesCellTypeProcesses")

    NBonferroni = len([GeneName for List in ListsOfListsOfGenesCellTypesProcessesActuallyinData for GeneName in List])
    
    for CurrentList in ListsOfListsOfGenesCellTypesProcessesActuallyinData:
        
        ListOfFoldChangesThisDay = []
        ListOfSTDEVforFoldChangesThisDay = []
        ListOfPvalues = []
        for i in range(len(CurrentList)): 
            GeneName = CurrentList[i]

            GenesExpression35 = C_35_70_norm[C_35_70.c_names.index(GeneName),0:N_cells] # !!! THIS IS WT DAY 35 !!!
            GenesExpression70 = C_35_70_norm[C_35_70.c_names.index(GeneName),N_cells:N_cells+N_cells] # !!! THIS IS WT DAY 70 !!!
            
            FoldChangeCurrentGene, STDDEVofFoldChangeCurrentGene = ComputeFoldChange(
                    GenesExpression35,GenesExpression70)
            
            ListOfFoldChangesThisDay.append(FoldChangeCurrentGene)
            ListOfSTDEVforFoldChangesThisDay.append(STDDEVofFoldChangeCurrentGene)
            
            ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
            # of 2 populations are statistically significantly different ######
            pvalue = ApplyZtestOnMeansOfDistributions(GenesExpression35, GenesExpression70)

            # Remember Bonferroni correction for multiple testing, 
            # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
            # by the number of repetitions of the test, 4 days x 7 lists = 28
            # The 7 lists are Stem, DAneur, Mito, Cell Cycle, Pro-apoptosis, Anti-apoptosis, Caspases
            print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
            print(' ')
            ListOfPvalues.append(pvalue)

  
        ########## NOW LET'S PLOT THIS FOLD-CHANGES ##########
        My_ylabel = r'Fold change ($log_2(\langle Day 35 \rangle / \langle Day 70 \rangle)$)'
        LabelUpBy = -0.2 * np.ones(100)
        LabelDownBy = 0.5 * np.ones(100)
        PlotFoldChangesManyGenesOneComparison(CurrentList, 
                                                  ListOfFoldChangesThisDay, 
                                                  ListOfSTDEVforFoldChangesThisDay, 
                                                  ListOfPvalues, 
                                                  NBonferroni, 
                                                  My_ylabel,
                                                  LabelUpBy,
                                                  LabelDownBy)


def ComputeCumulativeGeneExpressionFromListAndCompareObject(C, C_norm, N_cells, CurrentGenesList):
    
    CumulativeGeneExpressionA = np.zeros(N_cells)
    CumulativeGeneExpressionB = np.zeros(N_cells)
    
    for i in range(len(CurrentGenesList)): 
        GeneName = CurrentGenesList[i]
        if GeneName in C.c_names:
            CumulativeGeneExpressionA = CumulativeGeneExpressionA + C_norm[C.c_names.index(GeneName),0:N_cells]
            CumulativeGeneExpressionB = CumulativeGeneExpressionB + C_norm[C.c_names.index(GeneName),N_cells:N_cells+N_cells]
        elif GeneName not in C.c_names:
            print('Gene ' + GeneName + ' was not found in compare object C, when computing Cumulative Gene Expression (on list with ' + str(len(CurrentGenesList)) + ' genes)')
            
    return CumulativeGeneExpressionA, CumulativeGeneExpressionB
        

def PlotHistogramsCumulatives(Nlines, Ncols, ListOfCumulatives35, ListOfCumulatives70, ListOfListsNames, N_cells, NBonferroni, ylimVector=[], FontListLabel=font2):

    fig, axes = plt.subplots(nrows=Nlines, ncols=Ncols)
    fig.subplots_adjust(hspace=0.3, wspace=0.05)

    for ax in axes.flat:
        # Hide ticks and labels on X
        ax.xaxis.set_visible(True)
        ax.yaxis.set_visible(True)

    # Set up ticks only on one side for the "edge" subplots...
    if ax.is_first_col():
        ax.yaxis.set_ticks_position('left')
        ax.yaxis.set_visible(True)
    
    for i in range(Nlines):            
        MaxX = max(max(ListOfCumulatives35[i]),max(ListOfCumulatives70[i]))
        CumulativeGeneExpression35_NORM = ListOfCumulatives35[i] / MaxX
        CumulativeGeneExpression70_NORM = ListOfCumulatives70[i] / MaxX
    
        binwidth = 0.05 # 0.010 # 0.025
        axes[i].hist(CumulativeGeneExpression35_NORM,bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color='Red', normed=0, label='35d DA dif', alpha=0.5)        
        axes[i].hist(CumulativeGeneExpression70_NORM,bins=np.arange(0, 1 + binwidth, binwidth), histtype='stepfilled', color='Blue', normed=0, label='70d DA dif', alpha=0.5)
            
        axes[i].set_ylabel("N Cells")
        axes[Nlines-1].set_xlabel("Cumulative Gene Expression (norm. to max)")
        # axes[i].set_ylim(0,220)
        # axes[2].set_ylim(0,450)
        if ylimVector: # this returns True if the vecor is not empty
            axes[i].set_ylim(0,ylimVector[i])
        axes[i].text(1.1, 0, ListOfListsNames[i], fontdict=FontListLabel) 

        if i == Nlines-1:
            plt.legend()
                            
        BONF = 'Cumulative Gene Expressions, each distribution has '+ str(N_cells) + ' cells,\n Bonferroni correction applied for ' + str(round(NBonferroni)) + ' z-tests repetitions'
        plt.suptitle(BONF)
        
        if ListOfPvalues[i] < 0.0001 / NBonferroni:
            PValueLabel = 'p-val < 0.0001'
            font = font0001
        elif ListOfPvalues[i] < 0.001 / NBonferroni:
            PValueLabel = 'p-val < 0.001'
            font = font001
        elif ListOfPvalues[i] < 0.01 / NBonferroni:
            PValueLabel = 'p-val < 0.01'
            font = font01
        elif ListOfPvalues[i] < 0.05 / NBonferroni:
            PValueLabel = 'p-val < 0.05'
            font = font05
        else:
            PValueLabel = ' '
            font = font05
            
        axes[i].text(0.65, 50, PValueLabel, fontdict=font)

    plt.show()
    
    
if FoldChangeGeneListsNeuronalSubtypes == 1:
    print("Step 9: FoldChangeGeneListsNeuronalSubtypes")

    ListOfListsOfGenesForCumulatives = ListOfListsOfLisaGenesActuallyPresentInData

    NBonferroni = 6 # len(ListOfListsOfGenesForCumulatives)
    # This number 6 is only to account for the lists which actually go into the paper. 
    # For a more general use, keep len(...)
    
    ListOfFoldChangesThisDay = []
    ListOfSTDEVforFoldChangesThisDay = []
    ListOfPvalues = []
    
    ListOfCumulatives35 = []
    ListOfCumulatives70 = []

    for CurrentList in ListOfListsOfGenesForCumulatives:

        CumulativeGeneExpression35, CumulativeGeneExpression70 = ComputeCumulativeGeneExpressionFromListAndCompareObject(
                C_35_70, C_35_70_norm, N_cells, CurrentList)
        
        ListOfCumulatives35.append(CumulativeGeneExpression35)
        ListOfCumulatives70.append(CumulativeGeneExpression70)
        
        FoldChangeCurrentGene, STDDEVofFoldChangeCurrentGene = ComputeFoldChange(
                CumulativeGeneExpression35,CumulativeGeneExpression70)
        
        ListOfFoldChangesThisDay.append(FoldChangeCurrentGene)
        ListOfSTDEVforFoldChangesThisDay.append(STDDEVofFoldChangeCurrentGene)
        
        ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
        # of 2 populations are statistically significantly different ######
        pvalue = ApplyZtestOnMeansOfDistributions(CumulativeGeneExpression35, CumulativeGeneExpression70)

        # Remember Bonferroni correction for multiple testing, 
        # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
        # by the number of repetitions of the test, 4 days x 7 lists = 28
        # The 7 lists are Stem, DAneur, Mito, Cell Cycle, Pro-apoptosis, Anti-apoptosis, Caspases
        print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
        print(' ')
        ListOfPvalues.append(pvalue)
            
    ########## NOW LET'S PLOT THIS FOLD-CHANGES ##########
    My_ylabel = r'Fold change ($log_2(\langle Day 35 \rangle / \langle Day 70 \rangle)$)'
    ListOfListsNames = ListOfListsNames_NeuronalSubtypes
    LabelUpBy = [0,0,0,0]
    LabelDownBy = [0.2,0.2,0.2,0.2]
    PlotFoldChangesManyGenesOneComparison(ListOfListsNames, 
                                              ListOfFoldChangesThisDay, 
                                              ListOfSTDEVforFoldChangesThisDay, 
                                              ListOfPvalues, 
                                              NBonferroni, 
                                              My_ylabel,
                                              LabelUpBy,
                                              LabelDownBy)
    
    ########## MAKE HISTOGRAMS ##########            
    Nlines = 4 
    Ncols = 1
    ylimVector = [220, 220, 220, 450]
    PlotHistogramsCumulatives(Nlines, Ncols, ListOfCumulatives35, ListOfCumulatives70, ListOfListsNames, N_cells, NBonferroni, ylimVector)
  
   
if FoldChangeGeneListsGeneralCellTypesProcesses == True:
    print("Step 12: FoldChangeGeneListsGeneralCellTypesProcesses")

    ListOfListsOfGenesForCumulatives = ListsOfListsOfGenesCellTypesProcesses

    NBonferroni = 6 # len(ListOfListsOfGenesForCumulatives)
    # This number 6 is only to account for the lists which actually go into the paper. 
    # For a more general use, keep len(...)
    
    ListOfFoldChangesThisDay = []
    ListOfSTDEVforFoldChangesThisDay = []
    ListOfPvalues = []
    
    ListOfCumulatives35 = []
    ListOfCumulatives70 = []

    for CurrentList in ListOfListsOfGenesForCumulatives:

        CumulativeGeneExpression35, CumulativeGeneExpression70 = ComputeCumulativeGeneExpressionFromListAndCompareObject(
                C_35_70, C_35_70_norm, N_cells, CurrentList)
        
        ListOfCumulatives35.append(CumulativeGeneExpression35)
        ListOfCumulatives70.append(CumulativeGeneExpression70)
        
        FoldChangeCurrentGene, STDDEVofFoldChangeCurrentGene = ComputeFoldChange(
                CumulativeGeneExpression35,CumulativeGeneExpression70)
        
        ListOfFoldChangesThisDay.append(FoldChangeCurrentGene)
        ListOfSTDEVforFoldChangesThisDay.append(STDDEVofFoldChangeCurrentGene)
        
        ########## APPLY Z-TEST (FOR DIFFERENT VARIANCES to verify if means 
        # of 2 populations are statistically significantly different ######
        pvalue = ApplyZtestOnMeansOfDistributions(CumulativeGeneExpression35, CumulativeGeneExpression70)

        # Remember Bonferroni correction for multiple testing, 
        # computed by dividing the thrashold you wish, e.g. 0.01 pval, 
        # by the number of repetitions of the test, 4 days x 7 lists = 28
        # The 7 lists are Stem, DAneur, Mito, Cell Cycle, Pro-apoptosis, Anti-apoptosis, Caspases
        print('p-value = ' + str(pvalue) + '   (Remember to use Bonferroni Correction)')
        print(' ')
        ListOfPvalues.append(pvalue)
            
    ########## NOW LET'S PLOT THIS FOLD-CHANGES ##########
    My_ylabel = r'Fold change ($log_2(\langle Day 35 \rangle / \langle Day 70 \rangle)$)'
    ListOfListsNames = ListOfListsNames_CellTypesProcesses
    LabelUpBy = [0,0,0,0,0,0,0,0]
    LabelDownBy = [0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2]
    PlotFoldChangesManyGenesOneComparison(ListOfListsNames, 
                                              ListOfFoldChangesThisDay, 
                                              ListOfSTDEVforFoldChangesThisDay, 
                                              ListOfPvalues, 
                                              NBonferroni, 
                                              My_ylabel,
                                              LabelUpBy,
                                              LabelDownBy)
    
    ########## MAKE HISTOGRAMS ##########            
    Nlines = 8 
    Ncols = 1
    PlotHistogramsCumulatives(Nlines, Ncols, ListOfCumulatives35, ListOfCumulatives70, ListOfListsNames, N_cells, NBonferroni, FontListLabel=font2s) 
    
    Nlines = 2 
    NBonferroni = Nlines
    PlotHistogramsCumulatives(Nlines, Ncols, ListOfCumulatives35[0:Nlines], ListOfCumulatives70[0:Nlines], ListOfListsNames[0:Nlines], N_cells, NBonferroni, FontListLabel=font2s) 
      
    
if DoComparisonDR == 1:
    print("Step 10a: DoComparisonDR")

    Myt1 = 'COMPARISON,sca - PCA'
    Myt2 = 'COMPARISON,sca - tSNE'
    MySaveAs1 = './COMPARISON_PCA_sca'+ '_' + str(N_cells) + '.pdf'
    MySaveAs2 = './COMPARISON_tSNE_sca'+ '_' + str(N_cells) + '.pdf'
    a,b = get_DR_bis(C_35_70_sca[:,0:N_cells],C_35_70_sca[:,N_cells:N_cells+N_cells],t1=Myt1,t2=Myt2,SaveAs1=MySaveAs1,SaveAs2=MySaveAs2)

    Myt1 = 'COMPARISON,norm - PCA'
    Myt2 = 'COMPARISON,norm - tSNE'
    MySaveAs1 = './COMPARISON_PCA_norm'+ '_' + str(N_cells) + '.pdf'
    MySaveAs2 = './COMPARISON_tSNE_norm'+ '_' + str(N_cells) + '.pdf'
    a,b = get_DR_bis(C_35_70_norm[:,0:N_cells],C_35_70_norm[:,N_cells:N_cells+N_cells],t1=Myt1,t2=Myt2,SaveAs1=MySaveAs1,SaveAs2=MySaveAs2)

#    Myt1 = 'COMPARISON,norm - PCA 3D'
#    Myt2 = 'COMPARISON,norm - tSNE 3D'
#    MySaveAs1 = './COMPARISON_PCA_3D_norm'+ '_' + str(N_cells) + '.pdf'
#    MySaveAs2 = './COMPARISON_tSNE_3D_norm'+ '_' + str(N_cells) + '.pdf'
#    a,b = get_DR_bis(C_35_70_norm[:,0:N_cells],C_35_70_norm[:,N_cells:N_cells+N_cells],com=3,t1=Myt1,t2=Myt2,SaveAs1=MySaveAs1,SaveAs2=MySaveAs2)

if DoComparisonDR_ColoredByScores == 1:
    print("Step 10b: DoComparisonDR_ColoredByScores")
    
    k = 0
    for ListOfListsOfGenes in [LisOfListsOfGenesNeuronalSubtypes, ListsOfListsOfGenesCellTypesProcesses]:
    
        ListOfListsOfGenesForCumulatives = ListOfListsOfGenes
        ListOfCumulatives35 = []
        ListOfCumulatives70 = []
        
        if k == 0:
            ListNames_Short = ListOfListsNames_NeuronalSubtypes_Short
        elif k == 1:
            ListNames_Short = ListOfListsNames_CellTypesProcesses_Short
    
        for CurrentList in ListOfListsOfGenesForCumulatives:
    
            CumulativeGeneExpression35, CumulativeGeneExpression70 = ComputeCumulativeGeneExpressionFromListAndCompareObject(
                    C_35_70, C_35_70_norm, N_cells, CurrentList)
            
            ListOfCumulatives35.append(CumulativeGeneExpression35)
            ListOfCumulatives70.append(CumulativeGeneExpression70)
        
        i = 0
        for i in range(len(ListOfListsOfGenesForCumulatives)):
            #GeneNameToColor = 'CD4'
            Myt1 = 'COMPARISON,norm - PCA\npoints colored by expression of ' + str(ListNames_Short[i])
            Myt2 = 'COMPARISON,norm - tSNE\npoints colored by expression of ' + str(ListNames_Short[i])
            MySaveAs1 = './COMPARISON_PCA_norm_' + str(ListNames_Short[i]) + '_' + str(N_cells) + '.pdf'
            MySaveAs2 = './COMPARISON_tSNE_norm_' + str(ListNames_Short[i]) + '_' + str(N_cells) + '.pdf'
            ColorMap = 'rainbow'
            CumulativeGeneExpressionOneList = np.concatenate((ListOfCumulatives35[i], ListOfCumulatives70[i]), axis=0)
            a,b = get_DR_Colored_By_User_Provided_Vector(C_35_70_norm,CumulativeGeneExpressionOneList,ColorMap,2,[],Myt1,Myt2,MySaveAs1,MySaveAs2)

        k = k + 1

    ListOfCumulatives35 = []
    ListOfCumulatives70 = []
    
#    for GeneNameToColor in ['FOXA2','TH']:
#        #GeneNameToColor = 'CD4'
#        Myt1 = 'COMPARISON,norm - PCA\npoints colored by expression of ' + str(GeneNameToColor)
#        Myt2 = 'COMPARISON,norm - tSNE\npoints colored by expression of ' + str(GeneNameToColor)
#        MySaveAs1 = './COMPARISON_PCA_norm_' + str(GeneNameToColor) + '_' + str(N_cells) + '.pdf'
#        MySaveAs2 = './COMPARISON_tSNE_norm_' + str(GeneNameToColor) + '_' + str(N_cells) + '.pdf'
#        ColorMap = 'rainbow'
#        a,b = get_DR_Colored_By_User_Provided_Vector(C_35_70_norm,C_35_70_norm[C_35_70.c_names.index(GeneNameToColor),:],ColorMap,2,[],Myt1,Myt2,MySaveAs1,MySaveAs2)


if DoComparisonDR_PCA_before_tSNE == 1:
    print("Step 10c: DoComparisonDR_PCA_before_tSNE")
    
    NumberOfPCAbeforetSNE = 10 #20 #10 #6
    
    y_pca_temp, y_tsne_temp = get_DR(C_35_70_norm,g_ind=-1, d=NumberOfPCAbeforetSNE, off=[])
    y_pca_temp = np.transpose(y_pca_temp)
    Myt1 = 'COMPARISON,norm - PCA_on_' + str(NumberOfPCAbeforetSNE) + '_PCA_'
    Myt2 = 'COMPARISON,norm - tSNE_on_' + str(NumberOfPCAbeforetSNE) + '_PCA_'
    MySaveAs1 = './COMPARISON_PCA_on_' + str(NumberOfPCAbeforetSNE) + '_PCA_on_norm'+ '_' + str(N_cells) + '.pdf'
    MySaveAs2 = './COMPARISON_tSNE_on_' + str(NumberOfPCAbeforetSNE) + '_PCA_on_norm'+ '_' + str(N_cells) + '.pdf'
    a,b = get_DR_bis(y_pca_temp[:,0:N_cells],y_pca_temp[:,N_cells:N_cells+N_cells],com=3,t1=Myt1,t2=Myt2,SaveAs1=MySaveAs1,SaveAs2=MySaveAs2)

    k = 0
    for ListOfListsOfGenes in [LisOfListsOfGenesNeuronalSubtypes, ListsOfListsOfGenesCellTypesProcesses]:
    
        ListOfListsOfGenesForCumulatives = ListOfListsOfGenes
        ListOfCumulatives35 = []
        ListOfCumulatives70 = []
        
        if k == 0:
            ListNames_Short = ListOfListsNames_NeuronalSubtypes_Short
        elif k == 1:
            ListNames_Short = ListOfListsNames_CellTypesProcesses_Short
    
        for CurrentList in ListOfListsOfGenesForCumulatives:
    
            CumulativeGeneExpression35, CumulativeGeneExpression70 = ComputeCumulativeGeneExpressionFromListAndCompareObject(
                    C_35_70, C_35_70_norm, N_cells, CurrentList)
            
            ListOfCumulatives35.append(CumulativeGeneExpression35)
            ListOfCumulatives70.append(CumulativeGeneExpression70)
        
        i = 0
        for i in range(len(ListOfListsOfGenesForCumulatives)):
            #GeneNameToColor = 'CD4'
            Myt1 = 'COMPARISON,norm - PCA_on_' + str(NumberOfPCAbeforetSNE) + '_PCA_\npoints colored by expression of ' + str(ListNames_Short[i])
            Myt2 = 'COMPARISON,norm - tSNE_on_' + str(NumberOfPCAbeforetSNE) + '_PCA_\npoints colored by expression of ' + str(ListNames_Short[i])
            MySaveAs1 = './COMPARISON_PCA_on_' + str(NumberOfPCAbeforetSNE) + '_PCA_norm_' + str(ListNames_Short[i]) + '_' + str(N_cells) + '.pdf'
            MySaveAs2 = './COMPARISON_tSNE_on_' + str(NumberOfPCAbeforetSNE) + '_PCA_norm_' + str(ListNames_Short[i]) + '_' + str(N_cells) + '.pdf'
            ColorMap = 'rainbow'
            CumulativeGeneExpressionOneList = np.concatenate((ListOfCumulatives35[i], ListOfCumulatives70[i]), axis=0)
            a,b = get_DR_Colored_By_User_Provided_Vector(y_pca_temp,CumulativeGeneExpressionOneList,ColorMap,3,[],Myt1,Myt2,MySaveAs1,MySaveAs2)

        k = k + 1

    ListOfCumulatives35 = []
    ListOfCumulatives70 = []
    
def ComputeGeneGeneCorrelationOn2SamplesOfCompareObject(C,
                                                        C_norm,
                                                        List_Gene_Names_for_Correlation,
                                                        N_cells,
                                                        ThrasoldForPval,
                                                        RemoveNotSignificantCorr,
                                                        WhiteDiagonal):

    CurrentList = []
    N = len(List_Gene_Names_for_Correlation)
    pvalPearsonCorr_Matrix = 100 * np.ones([N, N])
    MyCorrelationMatrix_Full = 100 * np.ones([N, N])
    MyCorrelationMatrix_Masked = 100 * np.ones([N, N])
    
    i = 0
    for GeneName1 in List_Gene_Names_for_Correlation:
        j = 0
        CurrentList.append(GeneName1)
        for GeneName2 in List_Gene_Names_for_Correlation:
            
            Gene_1_Expression = C_norm[C.c_names.index(GeneName1),0+SampleIndex*N_cells:N_cells+SampleIndex*N_cells]
            Gene_2_Expression = C_norm[C.c_names.index(GeneName2),0+SampleIndex*N_cells:N_cells+SampleIndex*N_cells]
            
            PearsonCorrCoeff, pvalPearsonCorr = pearsonr(Gene_1_Expression, Gene_2_Expression)
            pvalPearsonCorr_Matrix[i,j] = pvalPearsonCorr
            
#            if RemoveNotSignificantCorr == False:
#                MyCorrelationMatrix_Full[i,j] = PearsonCorrCoeff
#                              
#                if WhiteDiagonal == True:
#                    MyCorrelationMatrix_Masked[i,j] = PearsonCorrCoeff
#                    if i == j:
#                        MyCorrelationMatrix_Masked[i,j] = 0
                        
#            elif RemoveNotSignificantCorr == True:
            MyCorrelationMatrix_Full[i,j] = PearsonCorrCoeff

            if pvalPearsonCorr <= ThrasoldForPval:
                MyCorrelationMatrix_Masked[i,j] = PearsonCorrCoeff
            elif pvalPearsonCorr > ThrasoldForPval:
                MyCorrelationMatrix_Masked[i,j] = 0
                          
            if WhiteDiagonal == True:
#                MyCorrelationMatrix_Masked[i,j] = MyCorrelationMatrix_Full[i,j]
                if i == j:
                    MyCorrelationMatrix_Masked[i,j] = 0 
                    MyCorrelationMatrix_Full[i,j] = 0 
            j = j + 1
        i = i + 1
    
    return MyCorrelationMatrix_Full, pvalPearsonCorr_Matrix, MyCorrelationMatrix_Masked

    
    return MyCorrelationMatrix_Full, pvalPearsonCorr_Matrix, MyCorrelationMatrix_Masked


if DoGeneGeneCorrelation_GenesFrom4Lists_cells == 1:
    print("Step 11: DoGeneGeneCorrelation_GenesFrom4Lists_cells")
        
    WhiteDiagonal = True # True --> mask values (1s) from the diagonal, False --> do not mask
    RemoveNotSignificantCorr = True
    RemoveNANs = True
    
    # Let's put the 4 lists of genes typical from neuronal subtypes from Lisa all together, making up a total of 60 genes 
    List_Gene_Names_for_Correlation = [gene for GenesList in LisOfListsOfGenesNeuronalSubtypes_NoAbsentGenes for gene in GenesList]
    
    Nsamples = 2  # Here 2 samples, in fact the samples in C object are 0 = Day 35, 1 = Day 70
    Ngenes = len(List_Gene_Names_for_Correlation)
    N_Bonf_Corr = (Ngenes*Ngenes-Ngenes)/2 * Nsamples 
    ThrasoldForPval = 0.05/N_Bonf_Corr
        
    MyAverageCorrVector = []
    for SampleIndex in range(Nsamples): # Here samples in C object are 0 = Day 35, 1 = Day 70
        
        ######### COMPUTE GENE GENE CORRELATION MATRIX #########
        MyCorrelationMatrix_Full, pvalPearsonCorr_Matrix, MyCorrelationMatrix_Masked = ComputeGeneGeneCorrelationOn2SamplesOfCompareObject(C_35_70,
                                                            C_35_70_norm,
                                                            List_Gene_Names_for_Correlation,
                                                            N_cells,
                                                            ThrasoldForPval,
                                                            RemoveNotSignificantCorr,
                                                            WhiteDiagonal)
        ######### PLOT GENE GENE CORRELATION MATRIX #########
        fig, ax = plt.subplots()
        
        if WhiteDiagonal == False:
            MyCorrelationMatrix = MyCorrelationMatrix_Full
        elif WhiteDiagonal == True:
            MyCorrelationMatrix = MyCorrelationMatrix_Masked
        else:
            print("Error in Masking Correlations!")
        
        MyCorrelationMatrix_Triangle = copy.deepcopy(MyCorrelationMatrix)
        for i in range(len(List_Gene_Names_for_Correlation)):
            for j in range(len(List_Gene_Names_for_Correlation)):
                if j>i:
                    MyCorrelationMatrix_Triangle[i,j] = 0
                if RemoveNANs and numpy.isnan(MyCorrelationMatrix_Triangle[i,j]):
                    print("nan found at element with indexes: " + str(i) + "," + str(j) + " of correlation matrix, SUBSTITUTED with 0")
                    MyCorrelationMatrix_Triangle[i,j] = 0
        
        heatmap = ax.pcolor(MyCorrelationMatrix_Triangle, cmap=plt.cm.seismic_r, vmin=-1, vmax=1) # RdGy_r  # PiYG  # seismic   # PRGn_r
        ax.invert_yaxis()
        for i in range(len(List_Gene_Names_for_Correlation)):
            for j in range(len(List_Gene_Names_for_Correlation)):
                ax.text(- 6, 0.8 + i, List_Gene_Names_for_Correlation[i], fontdict=font3ss)
                ax.text(0.1 + j, - 7, List_Gene_Names_for_Correlation[j], fontdict=font3ss, rotation=90)
                
        #Spacing between each line
        intervals = 1
        loc = plticker.MultipleLocator(base=intervals)
        ax.xaxis.set_major_locator(loc)
        ax.yaxis.set_major_locator(loc)
        ax.grid(which="major", color="black", linestyle='-', linewidth=1)
        ax.tick_params(labelbottom='off')    
        ax.tick_params(labelleft='off')    
        cbar = fig.colorbar(heatmap)
        cbar.set_label('Pearson Correlation', rotation=90)
        """
        i = 0
        j = 0
        for i in [0,1,2,3,4]:
            for j in [0,1,2,3,4]:
                if j>=i:
                    if i != j:
                        if MyCorrelationMatrix[i,j] != 0:
                            Current_pvalPearsonCorr = pvalPearsonCorr_Matrix[i,j]
                            if Current_pvalPearsonCorr <= 0.0001/N_Bonf_Corr:
                                Significance = ' ****'
                            elif Current_pvalPearsonCorr <= 0.001/N_Bonf_Corr:
                                Significance = ' ***'
                            elif Current_pvalPearsonCorr <= 0.01/N_Bonf_Corr:
                                Significance = ' **'
                            elif Current_pvalPearsonCorr <= 0.05/N_Bonf_Corr:
                                Significance = ' *'
                            MyText = str(round(MyCorrelationMatrix[i,j],2)) + Significance
                            ax.text(i+0.5, j+0.5, MyText, 
                                horizontalalignment='center', verticalalignment='center', fontdict=font3)
                        elif MyCorrelationMatrix[i,j] == 0:
                            ax.text(i+0.5, j+0.5, r'$\approx 0$', 
                                horizontalalignment='center', verticalalignment='center', fontdict=font3)
                    elif i == j:
                        if WhiteDiagonal == 0:
                            ax.text(i+0.5, j+0.5, round(MyCorrelationMatrix[i,j],2), 
                                    horizontalalignment='center', verticalalignment='center', fontdict=font3)
                        elif WhiteDiagonal == 1:
                            ax.text(i+0.5, j+0.5,' ', 
                                    horizontalalignment='center', verticalalignment='center', fontdict=font3)
        """
        plt.show()
        """
        Dim = len(MyCorrelationMatrix)
        if WhiteDiagonal == 0:
            AverageCorr = (sum(sum(abs(MyCorrelationMatrix)))-Dim)/(Dim*Dim-Dim)
        elif WhiteDiagonal == 1: # Note that MyCorrelationMatrix contains the diagonal here above, doesn't here below.
            AverageCorr = (sum(sum(abs(MyCorrelationMatrix))))/(Dim*Dim-Dim)
        print(AverageCorr)
        print(' ')
        MyAverageCorrVector.append(AverageCorr)            
        """
    """        
    plt.figure()      
    plt.plot([35,70], MyAverageCorrVector[0:2], marker='o', linestyle='--', color=(255.0/256,0.,0.), label='LRRK2G2019S')
    plt.plot([35,70], MyAverageCorrVector[2:4], marker='s', linestyle='-', color=(128.0/256,128.0/256,128.0/256), label='LRRK2WT')
    plt.xlabel('Time (Days after differentiation)')
    plt.ylabel(r'Average Gene-Gene Correlation')
    plt.ylim(0,0.2)
    plt.title('Average of the absolute value of the Gene-Gene Pearson Correlation Coefficient')
    plt.legend()
    plt.show()
    """
    
    
if DoGeneGeneCorrelation_GenesFromCellTypesProcessesLists == 1:
    print("Step 11b: DoGeneGeneCorrelation_GenesFrom4Lists_cells")
        
    WhiteDiagonal = True # True --> mask values (1s) from the diagonal, False --> do not mask
    RemoveNotSignificantCorr = True
    RemoveNANs = True
    
    MyListOfLists = [ListsOfListsOfGenesCellTypesProcessesActuallyinData[0], 
                     ListsOfListsOfGenesCellTypesProcessesActuallyinData[1], 
                     ListsOfListsOfGenesCellTypesProcessesActuallyinData[2]]
#                     ListsOfListsOfGenesCellTypesProcessesActuallyinData[3], 
#                     ListsOfListsOfGenesCellTypesProcessesActuallyinData[4], 
#                     ListsOfListsOfGenesCellTypesProcessesActuallyinData[5],
#                     ListsOfListsOfGenesCellTypesProcessesActuallyinData[6], 
#                     ListsOfListsOfGenesCellTypesProcessesActuallyinData[7]]
    List_Gene_Names_for_Correlation = [gene for GenesList in MyListOfLists for gene in GenesList]
    
    Nsamples = 2  # Here 2 samples, in fact the samples in C object are 0 = Day 35, 1 = Day 70
    Ngenes = len(List_Gene_Names_for_Correlation)
    N_Bonf_Corr = (Ngenes*Ngenes-Ngenes)/2 * Nsamples 
    ThrasoldForPval = 0.05/N_Bonf_Corr
        
    MyAverageCorrVector = []
    for SampleIndex in range(Nsamples): # Here samples in C object are 0 = Day 35, 1 = Day 70
        
        ######### COMPUTE GENE GENE CORRELATION MATRIX #########
        MyCorrelationMatrix_Full, pvalPearsonCorr_Matrix, MyCorrelationMatrix_Masked = ComputeGeneGeneCorrelationOn2SamplesOfCompareObject(C_35_70,
                                                            C_35_70_norm,
                                                            List_Gene_Names_for_Correlation,
                                                            N_cells,
                                                            ThrasoldForPval,
                                                            RemoveNotSignificantCorr,
                                                            WhiteDiagonal)
        ######### PLOT GENE GENE CORRELATION MATRIX #########
        fig, ax = plt.subplots()
        
        if WhiteDiagonal == False:
            MyCorrelationMatrix = MyCorrelationMatrix_Full
        elif WhiteDiagonal == True:
            MyCorrelationMatrix = MyCorrelationMatrix_Masked
        else:
            print("Error in Masking Correlations!")
        
        MyCorrelationMatrix_Triangle = copy.deepcopy(MyCorrelationMatrix)
        for i in range(len(List_Gene_Names_for_Correlation)):
            for j in range(len(List_Gene_Names_for_Correlation)):
                if j>i:
                    MyCorrelationMatrix_Triangle[i,j] = 0
                if RemoveNANs and numpy.isnan(MyCorrelationMatrix_Triangle[i,j]):
                    print("nan found at element with indexes: " + str(i) + "," + str(j) + " of correlation matrix, SUBSTITUTED with 0")
                    MyCorrelationMatrix_Triangle[i,j] = 0
        
        heatmap = ax.pcolor(MyCorrelationMatrix_Triangle, cmap=plt.cm.seismic_r, vmin=-1, vmax=1) # RdGy_r  # PiYG  # seismic   # PRGn_r
        ax.invert_yaxis()
        for i in range(len(List_Gene_Names_for_Correlation)):
            for j in range(len(List_Gene_Names_for_Correlation)):
                ax.text(- 6, 0.8 + i, List_Gene_Names_for_Correlation[i], fontdict=font3ss)
                ax.text(0.1 + j, - 7, List_Gene_Names_for_Correlation[j], fontdict=font3ss, rotation=90)
                
        #Spacing between each line
        intervals = 1
        loc = plticker.MultipleLocator(base=intervals)
        ax.xaxis.set_major_locator(loc)
        ax.yaxis.set_major_locator(loc)
        ax.grid(which="major", color="black", linestyle='-', linewidth=1)
        ax.tick_params(labelbottom='off')    
        ax.tick_params(labelleft='off')    
        cbar = fig.colorbar(heatmap)
        cbar.set_label('Pearson Correlation', rotation=90)
        

if DoGeneGeneCorrelation_GenesFromCellTypesProcessesLists_OnlyNeuroStem == 1:
    print("Step 11d: DoGeneGeneCorrelation_GenesFromCellTypesProcessesLists_OnlyNeuroStem")
        
    WhiteDiagonal = True # True --> mask values (1s) from the diagonal, False --> do not mask
    RemoveNotSignificantCorr = False
    RemoveNANs = True
    DoubleFace = True
    
    MyListOfLists = [ListsOfListsOfGenesCellTypesProcessesActuallyinData[0], 
                     ListsOfListsOfGenesCellTypesProcessesActuallyinData[1]]
    List_Gene_Names_for_Correlation = [gene for GenesList in MyListOfLists for gene in GenesList]
    
    Nsamples = 2  # Here 2 samples, in fact the samples in C object are 0 = Day 35, 1 = Day 70
    Ngenes = len(List_Gene_Names_for_Correlation)
    N_Bonf_Corr = (Ngenes*Ngenes-Ngenes)/2 * Nsamples 
    ThrasoldForPval = 0.05/N_Bonf_Corr
        
    MyAverageCorrVector = []
    for SampleIndex in range(Nsamples): # Here samples in C object are 0 = Day 35, 1 = Day 70
        
        ######### COMPUTE GENE GENE CORRELATION MATRIX #########
        MyCorrelationMatrix_Full, pvalPearsonCorr_Matrix, MyCorrelationMatrix_Masked = ComputeGeneGeneCorrelationOn2SamplesOfCompareObject(C_35_70,
                                                            C_35_70_norm,
                                                            List_Gene_Names_for_Correlation,
                                                            N_cells,
                                                            ThrasoldForPval,
                                                            RemoveNotSignificantCorr,
                                                            WhiteDiagonal)
        ######### PLOT GENE GENE CORRELATION MATRIX #########
        fig, ax = plt.subplots()
        
        if DoubleFace == False:
            
            if WhiteDiagonal == False:
                MyCorrelationMatrix = MyCorrelationMatrix_Full
            elif WhiteDiagonal == True:
                MyCorrelationMatrix = MyCorrelationMatrix_Masked
            else:
                print("Error in Masking Correlations!")
            
            MyCorrelationMatrix_Triangle = copy.deepcopy(MyCorrelationMatrix)
            for i in range(len(List_Gene_Names_for_Correlation)):
                for j in range(len(List_Gene_Names_for_Correlation)):
                    if j>i:
                        MyCorrelationMatrix_Triangle[i,j] = 0
                    if RemoveNANs and numpy.isnan(MyCorrelationMatrix_Triangle[i,j]):
                        print("nan found at element with indexes: " + str(i) + "," + str(j) + " of correlation matrix, SUBSTITUTED with 0")
                        MyCorrelationMatrix_Triangle[i,j] = 0
                        
        elif DoubleFace == True:
            MyCorrelationMatrix_Triangle = copy.deepcopy(MyCorrelationMatrix_Full)
            for i in range(len(List_Gene_Names_for_Correlation)):
                for j in range(len(List_Gene_Names_for_Correlation)):
                    if j>i:
                        MyCorrelationMatrix_Triangle[i,j] = MyCorrelationMatrix_Masked[i,j]
                    if RemoveNANs and numpy.isnan(MyCorrelationMatrix_Triangle[i,j]):
                        print("nan found at element with indexes: " + str(i) + "," + str(j) + " of correlation matrix, SUBSTITUTED with 0")
                        MyCorrelationMatrix_Triangle[i,j] = 0
                        

        
        heatmap = ax.pcolor(MyCorrelationMatrix_Triangle, cmap=plt.cm.seismic_r, vmin=-1, vmax=1) # RdGy_r  # PiYG  # seismic   # PRGn_r
        ax.invert_yaxis()
        for i in range(len(List_Gene_Names_for_Correlation)):
            for j in range(len(List_Gene_Names_for_Correlation)):
                ax.text(- 6, 0.8 + i, List_Gene_Names_for_Correlation[i], fontdict=font3ss)
                ax.text(0.1 + j, - 7, List_Gene_Names_for_Correlation[j], fontdict=font3ss, rotation=90)
                
        #Spacing between each line
        intervals = 1
        loc = plticker.MultipleLocator(base=intervals)
        ax.xaxis.set_major_locator(loc)
        ax.yaxis.set_major_locator(loc)
        ax.grid(which="major", color="black", linestyle='-', linewidth=1)
        ax.tick_params(labelbottom='off')    
        ax.tick_params(labelleft='off')    
        cbar = fig.colorbar(heatmap)
        cbar.set_label('Pearson Correlation', rotation=90)