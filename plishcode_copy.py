#!/usr/bin/python

#orientation note:
#input mRNA is oriented 3' to 5'; therefore probes produced are oriented 5' to 3'

mrna_sequence_amanda = "CGTTGTTGTTCACATAAACCGTCCTTAGTATGATGACAATATAGCCAGGCGTTGTTGTTCACATAAACCGTCCTTAGTATGATGACAATATAGCCAGG"
mrna_sequence = mrna_sequence_amanda.upper()    #Use a preloaded mRNA sequence

#THIS IS A COMMENT TO CHECK BRANCH SWITCHING

#import libraries: 
#import cmath    #use for mathematics

#import plotting functions
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

import copy
import math
import pandas as pd
import sys
import numpy
import xlsxwriter as xl  #write to excel files
import os
from Bio import Entrez 
from Bio.Blast.Applications import NcbiblastxCommandline    
import PLISHProbeDesigner_copy as probe_designer
import PLISHDesigner_support_copy
from Bio.Blast import NCBIXML
import primer3                    #use this library to calculate thermodynamic values

#used to calculate melting temp
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq


import Tkinter as tk

'''Default Probe Constants: Set by GUI
-------------------------------------------------------------'''

#Provide point of contact email to Entrez
Entrez.email = "tom5@stanford.edu"

species_name = "Gallus gallus"
gene_name = "TECTA"

desktop_path = "/Users/tompritsky/Desktop"

BLAST_xml_path = "/Users/tompritsky/Desktop/HellerLab/PLISH_SCRIPTS/blast_RUN.xml"

#Path to FASTA file for mrna of gene of interest (please set)
FASTA_file = "/Users/tompritsky/Desktop/HellerLab/dataSets/gene_seqs/gg_c14orf180_full_seq.txt"

#Nucleotide database
#nucleotideDatabase="/Users/tompritsky/Desktop/HellerLab/gallus_gallus.fasta"
nucleotideDatabase = "/Users/tompritsky/Desktop/HellerLab/dataSets/Gallus_Gallus_Exons_and_UTR/-Gallus_Gallus_Exons_and_UTR"

#BLAST sequence alignment thresholds (please set):
E_value_threshold = 100

#melt_temp values (please set):
salt_conc = .825
percent_formamide = 30
#melt_temp = 0

#percentGC values (please set):
#min_percent_GC = 55.0149253731343 - 1.5*8.18487976462124     #minimum percentage of GC in a single Hprobe
#max_percent_GC = 55.0149253731343 + 1.5*8.18487976462124    #maximum percentage of GC in a single Hprobe

min_percent_GC = 0     #minimum percentage of GC in a single Hprobe
max_percent_GC = 100    #maximum percentage of GC in a single Hprobe

#end annealment values (please set):
max_end_annealment = 2      #maximum number of overlaps between start and end of h_probe sequence

#Initiator Sequence: Alexa 647
initiator_sequence_dictionary = {'B1':'gAggAgggCAgCAAACgggAAgAgTCTTCCTTTACg', 'B2':'CCTCgTAAATCCTCATCAATCATCCAgTAAACCgCC', 'B3':'gTCCCTgCCTCTATATCTCCACTCAACTTTAACCCg', 'B4':'CCTCAACCTACCTCCAACTCTCACCATATTCgCTTC', 'B5':'CTCACTCCCAATCTCTATCTACCCTACAAATCCAAT'}
initiator_label_list = (initiator_sequence_dictionary.keys())
initiator_label = 'B3'
initiator_sequence = initiator_sequence_dictionary.get(initiator_label, "The initiator label " + str(initiator_label) + " has not been defined, please specify another.")
print("Initiator Sequence: " + initiator_sequence)
initiator_sequence = initiator_sequence.upper() #TROJAN EDIT: Why doesn't .upper on line 587/588 work?

'''Default Probe Constants: Aren't Set by GUI
-------------------------------------------------------------'''

#Spacer Sequence
spacer_sequence = "TT"
spacer_sequence = spacer_sequence.upper()

#Length of probe alignment sequence
alignment_probe_length = 25    #stores length of alignment region on each probe
alignment_sequence_length = alignment_probe_length*2 + 2     #MUST BE EVEN; stores complete length of alignment with mrna sequence, plus two extra nucleotides between probe halves 

#Initiator Sequence List
poss_initiator_sequences = {}

'''Primer3 thermodynamic consts: Will be set by GUI; currently pre-set
---------------------------------------------------------------------------'''
#hairpinDG_min = -0.49763270587821-0.978872719822112;      #minimum value for the hairpin DG formation, in Kcal/mol (more negative DG, more stable hairpin)
hairpinDG_min = -100
homodimerDG_min = -5.19 -2.67;
heterodimerDG_min = -6;
#melt_temp_min = 0
#melt_temp_max = 100

#melt_temp_min = 64.3692549489138 - 0*3.35580070349471      #need to fill in 
#melt_temp_max = 64.3692549489138 + 2*3.35580070349471

melt_temp_min = 0      #need to fill in 
melt_temp_max = 100


'''Internal variables:
--------------------------------------------------------------'''
path = desktop_path + "/PLISH"     #stores path to directory where plish output will be stored.
                                   #The directory is a desktop file called PLISH. 

zero_alignment_full_binding_sequences = {}  #store full length binding sequences with no BLAST alignments (reverse complement)
one_alignment_full_binding_sequences = {}   #store binding sequences with one BLAST alignment (reverse complement)
two_alignment_full_binding_sequences = {}   #store binding sequences with two BLAST alignments (reverse complement)

mult_alignment_full_binding_sequences = {}  #store sequences with multiple BLAST alignments 
specific_40bp_sequences_reverseComp = {}        #store all sequences with 0, 1, or 2 alignments (reverse complement)
potential_Hprobe_pairs = []     #stores potential Hprobe pairs
valid_Hprobe_pairs = []     #stores valid Hprobe pairs

'''probe offset preparation (NEW):
_________________________________________________________________'''

potential_overlap_dict = {}       #store potential overlap dictionary for each probe
valid_key_list = []     #stores indices of all successful probes, after filtering and BLAST alignment (reported from 3' end of gene of interest)
actual_overlap_dict = {}    #stores dictionary with key = probe indice and value = num_overlaps
adjusted_overlap_dict = {} #stores dictionary with key = probe indice and value = list of valid overlapping probes
valid_key_dictionary = {}   #stores each valid probe indice, as well as the number of BLAST alignments for that probe

#create dictionary of every potential valid indice and all it's potential overlaps. Receives arguments of 
#length (length of alignment mrna sequence), size (length of left and right probes and two nucleotide spacer),
#and overlap_max_nucleotides (the maximum number of valid nucleotide overlaps between two probe sets)
def create_potential_overlap_dict(length, size, overlap_max_nucleotides):
  overlap_dict = {}
  valid_indices = []
  for indice in range(0, length):
        if int(indice) + size <= length:    #if int(indice) - 25 >= 0 and len(mrna_sequence) - int(indice) >=25: 
            valid_indices.append(indice)
  for i in valid_indices:
    overlap_list = list(range(i + 1, i+size - overlap_max_nucleotides))
    overlap_list_to_keep = []
    for j in overlap_list:
      if(j in valid_indices):
        overlap_list_to_keep.append(j)
    overlap_dict.update({i:overlap_list_to_keep})
  return overlap_dict

#create dictionary of each valid indice and a list of all the valid indices it overlaps
def create_adjusted_overlap_dict(pot_overlap_dict, val_key_list):
  adjusted_overlap_dict = copy.deepcopy(pot_overlap_dict)
  for key, value in pot_overlap_dict.items():
    for indice in value:
      if(indice not in val_key_list):
        adjusted_overlap_dict[key].remove(indice)
  return adjusted_overlap_dict

#create dictionary storing every valid indice and the number of overlaps it has with other valid indices
def create_actual_overlap_dict(pot_overlap_dict, val_key_list):
  overlap_dict = {}
  for key in val_key_list:
    num_overlaps = 0
    for value in pot_overlap_dict[key]:
      if value in val_key_list:
        num_overlaps += 1
    overlap_dict.update({key: num_overlaps}) 
  return overlap_dict

def filter_by_overlaps(valid_key_list):
    potential_overlap_dict = create_potential_overlap_dict(len(mrna_sequence), alignment_sequence_length, overlap_max_nucleotides = -2)
    actual_overlap_dict = create_actual_overlap_dict(potential_overlap_dict, valid_key_list)    #stores list of valid probe dictionaries with key = probe indice and value = num_overlaps
    adjusted_overlap_dict = create_adjusted_overlap_dict(potential_overlap_dict, valid_key_list)   #list of probe dictionaries with key = probe indice and value = valid probe overlap indices (invalid probe overlaps are removed)
    
    temp_valid_key_list = copy.deepcopy(valid_key_list)
    temp_valid_key_list.sort()
    current_indice = 0
    while(current_indice < len(temp_valid_key_list)):  
      #print("Iteration " + str(current_indice))
      current_key = temp_valid_key_list[current_indice]
      #print("current key " + str(current_key))
      #print("current indice: " +str(current_indice))
      if( not any(i in adjusted_overlap_dict[current_key]  for i in temp_valid_key_list)):
        current_indice += 1
        #print("Current indice has no overlaps")
        #print("\n")
      else:
        print("Temp valid key list: " + str(temp_valid_key_list))
        formated_section_temp_valid_key_list = temp_valid_key_list[0:current_indice]    #store segment of list that is already formatted
        unformatted_section_temp_valid_key_list = temp_valid_key_list[current_indice:]
        #print("unformatted_section_temp_valid_key_list: " +str(unformatted_section_temp_valid_key_list))
        unformatted_section_temp_valid_key_list = filter_by_overlaps_helper(unformatted_section_temp_valid_key_list, adjusted_overlap_dict, actual_overlap_dict, 0)
        #print("unformatted_section_temp_valid_key_list: " +str(unformatted_section_temp_valid_key_list))
        temp_valid_key_list = formated_section_temp_valid_key_list + unformatted_section_temp_valid_key_list
        
        #print("updated Temp_valid key list: " + str(temp_valid_key_list))
        if( set(adjusted_overlap_dict[temp_valid_key_list[0]]).isdisjoint(temp_valid_key_list)):         #
          #print("moved to next indice")  
          current_indice += 1 
    return temp_valid_key_list

#this function removes all overlapping probes from a set of overlapping probes, excluding the probe with least overlaps total
def filter_by_overlaps_helper(valid_key_list, adjusted_overlap_dict, actual_overlap_dict, current_indice):
  
  temp_list_new = copy.deepcopy(valid_key_list)
  
  current_key = valid_key_list[current_indice]
  #print("FBO: current key: " + str(current_key))
  min_overlaps = actual_overlap_dict[current_key]  #get number of overlaps of first valid probe
  min_overlap_indice = current_key      #get indice of first valid probe
  #print(min_overlaps)
  
  #iterate through all valid overlaps for first valid probe
  for valid_key_overlap in adjusted_overlap_dict[current_key]:
    if(actual_overlap_dict[valid_key_overlap] < min_overlaps):    #check if a given overlapping probe has less overlaps than the current probe with least overlaps
      min_overlaps = actual_overlap_dict[valid_key_overlap]       #adjust min overlaps
      min_overlap_indice = valid_key_overlap         #adjust indice of minimal overlap probe
  

  #print("temp_list_new " + str(temp_list_new))
  #print("adjusted_overlap_dict: " + str(adjusted_overlap_dict))
  #print("adjutsed_overlap_dict of current key: " + str(adjusted_overlap_dict[current_key]))
  #print("actual overlap dict: " + str(actual_overlap_dict))
 
  #print("Min overlap indice: " + str(min_overlap_indice))
  
  #remove all overlaps of min_overlap_indice 
  temp_list_new = list(set(temp_list_new)- set(adjusted_overlap_dict[min_overlap_indice]))
  temp_list_new.sort()
  #print("Temp List New: " + str(temp_list_new))
  
  #remove all indices before min_overlap_indice, unless min_overlap indice is first indice 
  if(min_overlap_indice != current_key):
    new_lowest_indice = temp_list_new.index(min_overlap_indice)
    temp_list_new = temp_list_new[new_lowest_indice:]
  
  #temp_list_new = list(set(temp_list_new) & set(set(temp_list_new)^set(adjusted_overlap_dict[min_overlap_indice])))
  #print("temp_list new: " + str(temp_list_new))
  
  #print(temp_list_new)
  temp_list_new.sort()
  return temp_list_new

'''Internal Classes: Contains the hprobe and Hprobe_pair classes, which are 
used to store data about the Hprobes designed by the PLISH algorithm
---------------------------------------------------------------------------'''
            
'''Hprobe: The class that defines an hprobe, a 20 base pair sequence that anneals to RNA regions of interest. 
           Contains the member variables index, sequence, percent_GC, melt_temp, and end_annealment. The hprobes
           are oriented from 5' to 3' and include an initiator sequence, a spacer, and a binding sequence.'''
class Hprobe:
    is_left_probe = None     #stores True for probe aligning closest to 3' end of mrna, False for probe aligning closer to 5' end
    complete_sequence = ""      #complete probe sequence
    start_index = None    #initial base pair indice of RNA sequence where hprobe aligns
    binding_sequence = ""         #stores target binding sequence for the hprobe
    percent_CG = None     #the percentage of 'G' or 'C' nucleotides in hprobe
    melt_temp = None      #the melting temperature for a probe
    end_annealment_count = None   #the number of sequence overlaps between the start and end of the hprobe
    hairpinDG = None
    homodimerDG = None
    
    #constructor for hprobe class
    def __init__(self, is_left_probe, complete_sequence, index, binding_sequence, percent_GC, melt_temp, end_annealment, hairpinDG, homodimerDG):
        self.is_left_probe = is_left_probe
        self.complete_sequence = complete_sequence
        self.start_index = index
        self.binding_sequence = binding_sequence
        self.percent_GC = percent_GC
        self.melt_temp = melt_temp
        self.end_annealment_count = end_annealment
        self.hairpinDG = hairpinDG
        self.homodimerDG = homodimerDG

#make_Hprobe: Used to create new instances of the hprobe class
def make_Hprobe(is_left_probe, complete_sequence, index, binding_sequence, percent_GC, melt_temp, end_annealment, hairpinDG, homodimerDG):
    hprobe = Hprobe(is_left_probe, complete_sequence, index, binding_sequence, percent_GC, melt_temp, end_annealment, hairpinDG, homodimerDG)
    return hprobe         
    
'''Hprobe_pair: A container class for Hprobe pairs. Stores the left and right hprobes of a pair, as 
                well as a num_alignments member variable that holds the number of alignments between 
                the hprobe pair and the mrna.'''     
class Hprobe_pair:
    left = Hprobe
    right = Hprobe
    num_alignments = None
    sequence_number = None #the index of the sequence, which indicates relative order in which sequences were processed
    heterodimerDG = None
    
    def __init__(self, left_probe, right_probe, num_alignments, sequence_number, heterodimerDG):
        self.left = left_probe
        self.right = right_probe
        self.num_alignments = num_alignments
        self.sequence_number = sequence_number
        self.heterodimerDG = heterodimerDG

#make_Hprobe_pair: Used to create new instances of the Hprobe_pair class
def make_Hprobe_pair(left_probe, right_probe, num_alignments, sequence_number, heterodimerDG):
    hprobe_pair = Hprobe_pair(left_probe, right_probe, num_alignments,sequence_number, heterodimerDG)
    return hprobe_pair


'''Calculate probe metrics: Calculates probe parameters and outputs as a graph. '''

def probeMetrics(hprobe_pairs): 
    counter = 0;
    listOfRowLists = []
    columnNames = ['Hprobe', 'Sequence Number', 'Number of Alignments', 'mRNA Indice', 'Sequence (initiator and spacer(tt) lowercase)', 'Percent GC', 'Melting Temperature', 'HairpinDG (Kcal/mol)']
    
    for pair in hprobe_pairs:
        rowListLeft = []
        rowListRight = []
        left = pair.left     #store left probe in left
        right = pair.right   #store right probe in right

        #set left probe's title
        left_probe_title = ("Probe " + str(counter) + " (left)")
        #set right probe's title
        right_probe_title = ("Probe " + str(counter) + " (right)")
        #set overall probe title
        
        #append data to left probe output array
        rowListLeft.append(left_probe_title)  #append left probe title
        rowListLeft.append(pair.sequence_number)  #append the sequence number
        rowListLeft.append(pair.num_alignments) #the number of mrna alignments (total number in probe)
        #rowListLeft.append(left.end_annealment_count) #the number of probe end overlaps
        rowListLeft.append(left.start_index)  #initial base pair indice of RNA sequence where hprobe aligns
        rowListLeft.append(left.complete_sequence)     #the base pair sequence for the hprobe
        rowListLeft.append(left.percent_GC)   #the percentage of 'G' or 'C' nucleotides
        rowListLeft.append(left.melt_temp)    #the melting temperature for a probe
        rowListLeft.append(left.hairpinDG)   #the hairpin delta G value for the hprobe
        #rowListLeft.append(left.homodimerDG)    #the homodimer delta G value
        #rowListLeft.append(pair.heterodimerDG)    #the heterodimer delta G value
        
        #append data to right probe output array
        rowListRight.append(right_probe_title)  #append right probe title
        rowListRight.append(pair.sequence_number)  #append the sequence number
        rowListRight.append(pair.num_alignments) #the number of mrna alignments (total number in probe)
        #rowListRight.append(right.end_annealment_count) #the number of probe end overlaps
        rowListRight.append(right.start_index)  #initial base pair indice of RNA sequence where hprobe aligns
        rowListRight.append(right.complete_sequence)     #the base pair sequence for the hprobe
        rowListRight.append(right.percent_GC)   #the percentage of 'G' or 'C' nucleotides
        rowListRight.append(right.melt_temp)    #the melting temperature for a probe
        rowListRight.append(right.hairpinDG)   #the hairpin delta G value for the hprobe
        #rowListRight.append(right.homodimerDG)    #the homodimer delta G value
        
        listOfRowLists.append(rowListLeft)
        print("probeMetrics running")
        print(rowListLeft)
        listOfRowLists.append(rowListRight)
        counter = counter + 1
    
    probe_results = pd.DataFrame(listOfRowLists, columns=columnNames)
    print(probe_results)
    
    PercentGCList = probe_results['Percent GC'].tolist()
    MeltingTempList = probe_results['Melting Temperature'].tolist()
    HairpinDGList = probe_results['HairpinDG (Kcal/mol)'].tolist()
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(PercentGCList, MeltingTempList, HairpinDGList)
    
    ax.set_xlabel('Percent_GC')
    ax.set_ylabel('Melting Temperature')
    ax.set_zlabel('HairpinDG (Kcal/mol)')
    
    #PercentGCList, MeltingTempList = np.meshgrid(PercentGCList, MeltingTempList)
    #ax.plot_surface(PercentGCList, MeltingTempList, HairpinDGList, rstride=1, cstride=1, cmap=cm.viridis)
    plt.title("Probe Distribution")
    plt.show()
    
    #plt.hist(PercentGCList)
    #plt.show()
    #plt.hist(MeltingTempList)
    #plt.title("c14orf Melting Temp")
    #plt.show()
    #plt.hist(HairpinDGList)
    #plt.title("c14orf Hairpin_DG")
    #plt.show()
    

'''Internal Functions: Used to perform probe operations and calculate probe parameters, such as percent
GC and end annealment. 
---------------------------------------------------------------------------------------------------'''

#reverseString: returns the reverse ordered string. Used to maintain 5' to 3' orientation for Hprobes
def reverseString(new_string):
    return new_string[::-1]

#return hairpin dG from base pair sequence. Can also return melting temperature if needed
def calcHairpinDG(bp_sequence):
    hairpin_parameters = primer3.calcHairpin(bp_sequence)
    return (hairpin_parameters.dg)/1000

#return homodimer dG from base pair sequence. Can also return melting temperature if needed
def calcHomodimerDG(bp_sequence):
    homodimer_parameters = primer3.calcHomodimer(bp_sequence)
    return (homodimer_parameters.dg)/1000

#return heterodimer dG from single base pair sequences. Can also return melting temperature if needed
def calcHeterodimerDG(sequence1, sequence2):
    heterodimer_parameters = primer3.calcHeterodimer(sequence1, sequence2)
    return (heterodimer_parameters.dg)/1000

#reverseComplement: returns the reverse complement of a string. Used to find complementary Hprobe sequences
def reverseComplement(oldString):
    string_list = oldString
    newstring = ""
    for index in string_list:
        if index == 'A':
            newstring+="T"
        elif index == 'T':
            newstring+="A"
        elif index == 'C':
            newstring+="G"
        elif index == 'G':
            newstring+="C"
    return newstring

'''percentGC: Calculates the percentage of G or C nucleotides in the hprobe'''
def percentGC(bp_sequence):
    GC_counter = 0
    percent_GC = 0
    for index in range(len(bp_sequence)):
        if bp_sequence[index].upper() == 'G' or bp_sequence[index].upper() == 'C':
            GC_counter = GC_counter + 1
    percent_GC = float(GC_counter)/len(bp_sequence)*100
    return percent_GC

'''valid_percentGC: Ensures that the proportion of G/C pairs is within acceptable limits. Limits are set
                    via global variables min_percent_GC and max_percent_GC'''
def valid_percentGC(Hprobe_pair):
    if (min_percent_GC <= Hprobe_pair.left.percent_GC and Hprobe_pair.left.percent_GC <= max_percent_GC and min_percent_GC <= Hprobe_pair.right.percent_GC and Hprobe_pair.right.percent_GC <= max_percent_GC):
        return True
    else:
        return False

'''end_annealment_count: Counts total number of start/end sequence overlaps for an hprobe'''
def end_annealment_count(bp_sequence):
    sequence_overlaps = 0
    for i in range(0,2):
        if bp_sequence[i].upper() == bp_sequence[len(bp_sequence)-1-i].upper():
            sequence_overlaps = sequence_overlaps + 1
        else:
           break
    return sequence_overlaps

'''valid_end_annealment: Verifies that end annealment is within acceptable bounds. Maximum end
                        annealment is set by global variable max_end_annealment'''
def valid_end_annealment(Hprobe_pair):
    if Hprobe_pair.left.end_annealment_count > max_end_annealment or Hprobe_pair.right.end_annealment_count > max_end_annealment:
        return False
    else:
        return True

#returns the melting temperature and delta G value for hairpin formation
def get_Hairpins(bp_sequence):
    hairpin_parameters = primer3.calcHairpin(bp_sequence)
    return hairpin_parameters.tm, hairpin_parameters.dg

'''get_melt_temp: returns the integer melting temperature for a given Hprobe. Global variables set for salt_conc,
                  percent_GC, bp_sequence, percent_formamide'''
def get_melt_temp(salt_conc, percent_GC, bp_sequence, percent_formamide):
   Tm = 81.5 + 16.6*(math.log10(salt_conc)) + .41*(percent_GC) - 500/len(bp_sequence) - .61*percent_formamide    #source: http://bioinfo.ut.ee/primer3-0.4.0/input-help.htm
   return Tm.real
   '''Old Formula'''
   #n = len(bp_sequence) #get length of BP sequence
   #salt_conc = numpy.float(salt_conc)
   #percent_formamide = numpy.float(percent_formamide)
   #percent_GC = numpy.float(percent_GC)
   #melt_temp = 81 + 16.6*(cmath.log(salt_conc)) + .41*(percent_GC) - 500/n - .61*(percent_formamide) #use Wahl paper formula

'''valid_melt_temp: Ensure that melting temperature of H_probes is within an allowable range. 
                    Please make sure to set the min_melt_temp and max_melt_temp global variables. '''
def valid_melt_temp(Hprobe_pair):
    if Hprobe_pair.left.melt_temp > 2 or Hprobe_pair.melt_temp > 2:
        return False
    else:
        return True

'''Produce Summary File: Produce excel output sheet storing probe alignment results.
-----------------------------------------------------------------------------------------------'''

'''initializeWorksheet: Produces worksheet to store PLISH output. Initializes the column 
    headers with labels. Stores the output file in the PLISH folder on the desktop,
    with the path to the desktop directory specified by the global variable desktop_path.
    Returns the workbook and the initialized worksheet.'''
def initializeWorkSheet():
    #create list of column labels for output file
    label_list1 = ('Hprobe', 'Sequence Number', 'Number of Alignments', 'mRNA Indice', 'Sequence (initiator and spacer(tt) lowercase)', 'Percent GC', 'Melting Temperature', 'HairpinDG (Kcal/mol)')
    
    #Create excel workbook and add worksheet
    workbook = xl.Workbook(path + '/PLISH_Workbook_2.xlsx')
    
    #create worksheet to store left/right h_probes
    worksheet1 = workbook.add_worksheet()
    
    '''IMPLEMENT'''
    #create worksheet to store full h_probe pairs 
    worksheet2 = workbook.add_worksheet()
    
    #set iterators for row, column
    row = 0
    column = 0

    #write labels to sheet column by column
    for label in label_list1:
        worksheet1.write(row, column, label)
        column += 1
    
    return workbook, worksheet1, worksheet2

#assign ranking to each probe
    
'''plotProbeData: Adds complete probe information for all probes with 0, 1, or 2
   overlaps to the output excel sheet. Takes in hprobe_pairs: the dataset of h_probes;
   worksheet1: stores information about left and right h_probes individually '''
def plotProbeData(hprobe_pairs, worksheet1, worksheet2, workbook):
    output_array_left = []   #store data for left probe
    output_array_right = []  #store data for right probe
    counter = 0              #count number of probe pairs
    row = 1     #set row counter
    column = 0  #set column counter
    
    #The cell_format parameter in the sub write methods is used to apply formatting to the cell. This parameter is optional but when present it should be a valid Format object)
    #cell_format = workbook.add_format({'bold': True, 'italic': True})
    
    # Create format: this corresponds to out of range probe variables
    format0 = workbook.add_format({'font_color': 'black'})
    format1 = workbook.add_format({'bg_color': '#FFC7CE', 'font_color': '#9C0006'})
    format2 = workbook.add_format({'bg_color': '#00FF00', 'font_color': 'black'})
    
    for pair in hprobe_pairs:
        left = pair.left     #store left probe in left
        right = pair.right   #store right probe in right

        #set left probe's title
        left_probe_title = ("Probe " + str(counter) + " (left)")
        #set right probe's title
        right_probe_title = ("Probe " + str(counter) + " (right)")
        #set overall probe title
        
        #append data to left probe output array
        output_array_left.append(left_probe_title)  #append left probe title
        output_array_left.append(pair.sequence_number)  #append the sequence number
        output_array_left.append(pair.num_alignments) #the number of mrna alignments (total number in probe)
        #output_array_left.append(left.end_annealment_count) #the number of probe end overlaps
        output_array_left.append(left.start_index)  #initial base pair indice of RNA sequence where hprobe aligns
        output_array_left.append(left.complete_sequence)     #the base pair sequence for the hprobe
        output_array_left.append(left.percent_GC)   #the percentage of 'G' or 'C' nucleotides
        output_array_left.append(left.melt_temp)    #the melting temperature for a probe
        output_array_left.append(left.hairpinDG)   #the hairpin delta G value for the hprobe
        #output_array_left.append(left.homodimerDG)    #the homodimer delta G value
        #output_array_left.append(pair.heterodimerDG)    #the heterodimer delta G value
        
        
        #append data to right probe output array
        output_array_right.append(right_probe_title)  #append right probe title
        output_array_right.append(pair.sequence_number)  #append the sequence number
        output_array_right.append(pair.num_alignments) #the number of mrna alignments (total number in probe)
        #output_array_right.append(right.end_annealment_count) #the number of probe end overlaps
        output_array_right.append(right.start_index)  #initial base pair indice of RNA sequence where hprobe aligns
        output_array_right.append(right.complete_sequence)     #the base pair sequence for the hprobe
        output_array_right.append(right.percent_GC)   #the percentage of 'G' or 'C' nucleotides
        output_array_right.append(right.melt_temp)    #the melting temperature for a probe
        output_array_right.append(right.hairpinDG)   #the hairpin delta G value for the hprobe
        #output_array_right.append(right.homodimerDG)    #the homodimer delta G value
        
        
        '''FIX MELT TEMP CHECK'''
        #Color first column in row green if hprobe is valid. NOT NEEDED SINCE PROBES PRETESTED
        '''
        worksheet1.write(row, 0, output_array_left[0], format2 
                        if(output_array_left[2]<=2
                           and output_array_left[3]<=max_end_annealment
                           and output_array_left[4]>=min_percent_GC
                           and output_array_left[4]<=max_percent_GC
                           and output_array_left[5]>=0)
                        else None)
        '''
        #column+=1
        #write data output to excel file. 
        for output in output_array_left[0:]:
            worksheet1.write(row, column, output)
            column += 1
        
        #reset the rows and columns
        column = 0
        row += 1
        
        #write data output to excel file. Color first column in row green if probe is valid: NOT NEEDED SINCE PROBES PRETESTED
        '''
        worksheet1.write(row, 0, output_array_right[0], format2 
                        if(output_array_right[2]<=2
                           and output_array_right[3]<=max_end_annealment
                           and output_array_right[4]>=min_percent_GC
                           and output_array_right[4]<=max_percent_GC
                           and output_array_right[5]>=0)
                        else None)
        '''
        #column+=1
        #write data to excel file
        for output in output_array_right[0:]:
            worksheet1.write(row, column, output)
            column += 1
        
        #iterate the rows, columns and counter
        column = 0
        row += 2
        counter += 1
        
        #clear left and right arrays
        output_array_left = []
        output_array_right = []
    
    
    
    workbook.close()
    
    # conditionally format worksheet1 to indicate erroneous hprobe values
    '''
    #NUM_ALIGNMENTS: format if too many alignments: NOT NEEDED SINCE REQUIREMENTS ARE TESTED BEFORE OUTPUT
    worksheet1.conditional_format("$C$2:$F$%d"% ((len(hprobe_pairs)*3)+1), {'type':'blanks', 'format': format0, 'stop_if_true': True})
    worksheet1.conditional_format("$C$2:$C$%d"% ((len(hprobe_pairs)*3)+1), {'type': 'cell',
        'criteria': '>',
        'value': 2,
        'format': format1})
    # MAX_END_ANNEAL: format if end annealment is greater than max_end_annealment
    worksheet1.conditional_format("$D$2:$D$%d"% ((len(hprobe_pairs)*3)+1), {'type': 'cell',
        'criteria': '>',
        'value': max_end_annealment,
        'format': format1})
    # MIN_GC: color if percent GC is less than max_percent_GC
    worksheet1.conditional_format("$E$2:$E$%d"% ((len(hprobe_pairs)*3)+1), {'type': 'cell',
        'criteria': '<',
        'value': min_percent_GC,
        'format': format1})
    # MAX_GC: color if percent GC is greater than max_percent_GC
    worksheet1.conditional_format("$E$2:$E$%d"% ((len(hprobe_pairs)*3)+1), {'type': 'cell',
        'criteria': '>',
        'value': max_percent_GC,
        'format': format1})
    # FIX: Conditionally format temperature
    worksheet1.conditional_format("$F$2:$F$%d"% ((len(hprobe_pairs)*3)+1), {'type': 'cell',
        'criteria': '>',
        'value': max_percent_GC,
        'format': format1})
    '''

'''CHECK IF NEEDED'''
#def addRowFormatting():
    
from Bio import SeqIO 
#handle = Entrez.efetch(db="gene", id="170609", rettype="db", retmode="txt")
#using this to download sequence from online

'''WRITE VALIDATE GUI INPUTS FUNCTION'''


'''Upload FASTA file: Upload a FASTA file for mrna of gene of interest to script,
   and print the description to the screen. The description includes the name and
   the sequence ID. Returns string containing the mrna sequence for gene of interest
   and a record ID from the original FASTA file.'''

'''EDIT: What is the order of input mrna. Must be 3' to 5' '''
   
def read_mrna_sequence():
    handle = FASTA_file   #set the FASTA file directory by adjusting global variable FASTA_file
    record = SeqIO.read(handle, "fasta")    #create record from FASTA file to allow better parsing
    
    #EDITED: store record ID of mRNA sequence
    print("my_Record")
    record_ID = record.id
    
    #handle.close()
    
    #print description of mRNA sequence
    print("mRNA Sequence description:")
    print(record.description)
    
    #convert the mrna sequence to string
    mrna_sequence = str(record.seq)
    #global mrna_sequence
    return mrna_sequence, record_ID

'''extractSequences: Returns a dictionary with 50BP sequences as values and their
indice locations as keys'''
def extractSequences(mrna_sequence):
    
    #create array to store valid start indices (so hprobes remain within bounds of mrna )
    valid_indices = []
    
    #test and print the start indices to ensure hprobes remain within bounds of mrna 
    '''EDIT'''
    for indice in range(0, len(mrna_sequence)-1):
        if int(indice) + alignment_sequence_length + 2 <= len(mrna_sequence):    #if int(indice) - 25 >= 0 and len(mrna_sequence) - int(indice) >=25: add 2 nucleotides to account for spacer between left and rigt probe halves
            valid_indices.append(indice)
    print("\nValid Start Indices:")   #prints valid start indices for potential hprobe pair
    print(valid_indices)     
    
    #create dictionary to store full length binding sequences for potential hprobe pair
    potential_full_binding_sequences = {}
    
    '''create potential 40 base pair hprobe sequences and add to dictionary:
        key equals the indice and value equals the sequence'''
    '''EDIT'''
    for i in valid_indices:
        potential_full_binding_sequences[i] = mrna_sequence[i:i+alignment_sequence_length + 2]  
    
    #Print all potential hprobe sequences
    print("\nPotential Hprobe Sequences:")
    for key, value in potential_full_binding_sequences.iteritems():
        print value   #print all potential h_probe sequences
    return potential_full_binding_sequences

'''Create all H_probe objects: Store potential h_probes as h_probe objects
   with member variables defined by the H_probe pair class. Two h_probes, 
   left and right h_probe, are stored in the container class hprobe_pair,
   which stores the left and right halves of an h_prober object.'''
def test_H_Probe(key, value, num_iterations):
    print("create all h probes running")
    print(str(100*key/num_iterations) + "% complete")
    
    #create reverse complement for testing
    value = reverseComplement(value)
    
    valid_probe = True
    #CHECK: if reverseComplement hprobe is necessary (currently commented out)
    #create h_probe objects corresponding to the left and right halves of an h_probe
    #edit
    left_hprobe_binding_sequence = value[0:alignment_probe_length]  #get corresponding mrna sequence (left probe)
    #left_hprobe_sequence = reverseComplement(left_hprobe_sequence)  #get complementary hprobe sequence (left probe)
    left_hprobe_binding_sequence = reverseString(left_hprobe_binding_sequence)  #reverse directionality (left probe)
    right_hprobe_binding_sequence = value[alignment_probe_length+2:alignment_sequence_length]
    #left_hprobe_sequence = reverseComplement(left_hprobe_sequence)
    right_hprobe_binding_sequence = reverseString(right_hprobe_binding_sequence)
    
    #calculate percent GC values for each hprobe binding sequence
    percent_GC_left = percentGC(left_hprobe_binding_sequence)
    percent_GC_right = percentGC(right_hprobe_binding_sequence)
    
    #test percent GC values for each hprobe binding sequence
    if(percent_GC_left < min_percent_GC or percent_GC_left > max_percent_GC or 
       percent_GC_right < min_percent_GC or percent_GC_right > max_percent_GC):
        valid_probe = False
    
    '''Generate full hprobe sequences ordered 5' to 3' (consisting of binding sequence, initiator sequence, and spacer). 
    These full hprobe are used in validity testing, including percent GC, end_annealment, etc.'''
    left_hprobe_sequence = initiator_sequence[0:len(initiator_sequence)/2] + spacer_sequence + left_hprobe_binding_sequence
    right_hprobe_sequence = right_hprobe_binding_sequence + spacer_sequence + initiator_sequence[len(initiator_sequence)/2:len(initiator_sequence)]
    
    #left_hprobe_sequence = left_hprobe_sequence.upper()     #uppercase entire sequence
    #right_hprobe_sequence = right_hprobe_sequence.upper()   #uppercase entire sequence
    
    #calculate hairpin delta G values for each hprobe
    hairpinDG_left = calcHairpinDG(left_hprobe_sequence)
    hairpinDG_right = calcHairpinDG(right_hprobe_sequence)
    
    #test hairpinDG values:
    if(hairpinDG_left < hairpinDG_min or hairpinDG_right < hairpinDG_min):
        valid_probe = False
    
    #calculate homodimer delta G values for each hprobe
    #homodimerDG_left = calcHomodimerDG(left_hprobe_sequence)
    #homodimerDG_right = calcHomodimerDG(right_hprobe_sequence)
    
    #test homodimerDG values:
    #if(homodimerDG_left < homodimerDG_min or homodimerDG_right < homodimerDG_min):
        #valid_probe = False
    
    #calculate heterodimer delta G values
    #heterodimerDG = calcHeterodimerDG(left_hprobe_sequence, right_hprobe_sequence)
    
    #test heterodimer delta G values
    #if(heterodimerDG < heterodimerDG_min):
        #valid_probe = False
    
    '''record and store the melt_temp for each h_probe. Ensure that the melt_temp equation
    is correct. Also verify that the global variables used to calculate the melt_temp are set 
    correctly.'''
    melt_temp_left = get_melt_temp(salt_conc, percent_GC_left, left_hprobe_sequence, percent_formamide)
    melt_temp_right = get_melt_temp(salt_conc, percent_GC_right, right_hprobe_sequence, percent_formamide)
    
    #need to implement this
    if(melt_temp_left < melt_temp_min or melt_temp_right < melt_temp_min or melt_temp_left > melt_temp_max or melt_temp_right > melt_temp_max):
        valid_probe = False
        
    return valid_probe      #return whether probe is valid (passes thermodynamic tests)
        
'''runBLAST: runs BLAST alignment on a nucleotide sequence and store the string and its location into the correct vector
based upon number of alignments. This function prints runtime information, but has no return values. The input is a string 'value',
storing the nucleotide sequence, an integer key storing its location, and an integer counter.'''
def runBLAST(key, value, counter, record_ID, FASTA_file_path = desktop_path+"/HellerLab/PLISH_SCRIPTS/tmp.fasta"):
    print("\nPotential Hprobe sequence " + str(counter+1) + ":")     #prints the sequence being aligned
    print(value)    #NOT Sure why this doesn't return the first value
    
    #Produce a temporary fasta file because the blast command needs a fasta file as input.
    out_fasta = open(FASTA_file_path, 'w')
    query = ">" + record_ID + "\n" + str(value)
    print("\nWriting to temp FASTA file...")
    out_fasta.write(query) # And I have this file with my sequence correctly filled
    out_fasta.close() # EDIT : It was out of my loop....
    
    #create command line argumnt
    print("\nRunning BLAST...")
    blastn_cline = NcbiblastxCommandline(cmd='/usr/local/ncbi/blast/bin/blastn', query=desktop_path+"/HellerLab/PLISH_SCRIPTS/tmp.fasta", db= nucleotideDatabase, evalue=0.001, outfmt=5, out=BLAST_xml_path)
    print(blastn_cline)
    #output_BLAST.close()
    
    #run blast
    print("BLAST Output")
    #stdout, stderr = blastn_cline()
    os.system(str(blastn_cline))
    sortBLASTOutput(key, value, counter, record_ID)
    
def sortBLASTOutput(key, value, counter, record_ID):
    result_handle = open(BLAST_xml_path, 'r')
    
    #Alternate result handle to check for two alignment sequences
    #result_handle = open("/Users/tompritsky/Desktop/PLISH_output/blast_RUN.xml", 'r')
    
    #loop through BLAST alignments and print them to screen
    blast_record = NCBIXML.read(result_handle) #how do I exclude specific sequences
    
    potential_alignments = [] #stores potential sequence alignments
    
    '''EDIT Validation to be more hollistic'''
    '''validate potential alignment: Ensure that the alignment is from the correct species but not
                                      from a homologous gene. Does this by ensuring that the species 
                                      name is in the alignment title and the gene name is not.'''                                    
    for alignment in blast_record.alignments:
        #if gene_name not in alignment.title and species_name in alignment.title:
        potential_alignments.append(alignment)
        
        print alignment
    
    #create tuple storing sequence [first element] and sequence number [second element]. The sequence number
    #is just a number to determine the order in which probes are passed to the BLAST alignment function after
    #the filtering step, and do not reflect the alignment location of the probe to the mRNA sequence of the gene 
    #of interest. Rather, the mRNA alignment location is given by the 'key'. 
    value = [value, counter]
    
    '''store sequences with 0,1,2, and multiple alignment counts. Sequences are stored
    in a dictionary ____alignments depending on the number of alignments. The key is 
    equal to the sequence start indice and the values equivalent to the nucleotide sequence.
    Additionally, store all valid probes into valid_key_dictionary dictionary so that overlaps can
    be reomved by filter_by_overlaps function. Keys in valid key dictionary are the sequence start
    indice and values are the number or alignments.'''
    if len(potential_alignments) == 0:
        zero_alignment_full_binding_sequences[key] = value
        valid_key_dictionary.update({key:len(potential_alignments)})
        print ("found zero alignment sequence")
    elif len(potential_alignments) == 1:
        one_alignment_full_binding_sequences[key] = value
        valid_key_dictionary.update({key:len(potential_alignments)})
        print ("found one alignment sequence")
    elif len(potential_alignments) == 2:
        two_alignment_full_binding_sequences[key] = value
        valid_key_dictionary.update({key:len(potential_alignments)})
        print ("found two alignment sequence")
    else:
        mult_alignment_full_binding_sequences[key] = value
    
    #write all alignments for a given sequence to file    
    
    '''Open file: Creates a new file to store all alignments for a given sequence. These
                  files are stored in a PLISH desktop folder. The desktop file path is 
                  given as a global variable: desktop_path'''
    text_file = open(path + "/" + "Sequence_" + str(counter) + "_alignment" + ".txt","w+")
    
    #for each sequence, write general information to alignments file 
    text_file.write("Sequence " + str(counter) + " in potential_40bp_sequences\n")
    text_file.write("Sequence: " + value[0] + "\n")
    text_file.write("Start_index: " + str(key) + "\n")

    #give information about each sequence alignment
    for alignment in potential_alignments:
        for hsp in alignment.hsps:
             if hsp.expect < E_value_threshold:     #should I test E-value again?
                     text_file.write('****Alignment****\n')
                     text_file.write('sequence:' + str(alignment.title) + "\n")
                     text_file.write('length:' + str(alignment.length) + "\n")
                     text_file.write('e value:' + str(hsp.expect) + "\n")
                     text_file.write(hsp.query[0:75] + '...\n')
                     text_file.write(hsp.match[0:75] + '...\n')
                     text_file.write(hsp.sbjct[0:75] + '...\n')
    
    #close the text file
    text_file.close()
    
'''Write summary file: This file summarizes the PLISH alignments and returns a list
   of all plish probes with 0, 1 or 2 alignments. This summary file is stored
   on the desktop, whose path is given by the global variable desktop_path, in the PLISH folder.
   NEW: The file reports the sequence of each probe, the number of alignments, and the mRNA alignment
   indice of each probe sequence. This mRNA alignment indice corresponds to the base pair indice from 
   the 3' end of the mRNA of the gene of interest where the probe would align.'''
def produceSummaryFile():
    
    #opens the summary file
    summary_file = open(path + "/alignment_summary.txt", "w+")
    
    #write title of summary file: Alignment Summary
    summary_file.write("Alignment Summary: \n\n")
    
    #Title and print all sequences with no alignments
    summary_file.write("Sequences with no alignment: \n")
    summary_file.write("Sequence" + 40*" " + "Sequence Number\n")
    for key, value in zero_alignment_full_binding_sequences.iteritems():
        summary_file.write(value[0] + 8*" " + str(value[1]) + "\n")
    
    #Title and print all sequences with one alignment to output file
    summary_file.write("\nSequences with one alignment: \n")
    summary_file.write("Sequence" + 40*" " + "Sequence Number\n")
    for key, value in one_alignment_full_binding_sequences.iteritems():
        summary_file.write(value[0] + 8*" " + str(value[1]) + "\n")
    
    #print all sequences with two alignments to output file
    summary_file.write("\nSequences with two alignments: \n")
    summary_file.write("Sequence" + 40*" " + "Sequence Number\n")
    for key, value in two_alignment_full_binding_sequences.iteritems():
        summary_file.write(value[0] + 8*" " + str(value[1]) + "\n")
    summary_file.close()

'''find the reverse complements of each sequence with 0, 1, or 2 overlaps. This allows us
to find the complementary hprobe base pairs for a given RNA segment. '''
def findReverseComplements():
    
    for key, value in zero_alignment_full_binding_sequences.iteritems():
        value[0] = reverseComplement(value[0])
    
    for key, value in one_alignment_full_binding_sequences.iteritems():
        value[0] = reverseComplement(value[0])
    
    for key, value in two_alignment_full_binding_sequences.iteritems():
        value[0] = reverseComplement(value[0])

'''Create H_probe objects: Store potential h_probes as h_probe objects
   with member variables defined by the H_probe pair class. Two h_probes, 
   left and right h_probe, are stored in the container class hprobe_pair,
   which stores the left and right halves of an h_prober object.'''
def create_H_Probes():
    for i in range(0, 3):
        if i == 0:
            list = zero_alignment_full_binding_sequences.iteritems() #create hprobes with no alignments
            num_alignments = 0
        elif i == 1:
            list = one_alignment_full_binding_sequences.iteritems() #create hprobes with one alignment
            num_alignments = 1
        else:
            list = two_alignment_full_binding_sequences.iteritems() #create hprobes with two alignments
            num_alignments = 2
        
        #CHECK: if reverseComplement hprobe is necessary (currently commented out)
        #create h_probe objects corresponding to the left and right halves of an h_probe
        for key, value in list:
            print("new hprobe produced")
            '''EDIT'''
            left_hprobe_binding_sequence = value[0][0:alignment_probe_length]  #get corresponding mrna sequence (left probe)
            #left_hprobe_sequence = reverseComplement(left_hprobe_sequence)  #get complementary hprobe sequence (left probe)
            left_hprobe_binding_sequence = reverseString(left_hprobe_binding_sequence)  #reverse directionality (left probe)
            right_hprobe_binding_sequence = value[0][alignment_probe_length+2:alignment_sequence_length]
            #left_hprobe_sequence = reverseComplement(left_hprobe_sequence)
            right_hprobe_binding_sequence = reverseString(right_hprobe_binding_sequence)
            
            #measure and store the percent of 'G' or 'C' nucleaotides in each h_probe binding sequence
            percent_GC_left = percentGC(left_hprobe_binding_sequence)
            percent_GC_right = percentGC(right_hprobe_binding_sequence) 
            
            '''Generate full hprobe sequences ordered 5' to 3' (consisting of binding sequence, initiator sequence, and spacer). 
            These full hprobe are used in validity testing, including percent GC, end_annealment, etc. 
            Note that probe initiator sequence and spacer sequences are capitalized for readability.'''
            left_hprobe_sequence = initiator_sequence[0:len(initiator_sequence)/2].lower() + spacer_sequence.lower() + left_hprobe_binding_sequence
            right_hprobe_sequence = right_hprobe_binding_sequence + spacer_sequence.lower() + initiator_sequence[len(initiator_sequence)/2:len(initiator_sequence)].lower()
            
            #left_hprobe_sequence = left_hprobe_sequence.upper()     #uppercase entire sequence
            #right_hprobe_sequence = right_hprobe_sequence.upper()   #uppercase entire sequence
            
            #get the start indices (from 3' end of mrna) of each h_probe binding region
            left_index = key
            '''EDIT'''
            right_index = key + alignment_probe_length + 2      #add two since there is a spacing of two base pairs between left and right probes of a pair
            
            #get the sequence number for the h_probe
            sequence_number = value[1]
            
            #calculate hairpin delta G values for each hprobe
            hairpinDG_left = calcHairpinDG(left_hprobe_sequence)
            hairpinDG_right = calcHairpinDG(right_hprobe_sequence)
            
            #calculate homodimer delta G values for each hprobe
            homodimerDG_left = calcHomodimerDG(left_hprobe_sequence)
            homodimerDG_right = calcHomodimerDG(right_hprobe_sequence)
            
            #calculate heterodimer delta G values
            heterodimerDG = calcHeterodimerDG(left_hprobe_sequence, right_hprobe_sequence)
            
            '''record and store the melt_temp for each h_probe. Ensure that the melt_temp equation
            is correct. Also verify that the global variables used to calculate the melt_temp are set 
            correctly.'''
            melt_temp_left = get_melt_temp(salt_conc, percent_GC_left, left_hprobe_sequence, percent_formamide)
            melt_temp_right = get_melt_temp(salt_conc, percent_GC_right, right_hprobe_sequence, percent_formamide)
            
            '''Record and store the end annealment values for each h_probe. Be sure that the global
            variable for maximum end annealment (max_end_annealment) is set.'''
            end_annealment_left = end_annealment_count(left_hprobe_sequence) 
            end_annealment_right = end_annealment_count(right_hprobe_sequence)     
                   
            '''Create the left and right h_probe objects once all member variables have been 
            set.'''          
            left_hprobe = Hprobe(True, left_hprobe_sequence, left_index, left_hprobe_binding_sequence, percent_GC_left,  melt_temp_left, end_annealment_left, hairpinDG_left, homodimerDG_left)    #create left probe
            right_hprobe = Hprobe(False, right_hprobe_sequence, right_index, right_hprobe_binding_sequence, percent_GC_right, melt_temp_right, end_annealment_right, hairpinDG_right, homodimerDG_right)    #create right probe
            
            #add left and right half h_probes to Hprobe_pair object, which stores the complete h_probe
            hprobe_pair = Hprobe_pair(left_hprobe, right_hprobe, num_alignments, sequence_number, heterodimerDG)
            potential_Hprobe_pairs.append(hprobe_pair)
            
'''Print Function: This function prints the results to an excel file, which displays 
                   which potential probes are most robust, as well as which have potential
                   concerns. Key is included in the excel file. Probes that pass with no 
                   concerns are green, those with one concern are yellow, and those with two
                   concerns are red.'''
def printToExcel():
    
    #Output probe metrics to console
    probeMetrics(potential_Hprobe_pairs)
    
    #print hairpinDG_min for debugging
    print(hairpinDG_min)
    print(hairpinDG_min)
    print(max_percent_GC)
    
    print("Length of potential_Hprobe_pairs in printToExcel: " + str(len(potential_Hprobe_pairs)))
    for hprobe_pair in potential_Hprobe_pairs:
        #print("printToExcel working")
        if (valid_percentGC(hprobe_pair)):
            print (hprobe_pair.left.complete_sequence + " " + hprobe_pair.right.complete_sequence + "\n")
    
    #initialize workbook and worksheet to store PLISH data
    workBook, worksheet1, worksheet2 = initializeWorkSheet()
    plotProbeData(potential_Hprobe_pairs, worksheet1, worksheet2, workBook)

#create PLISH directory for program output
if not os.path.exists(path):
    try:
        os.makedirs(path)
    except OSError:  
        print ("Creation of the directory %s failed" % path)
    else:  
        print ("Successfully created the directory %s " % path)
'''if __name__ == "__main__":
    mrna_sequence, record_ID = read_mrna_sequence()
    runBLAST(3, 'AAAAAAAA', 4, record_ID, '/Users/tompritsky/Desktop/BLAST_practice_files/two_alignment.fasta')
    sortBLASTOutput(3, 'AAAAAAAA', 4, record_ID)
    produceSummaryFile()
'''   
    

