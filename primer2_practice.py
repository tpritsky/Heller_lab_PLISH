#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 17:03:51 2019

@author: tompritsky
"""


import primer3                    #use this library to calculate thermodynamic values

#please note that one of the base pair sequences must be less than <=60bp in length. This bound is imposed by primer3

bp_sequence1 = 'ATTTCGCGCGCGTTTTTTTTT'
bp_sequence2 = 'ATTTGCTGGTCTTTTTTTTTT'

hairpin_parameters = primer3.calcHairpin(bp_sequence1)
print(hairpin_parameters.tm, hairpin_parameters.dg)
#primer3.calcTm()

#return hairpin dG from base pair sequence. Can also return melting temperature if needed
def calcHairpinDG(bp_sequence):
    hairpin_parameters = primer3.calcHairpin(bp_sequence)
    return hairpin_parameters.dg

#return homodimer dG from base pair sequence. Can also return melting temperature if needed
def calcHomodimerDG(bp_sequence):
    homodimer_parameters = primer3.calcHomodimer(bp_sequence)
    return homodimer_parameters.dg

#return heterodimer dG from a single base pair sequences. Can also return melting temperature if needed
def calcHeterodimerDG(sequence1, sequence2):
    heterodimer_parameters = primer3.calcHeterodimer(sequence1, sequence2)
    return heterodimer_parameters.dg
    
#Calculate the 3' end stability of DNA sequence 1 against DNA sequence 2
def calcEndStability(bp_sequence1, bp_sequence2):
    end_stability_parameters = primer3.calcEndStability(bp_sequence1, bp_sequence2)
    return end_stability_parameters.dg

hairpin_sequence1_dg = calcHairpinDG(bp_sequence1)
homodimer_sequence1_dg = calcHomodimerDG(bp_sequence1)
heterodimer_dg = calcHeterodimerDG('ATTTCGCGCGCGTTTTTTTTT', 'ATTTCGCGCGCGTTTTTTTTT')
#stability_dg = calcEndStability(bp_sequence1, bp_sequence2)

print ("For sequence 1, the hairpin delta g value is " + str(hairpin_sequence1_dg) + " and the homodimer dg value is " + str(homodimer_sequence1_dg))
print("For sequences 1 and 2, the heterodimer dg value is " + str(heterodimer_dg))
