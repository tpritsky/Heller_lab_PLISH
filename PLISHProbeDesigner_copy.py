#! /usr/bin/env python
#  -*- coding: utf-8 -*-
#
# GUI module generated by PAGE version 4.14
# In conjunction with Tcl version 8.6
#    Jul 23, 2018 03:25:53 AM
 

#NOTE:
'''1. Changed all top to self.top'''

import sys
import Tkinter as tk
import tkMessageBox
import ttk
import plishcode_copy
import sysconfig

#Boolean indicating whether submit has been toggled
submit_in_progress = False

try:
    from Tkinter import *
except ImportError:
    from tkinter import *

try:
    import ttk
    py3 = False
except ImportError:
    import tkinter.ttk as ttk
    py3 = True

import PLISHDesigner_support_copy

#TEST: increments the value in a label
counter = 0
def inc_label(label):
      def count():
        global counter
        counter += 1
        label['text'] = str(counter)
        label.after(30, count)
      count()

'''reinitialize_sequence_number_label: reinitializes the text in a label to the 
sequence number currently being aligned. Takes in three inputs: 1. the label object
to be altered, 2. a list of keys(corresponding to sequence locations) and a list of values
(corresponding to nucleotide sequence strings)'''
#def reinitialize_sequence_number_label(label, key, value):

#reinitializes the text in a textEntry field to the string in value
i = 0
def reinitialize_txtEntry(txtEntry, key, value):
    def count():
        global i
        txtEntry.delete(0, END)
        if(i != len(value)):
            txtEntry.insert(0, value[i])
            plishcode_copy.runBLAST(key[i], value[i], i, plishcode_copy.record_ID)
            i += 1
            txtEntry.after(30, count)
    count()
  
'''reinitialize_sequence_label: reinitializes the text in a label corresponding
to the sequence being aligned (stored in value). Takes in three inputs: 1. the label object
to be altered, 2. a list of keys(corresponding to sequence locations) and 3. a list of values
(corresponding to nucleotide sequence strings)'''
i = 0
def reinitialize_sequence_label(label, key, value):
    def count():
        global i
        #if(i <= 10):
        label['text'] = value[i]
        plishcode_copy.runBLAST(key[i], value[i], i, plishcode_copy.record_ID)
        i += 1
        label.after(30, count)
    count()

'''TEST: reinitialize_sequence_label: reinitializes the text in a label corresponding
to the sequence being aligned (stored in value). Takes in three inputs: 1. the label object
to be altered, 2. a list of keys(corresponding to sequence locations) and 3. a list of values
(corresponding to nucleotide sequence strings)'''
i = 0
progressBarValue = 0
def reinitialize_minGC_label(root, currentSequence, sequenceNumber, progressBar, key, value, finalize):
    progressBarValue = 0
    def count():
        global i
        if(i == len(value)):
            finalize()
            return
        if(i < len(value)):
            currentSequence['text'] = value[i]
            sequenceNumber['text'] = ("Sequence " + str(i + 1) + " of " + str(len(value)))
            progressBarValue = 100
            print("PROG BAR VAL")
            print(progressBarValue)
            plishcode_copy.runBLAST(key[i], value[i], i, plishcode_copy.record_ID)
            i += 1
            root.after(30, count)
        print("\nSequences with no alignments:")
    count()

#run single BLAST alignment
def runSingleBlast(key, value, txtEntry):
    #run BLAST alignment on potential probe sequences
    if(reinitialize_txtEntry(txtEntry, value)):
        plishcode_copy.runBLAST(key, value, plishcode_copy.counter, plishcode_copy.record_ID)
        plishcode_copy.counter+=1
    
def vp_start_gui():
    '''Starting point when module is the main routine.'''
    global val, w, root
    root = Tk()
    top = PLISH_Designer (root)
    #PLISHDesigner_support.init(root, top)
    root.mainloop()
    w = None

def create_PLISH_Designer(root, *args, **kwargs):
    '''Starting point when module is imported by another program.'''
    global w, w_win, rt
    rt = root
    w = Toplevel (root)
    top = PLISH_Designer (w)
    #PLISHDesigner_support.init(w, top, *args, **kwargs)
    return (w, top)

def destroy_PLISH_Designer():
    global w
    w.destroy()
    w = None

class PLISH_Designer:
    #Default Values
    SaltConcDefault = ".5"
    percentFormamideDefault = ".2"
    MinGCDefault = ".5"
    MaxGCDefault = ".65"
    EThreshDefault = ".9"
    MaxEndAnnealDefault = "2"
    FASTAPathDefault = "C:/Desktop/PLISH/FASTAPath"
    DBPathDefault = "C:/Desktop/PLISH/DBPath"
    DesktopPathDefault = "C:Users/Mirko/Desktop"
    GeneValueDefault= "Otoferlin"
    top = None
    counter = 0
    
    def reinitialize_txtSaltConc(self, value):
        self.txtSaltConc.delete(0, END)
        self.txtSaltConc.insert(0, value)
    def reinitialize_txtpercentFormamide(self, value):
        self.txtpercentFormamide.delete(0, END)
        self.txtpercentFormamide.insert(0, value)
    def reinitialize_txtMinGC(self, value):
        self.txtMinGC.delete(0, END)
        self.txtMinGC.insert(0, value)
    def reinitialize_txtMaxGC(self, value):
        self.txtMaxGC.delete(0, END)
        self.txtMaxGC.insert(0, value)
    def reinitialize_txtEThresh(self, value):
        self.txtEThresh.delete(0, END)
        self.txtEThresh.insert(0, value)
    def reinitialize_txtMaxEndAnneal(self, value):
        self.txtMaxEndAnneal.delete(0, END)
        self.txtMaxEndAnneal.insert(0, value)
    
    '''Need to check if this applies to lists
    def reinitialize_listboxSpecies(self, value):
        self.listboxSpecies.delete(0, END)
        self.listboxSpecies.insert(0, value)'''
    def reinitialize_txtGene(self, value):
        self.txtGene.delete(0, END)
        self.txtGene.insert(0, value)
    def reinitialize_txtFASTAPath(self, value):
        self.txtFASTAPath.delete(0, END)
        self.txtFASTAPath.insert(0, value)
    def reinitialize_txtDBPath(self, value):
        self.txtDBPath.delete(0, END)
        self.txtDBPath.insert(0, value)
    def reinitialize_txtCurrentSeq(self, value):
        self.txtCurrentSeq.delete(0, END)
        self.txtCurrentSeq.insert(0, value)
    def reinitialize_txtCurrentSeq(self, value):
        self.txtCurrentSeq.delete(0, END)
        self.txtCurrentSeq.insert(0, value)
    def reinitialize_txtSeqNum(self, value):
        self.txtSeqNum.delete(0, END)
        self.txtSeqNum.insert(0, value)
    def reinitialize_txtDesktopPath(self, value):
        self.txtDesktopPath.delete(0, END)
        self.txtDesktopPath.insert(0, value)
    '''Need to check if this also works for progress bar
    def reinitialize_TProgressbar1(self, value):
        self.txtSeqNum.delete(0, END)
        self.txtSeqNum.insert(0, value)'''
    
    def setDefaultEntries(self):
        self.reinitialize_txtSaltConc(self.SaltConcDefault)
        self.reinitialize_txtpercentFormamide(self.percentFormamideDefault)
        self.reinitialize_txtMinGC(self.MinGCDefault)
        self.reinitialize_txtMaxGC(self.MaxGCDefault)
        self.reinitialize_txtEThresh(self.EThreshDefault)
        self.reinitialize_txtMaxEndAnneal(self.MaxEndAnnealDefault)
        self.reinitialize_txtGene(self.GeneValueDefault)
        self.reinitialize_txtFASTAPath(self.FASTAPathDefault)
        self.reinitialize_txtDBPath(self.DBPathDefault)
        self.reinitialize_txtDesktopPath(self.DesktopPathDefault)
        
    '''Set User Parameters: This function reinitializes all parameters to user set values.
    If user did not provide a parameter value, the parameters are set to default values.
    Called by submitUserInput, '''
    
    def setUserParams(self):
        #species
        if(self.combolistSpecies.get()):
            plishcode_copy.species_name = self.combolistSpecies.get()
            print("SPECIES")
            print(plishcode.species_name)
        #gene
        if(self.txtGene.get()):
            plishcode_copy.gene_name = self.txtGene.get()
            print("GENE")
            print(plishcode_copy.gene_name)
        #Path to FASTA file for mrna of gene of interest 
        if(self.txtFASTAPath.get()):
            plishcode_copy.FASTA_file = self.txtFASTAPath.get()
            print("GENE")
            print(plishcode_copy.FASTA_file)
        #desktop path
        if(self.txtDesktopPath.get()):
            plishcode_copy.desktop_path = self.txtDesktopPath.get()
            print("GENE")
            print(plishcode_copy.desktop_path)
        #Nucleotide database
        if(self.txtDBPath.get()):
            plishcode_copy.nucleotideDatabase= self.txtDBPath.get()
        #BLAST sequence alignment thresholds (please set):
        if(self.txtEThresh.get()):
            plishcode_copy.E_value_threshold = self.txtEThresh.get() 
        #Salt concentration values:
        if(self.txtSaltConc.get()):
            plishcode_copy.salt_conc = self.txtSaltConc.get()
        if(self.txtpercentFormamide.get()):
            plishcode_copy.percent_formamide = self.txtpercentFormamide.get()
        #Min percentGC value:
        if(self.txtMinGC.get()):
            plishcode_copy.min_percent_GC = self.txtMinGC.get()     #minimum percentage of GC in a single Hprobe
        #Max percentGC value
        if(self.txtMaxGC.get()):
            plishcode_copy.max_percent_GC = self.txtMaxGC.get()     #maximum percentage of GC in a single Hprobe'''
        #end annealment values:
        if(self.txtMaxEndAnneal.get()):
            plishcode_copy.max_end_annealment = self.txtMaxEndAnneal.get()     #maximum number of overlaps between start and end of h_probe sequence 
    
    def submitUserInput(self):
        self.setUserParams()
        global submit_in_progress
        if(submit_in_progress == True):
            return
        else:
            submit_in_progress = True
        def finalize():
            #runSingleBlast(key, value, self.txtGene)
            '''WHY DOESN'T self.txtGene reinitialize more than once?'''
                
            #print all sequences with no alignments. CURRENTLY just a test
            print("\nSequences with no alignments:*")
            for key, value in plishcode_copy.zero_alignment_full_binding_sequences.iteritems():
                print(value[1])
            
            '''Write summary file: This file summarizes the PLISH alignments and returns a list
               a list of all plish probes with 0, 1 or 2 alignments. This summary file is stored
               on the desktop, whose path is given by the global variable desktop_path.'''
            #TEST
            plishcode_copy.produceSummaryFile()
            
            #find the reverse complements of each sequence with 0, 1, or 2 overlaps. This allows us
            #to find the complementary hprobe base pairs for a given RNA segment. 
            plishcode_copy.findReverseComplements()
            
            '''Create H_probe objects: Store potential h_probes as h_probe objects
               with member variables defined by the H_probe pair class.'''
            plishcode_copy.create_H_Probes()
            
            '''Print PLISH probe information to an excel file'''
            plishcode_copy.printToExcel()
            print("EXITING")
            
            #this allows another submission after a run
            global submit_in_progress
            submit_in_progress = False
        
        '''reads FASTA file and returns string containing mRNA sequence for alignment,
        as wellas record ID fromo= original FASTA file'''
        plishcode_copy.mrna_sequence, plishcode_copy.record_ID = plishcode_copy.read_mrna_sequence()
        
        #extract all potential 40BP sequences from target mRNA sequence, which contain 'TA' or 'AG' in center
        plishcode_copy.potential_40bp_sequences = plishcode_copy.extractSequences(plishcode_copy.mrna_sequence)
        
        '''THIS IS WHERE I AM CONFUSED:
            BLAST code: This script runs BLAST, which can align potential hprobe sequences to mrna. Returns all
            alignments above a minimal threshold value. Still don't know why a for loop doesn't work with after.'''
        plishcode_copy.counter = 0 #store sequence number
        
        #reinitialize i
        global i
        i = 0
        key_list = []
        value_list = []
        global progressBarValue
        progressBarValue = 0
        for key, value in plishcode_copy.potential_40bp_sequences.iteritems():
            key_list.append(key)
            value_list.append(value)
            #progressBarValue += 1
        
        #Runs BLAST alignment locally and displays current sequence being aligned to GUI
        reinitialize_minGC_label(self.top, self.txtCurrentSeq, self.txtSeqNum, self.TProgressbar1, key_list, value_list, finalize)
        #reinitialize_sequence_label(self.txtCurrentSeq, key_list, value_list)

        
    def Cancel(self):
      msg = tkMessageBox.askyesno("PLISH Designer", "Are you sure you want to cancel run? Probe Designer will exit.")
      if (msg):
          self.top.destroy()  
          #self.top.quit()
    
    #master=top before
    def __init__(self, top1):
        self.top = top1
        '''This class configures and populates the toplevel window.
           top is the toplevel containing window.'''
        _bgcolor = '#d9d9d9'  # X11 color: 'gray85'
        _fgcolor = '#000000'  # X11 color: 'black'
        _compcolor = '#d9d9d9' # X11 color: 'gray85'
        _ana1color = '#d9d9d9' # X11 color: 'gray85' 
        _ana2color = '#d9d9d9' # X11 color: 'gray85' 
        self.style = ttk.Style()
        if sys.platform == "win32":
            self.style.theme_use('winnative')
        self.style.configure('.',background=_bgcolor)
        self.style.configure('.',foreground=_fgcolor)
        self.style.map('.',background=
            [('selected', _compcolor), ('active',_ana2color)])

        self.top.geometry("841x687+657+122")
        self.top.title("PLISH Designer")
        self.top.configure(background="#d9d9d9")

        self.btnSubmit = Button(self.top)
        self.btnSubmit.place(relx=0.15, rely=0.79, height=33, width=89)
        self.btnSubmit.configure(activebackground="#d9d9d9")
        self.btnSubmit.configure(activeforeground="#000000")
        self.btnSubmit.configure(background="#d9d9d9")
        self.btnSubmit.configure(disabledforeground="#a3a3a3")
        self.btnSubmit.configure(foreground="#000000")
        self.btnSubmit.configure(highlightbackground="#d9d9d9")
        self.btnSubmit.configure(highlightcolor="black")
        self.btnSubmit.configure(pady="0")
        self.btnSubmit.configure(text='''Submit''')
        self.btnSubmit.configure(width=89)
        self.btnSubmit.configure(command=self.submitUserInput)

        self.btnCancel = Button(self.top)
        self.btnCancel.place(relx=0.62, rely=0.35, height=33, width=85)
        self.btnCancel.configure(activebackground="#d9d9d9")
        self.btnCancel.configure(activeforeground="#000000")
        self.btnCancel.configure(background="#d9d9d9")
        self.btnCancel.configure(disabledforeground="#a3a3a3")
        self.btnCancel.configure(foreground="#000000")
        self.btnCancel.configure(highlightbackground="#d9d9d9")
        self.btnCancel.configure(highlightcolor="black")
        self.btnCancel.configure(pady="0")
        self.btnCancel.configure(text='''Cancel Run''')
        #self.btnCancel.configure(width=89)
        self.btnCancel.configure(command=self.Cancel)

        self.Label1 = Label(self.top)
        self.Label1.place(relx=0.01, rely=0.01, height=26, width=244)
        self.Label1.configure(background="#d9d9d9")
        self.Label1.configure(disabledforeground="#a3a3a3")
        self.Label1.configure(foreground="#000000")
        self.Label1.configure(text='''Welcome to PLISH Probe Designer...''')

        self.Label2 = Label(self.top)
        self.Label2.place(relx=0.01, rely=0.1, height=26, width=121)
        self.Label2.configure(background="#d9d9d9")
        self.Label2.configure(disabledforeground="#a3a3a3")
        self.Label2.configure(foreground="#000000")
        self.Label2.configure(text='''Input Parameters:''')
        
        self.Label3 = Label(self.top)
        self.Label3.place(relx=0.01, rely=0.16, height=26, width=131)
        self.Label3.configure(background="#d9d9d9")
        self.Label3.configure(disabledforeground="#a3a3a3")
        self.Label3.configure(foreground="#000000")
        self.Label3.configure(text='''1. Probe Constants:''')

        self.Label4 = Label(self.top)
        self.Label4.place(relx=0.04, rely=0.19, height=26, width=224)
        self.Label4.configure(background="#d9d9d9")
        self.Label4.configure(disabledforeground="#a3a3a3")
        self.Label4.configure(foreground="#000000")
        self.Label4.configure(text='''a. Melt. Temperature Parameters:''')

        self.Label5 = Label(self.top)
        self.Label5.place(relx=0.06, rely=0.22, height=26, width=150)
        self.Label5.configure(background="#d9d9d9")
        self.Label5.configure(disabledforeground="#a3a3a3")
        self.Label5.configure(foreground="#000000")
        self.Label5.configure(text='''1. Salt Concentration:''')

        self.txtSaltConc = Entry(self.top)
        self.txtSaltConc.place(relx=0.3, rely=0.23,height=24, relwidth=0.1)
        self.txtSaltConc.configure(background="white")
        self.txtSaltConc.configure(disabledforeground="#a3a3a3")
        self.txtSaltConc.configure(font="TkFixedFont")
        self.txtSaltConc.configure(foreground="#000000")
        self.txtSaltConc.configure(insertbackground="black")
        self.txtSaltConc.configure(width=84)

        self.Label6 = Label(self.top)
        self.Label6.place(relx=0.06, rely=0.26, height=26, width=116)
        self.Label6.configure(background="#d9d9d9")
        self.Label6.configure(disabledforeground="#a3a3a3")
        self.Label6.configure(foreground="#000000")
        self.Label6.configure(text='''2. % Formamide:''')

        self.Label7 = Label(self.top)
        self.Label7.place(relx=0.04, rely=0.3, height=26, width=93)
        self.Label7.configure(background="#d9d9d9")
        self.Label7.configure(disabledforeground="#a3a3a3")
        self.Label7.configure(foreground="#000000")
        self.Label7.configure(text='''b. Thresholds''')

        self.txtpercentFormamide = Entry(self.top)
        self.txtpercentFormamide.place(relx=0.3, rely=0.27, height=24
                , relwidth=0.1)
        self.txtpercentFormamide.configure(background="white")
        self.txtpercentFormamide.configure(disabledforeground="#a3a3a3")
        self.txtpercentFormamide.configure(font="TkFixedFont")
        self.txtpercentFormamide.configure(foreground="#000000")
        self.txtpercentFormamide.configure(insertbackground="black")
        self.txtpercentFormamide.configure(width=84)

        self.Label10 = Label(self.top)
        self.Label10.place(relx=0.06, rely=0.34, height=26, width=128)
        self.Label10.configure(background="#d9d9d9")
        self.Label10.configure(disabledforeground="#a3a3a3")
        self.Label10.configure(foreground="#000000")
        self.Label10.configure(text='''1. Min. Percent GC:''')
        self.Label10.configure(width=128)

        self.Label11 = Label(self.top)
        self.Label11.place(relx=0.06, rely=0.38, height=26, width=131)
        self.Label11.configure(background="#d9d9d9")
        self.Label11.configure(disabledforeground="#a3a3a3")
        self.Label11.configure(foreground="#000000")
        self.Label11.configure(text='''2. Max. Percent GC:''')

        self.txtMinGC = Entry(self.top)
        self.txtMinGC.place(relx=0.3, rely=0.33,height=24, relwidth=0.1)
        self.txtMinGC.configure(background="white")
        self.txtMinGC.configure(disabledforeground="#a3a3a3")
        self.txtMinGC.configure(font="TkFixedFont")
        self.txtMinGC.configure(foreground="#000000")
        self.txtMinGC.configure(insertbackground="black")
        self.txtMinGC.configure(width=84)

        self.txtMaxGC = Entry(self.top)
        self.txtMaxGC.place(relx=0.3, rely=0.37,height=24, relwidth=0.1)
        self.txtMaxGC.configure(background="white")
        self.txtMaxGC.configure(disabledforeground="#a3a3a3")
        self.txtMaxGC.configure(font="TkFixedFont")
        self.txtMaxGC.configure(foreground="#000000")
        self.txtMaxGC.configure(insertbackground="black")
        self.txtMaxGC.configure(width=84)

        self.Label12 = Label(self.top)
        self.Label12.place(relx=0.01, rely=0.48, height=26, width=205)
        self.Label12.configure(background="#d9d9d9")
        self.Label12.configure(disabledforeground="#a3a3a3")
        self.Label12.configure(foreground="#000000")
        self.Label12.configure(text='''2. Alignment MRNA Sequence:''')

        self.Label13 = Label(self.top)
        self.Label13.place(relx=0.06, rely=0.52, height=26, width=118)
        self.Label13.configure(background="#d9d9d9")
        self.Label13.configure(disabledforeground="#a3a3a3")
        self.Label13.configure(foreground="#000000")
        self.Label13.configure(text='''1. Species Name:''')

        '''
        self.listboxSpecies = Listbox(self.top)
        self.listboxSpecies.place(relx=0.25, rely=0.52, relheight=0.04
                , relwidth=0.15)
        self.listboxSpecies.configure(background="white")
        self.listboxSpecies.configure(disabledforeground="#a3a3a3")
        self.listboxSpecies.configure(font="TkFixedFont")
        self.listboxSpecies.configure(foreground="#000000")
        self.listboxSpecies.configure(width=124)
        self.listboxSpecies.insert(END, "Species Name")
        for item in ["Gallus Gallus", "Mus Musculus"]:
            self.listboxSpecies.insert(END, item)
        '''

        self.combolistSpecies = ttk.Combobox(self.top)
        self.combolistSpecies.place(relx=0.25, rely=0.52, relheight=0.03
                , relwidth=0.15)
        self.combolistSpecies['values'] = ('Chicken', 'Mouse')

        
        self.Label8 = Label(self.top)
        self.Label8.place(relx=0.06, rely=0.56, height=26, width=102)
        self.Label8.configure(background="#d9d9d9")
        self.Label8.configure(disabledforeground="#a3a3a3")
        self.Label8.configure(foreground="#000000")
        self.Label8.configure(text='''2. Gene Name:''')

        self.Label9 = Label(self.top)
        self.Label9.place(relx=0.06, rely=0.6, height=26, width=143)
        self.Label9.configure(background="#d9d9d9")
        self.Label9.configure(disabledforeground="#a3a3a3")
        self.Label9.configure(foreground="#000000")
        self.Label9.configure(text='''3. Path to FASTA file:''')

        self.txtGene = Entry(self.top)
        self.txtGene.place(relx=0.25, rely=0.57,height=24, relwidth=0.15)
        self.txtGene.configure(background="white")
        self.txtGene.configure(disabledforeground="#a3a3a3")
        self.txtGene.configure(font="TkFixedFont")
        self.txtGene.configure(foreground="#000000")
        self.txtGene.configure(insertbackground="black")
        self.txtGene.configure(width=124)

        self.txtFASTAPath = Entry(self.top)
        self.txtFASTAPath.place(relx=0.25, rely=0.61,height=24, relwidth=0.15)
        self.txtFASTAPath.configure(background="white")
        self.txtFASTAPath.configure(disabledforeground="#a3a3a3")
        self.txtFASTAPath.configure(font="TkFixedFont")
        self.txtFASTAPath.configure(foreground="#000000")
        self.txtFASTAPath.configure(insertbackground="black")
        self.txtFASTAPath.configure(width=124)

        self.Label14 = Label(self.top)
        self.Label14.place(relx=0.01, rely=0.65, height=26, width=124)
        self.Label14.configure(background="#d9d9d9")
        self.Label14.configure(disabledforeground="#a3a3a3")
        self.Label14.configure(foreground="#000000")
        self.Label14.configure(text='''3. Input File Paths:''')

        self.Label15 = Label(self.top)
        self.Label15.place(relx=0.06, rely=0.68, height=26, width=147)
        self.Label15.configure(background="#d9d9d9")
        self.Label15.configure(disabledforeground="#a3a3a3")
        self.Label15.configure(foreground="#000000")
        self.Label15.configure(text='''1. Database File Path:''')

        self.menubar = Menu(self.top,font="TkMenuFont",bg=_bgcolor,fg=_fgcolor)
        self.top.configure(menu = self.menubar)
        
        self.Label16 = Label(self.top)
        self.Label16.place(relx=0.06, rely=0.72, height=26, width=139)
        self.Label16.configure(background="#d9d9d9")
        self.Label16.configure(disabledforeground="#a3a3a3")
        self.Label16.configure(foreground="#000000")
        self.Label16.configure(text='''2. Desktop File Path:''')

        self.txtDBPath = Entry(self.top)
        self.txtDBPath.place(relx=0.25, rely=0.68,height=24, relwidth=0.15)
        self.txtDBPath.configure(background="white")
        self.txtDBPath.configure(disabledforeground="#a3a3a3")
        self.txtDBPath.configure(font="TkFixedFont")
        self.txtDBPath.configure(foreground="#000000")
        self.txtDBPath.configure(insertbackground="black")
        self.txtDBPath.configure(width=124)

        self.txtDesktopPath = Entry(self.top)
        self.txtDesktopPath.place(relx=0.25, rely=0.72, height=24, relwidth=0.15)

        self.txtDesktopPath.configure(background="white")
        self.txtDesktopPath.configure(disabledforeground="#a3a3a3")
        self.txtDesktopPath.configure(font="TkFixedFont")
        self.txtDesktopPath.configure(foreground="#000000")
        self.txtDesktopPath.configure(insertbackground="black")
        self.txtDesktopPath.configure(width=124)
        
        global progressBarValue
        self.TProgressbar1 = ttk.Progressbar(self.top, mode = 'indeterminate')
        #progress= ttk.Progressbar(root, orient = 'horizontal', maximum = 10000, variable=downloaded, mode = 'determinate')

        self.TProgressbar1.place(relx=0.68, rely=0.27, relwidth=0.12
                , relheight=0.0, height=17)

        self.Label17 = Label(self.top)
        self.Label17.place(relx=0.06, rely=0.41, height=26, width=142)
        self.Label17.configure(background="#d9d9d9")
        self.Label17.configure(disabledforeground="#a3a3a3")
        self.Label17.configure(foreground="#000000")
        self.Label17.configure(text='''3. E-value Threshold:''')

        self.txtEThresh = Entry(self.top)
        self.txtEThresh.place(relx=0.3, rely=0.41,height=24, relwidth=0.1)
        self.txtEThresh.configure(background="white")
        self.txtEThresh.configure(disabledforeground="#a3a3a3")
        self.txtEThresh.configure(font="TkFixedFont")
        self.txtEThresh.configure(foreground="#000000")
        self.txtEThresh.configure(insertbackground="black")
        self.txtEThresh.configure(width=84)

        self.Label18 = Label(self.top)
        self.Label18.place(relx=0.54, rely=0.1, height=26, width=117)
        self.Label18.configure(background="#d9d9d9")
        self.Label18.configure(disabledforeground="#a3a3a3")
        self.Label18.configure(foreground="#000000")
        self.Label18.configure(text='''Progress Metrics:''')

        self.Label19 = Label(self.top)
        self.Label19.place(relx=0.55, rely=0.16, height=26, width=140)
        self.Label19.configure(background="#d9d9d9")
        self.Label19.configure(disabledforeground="#a3a3a3")
        self.Label19.configure(foreground="#000000")
        self.Label19.configure(text='''1. Current Sequence:''')

        self.Label20 = Label(self.top)
        self.Label20.place(relx=0.55, rely=0.23, height=25, width=150)
        self.Label20.configure(background="#d9d9d9")
        self.Label20.configure(disabledforeground="#a3a3a3")
        self.Label20.configure(foreground="#000000")
        self.Label20.configure(text='''2. Sequence Number:''')

        self.txtCurrentSeq = Label(self.top)
        self.txtCurrentSeq.place(relx=.54, rely=0.19, height=26, width=400)
        self.txtCurrentSeq.configure(background="#d9d9d9")
        self.txtCurrentSeq.configure(disabledforeground="#a3a3a3")
        self.txtCurrentSeq.configure(foreground="#000000")
        self.txtCurrentSeq.configure(text='''waiting for run...''')
        self.txtCurrentSeq.configure(font=("Courier", 13))
        self.txtCurrentSeq.configure(width= 210)

        self.txtSeqNum = Label(self.top)
        self.txtSeqNum.place(relx=0.74, rely=0.23, height=25, width=202)
        self.txtSeqNum.configure(background="#d9d9d9")
        self.txtSeqNum.configure(disabledforeground="#a3a3a3")
        self.txtSeqNum.configure(foreground="#000000")
        self.txtSeqNum.configure(highlightcolor="#646464")
        self.txtSeqNum.configure(width=202)
        self.txtSeqNum.configure(text='''Display Run''')

        self.Label23 = Label(self.top)
        self.Label23.place(relx=0.55, rely=0.27, height=17, width=90)
        self.Label23.configure(background="#d9d9d9")
        self.Label23.configure(disabledforeground="#a3a3a3")
        self.Label23.configure(foreground="#000000")
        self.Label23.configure(text='''3. Progress:''')
        self.Label23.configure(width=77)

        #Fix: currently empty, but displays runtime information
        self.Label24 = Label(self.top)
        self.Label24.place(relx=0.56, rely=0.44, height=126, width=202)
        self.Label24.configure(background="#d9d9d9")
        self.Label24.configure(disabledforeground="#a3a3a3")
        self.Label24.configure(foreground="#000000")
        self.Label24.configure(text='''''')
        self.Label24.configure(width=202)

        self.Label25 = Label(self.top)
        self.Label25.place(relx=0.06, rely=0.45, height=26, width=131)
        self.Label25.configure(background="#d9d9d9")
        self.Label25.configure(disabledforeground="#a3a3a3")
        self.Label25.configure(foreground="#000000")
        self.Label25.configure(text='''4. Max End Anneal:''')

        self.txtMaxEndAnneal = Entry(self.top)
        self.txtMaxEndAnneal.place(relx=0.3, rely=0.45,height=24, relwidth=0.1)
        self.txtMaxEndAnneal.configure(background="white")
        self.txtMaxEndAnneal.configure(disabledforeground="#a3a3a3")
        self.txtMaxEndAnneal.configure(font="TkFixedFont")
        self.txtMaxEndAnneal.configure(foreground="#000000")
        self.txtMaxEndAnneal.configure(insertbackground="black")
        self.txtMaxEndAnneal.configure(width=84)

        self.btnSetDefaults = Button(self.top)
        self.btnSetDefaults.place(relx=0.18, rely=0.1, height=33, width=92)
        self.btnSetDefaults.configure(activebackground="#d9d9d9")
        self.btnSetDefaults.configure(activeforeground="#000000")
        self.btnSetDefaults.configure(background="#d9d9d9")
        self.btnSetDefaults.configure(disabledforeground="#a3a3a3")
        self.btnSetDefaults.configure(foreground="#000000")
        self.btnSetDefaults.configure(highlightbackground="#d9d9d9")
        self.btnSetDefaults.configure(highlightcolor="black")
        self.btnSetDefaults.configure(pady="0")
        self.btnSetDefaults.configure(text='''Set Defaults''')
        self.btnSetDefaults.configure(command=self.setDefaultEntries)

if __name__ == '__main__':
    vp_start_gui()


