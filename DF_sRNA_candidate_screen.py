#Mia Mihailovic 11-21-22 Dual-Function sRNA Filtering
#This code extracts dual-function sRNA candidates from a dataframe of Salis Lab RBS calculator results for a library of sRNAs
# and outputs relevant information including: sequences (DNA and protein) corresponding to the putative small peptide as well as 
# sequences useful for experimental design (eg. sRNA sequence (no stop) to be fused to SPA tag for Western)

import pandas as pd
import numpy as np
from Bio.Seq import Seq
#input = sRNA_TIRs, matrix of quantitative Translation Initiation Rate (TIR) for every represented nucleotide position for which there was an 
# RBS (columns, beginning at column 4) for each sRNA under investigation (rows) in Excel (direct output of RBS Calculator); 
# columns 2 and 3 are the sRNA sequence and extended sRNA sequence (to fit minimim length constraint of RBS Calculator), respectively

#def DF_sRNA_candidate_screen("sRNA_all_TIR.csv", "DF_sRNA_candidates")
def DF_sRNA_candidate_screen(TIR_csv, output_DF_candidates_csv):

    sRNA_TIRs = pd.read_csv(TIR_csv, header= None)
    sRNAs=[]
    sequences=[]
    passed_TIR=[]
    positions=[]
    passed_TIR_index=[]
    TIR=[]
    sRNA_seq=[]
    DF_DNA_seq = []
    DF_protein_seq = []
    sRNA_seq_from_start_codon=[]
    DF_sRNA_from_transc_start_no_stop=[]
    DF_sRNA_plus_upstream_no_stop=[]

    #loop through all positions of quantifiable TIR to extract putative peptide(s) for sRNAs
    for s, row in sRNA_TIRs.iterrows():
    sRNA_seq=sRNA_TIRs.iloc[s,1] 
    len_sRNA=int(len(sRNA_seq))
    #    print(len_sRNA)
    # print(pd.isna(sRNA_TIRs.iloc[1,39]))
    for n in range(3,204): #can't do till len_sRNA because not all positions are represented in the output csv.. had to hard-code
        if s>0: #exclude header
            if pd.isna(sRNA_TIRs.iloc[s,n]) =="True":
                print("not a number")
            else:
                if pd.to_numeric(sRNA_TIRs.iloc[s,n]) > 26: #empirical translation initiation rate cutoff based on positive hits to date
                    sequences.append(sRNA_seq)
                    sRNAs.append(sRNA_TIRs.iloc[s,0])
                    passed_TIR = sRNA_TIRs.iloc[0,n] #get relevant column name
                    passed_TIR_index.append(passed_TIR)
                    TIR.append(sRNA_TIRs.iloc[s,n])
                    position=''.join(filter(str.isdigit,passed_TIR)) #strip column name down to digits to get the right position eg pos 17 TIR --> 17
                    positions.append(position)
                    sRNA_seq_from_start_codon=str(sRNA_seq[int(position):int(len_sRNA)])
                    protein=Seq(sRNA_seq_from_start_codon).translate(to_stop=True) #translate the peptide from the putative start using Bio.Seq
                    len_prot=int(len(protein))
                    DF_sRNA_from_transc_start_no_stop = str(sRNA_seq[0:int(position)+(len_prot*3)]) #exclude stop
                    DF_sRNA_plus_upstream_no_stop.append(DF_sRNA_from_transc_start_no_stop)
                    DF_sRNA_start_to_stop_codon_no_stop = str(sRNA_seq[int(position):int(position)+(len_prot*3)]) #exclude stop
                    DF_DNA_seq.append(DF_sRNA_start_to_stop_codon_no_stop)
                    DF_protein_seq.append(protein)
    #output = details corresponding to each putative small peptide encoded by a sRNA (can be more than one for each sRNA): 
    import csv
    with open(output_DF_candidates_csv,"w") as out:
        writer = csv.writer(out)
        writer.writerows(zip(sRNAs,sequences, passed_TIR_index, positions, TIR, DF_DNA_seq, DF_sRNA_plus_upstream_no_stop, DF_protein_seq))
    out.close() 

