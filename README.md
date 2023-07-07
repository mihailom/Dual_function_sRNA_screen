# Dual_function_sRNA_screen

This repo is meant for compiling likely small proteins (based on translation initiation rate cutoffs) encoded within regulatory RNA sequences (eg bacterial sRNAs) that may indicate dual-functionality of sRNAs. 
These 2 functions should be used in tandem as such: 
1. extract RNA sequences of interest from a genome file given coordinates, orientation using [get_sequences.py](https://github.com/mihailom/Dual_function_sRNA_screen/blob/main/get_sequences.py)
2. run sequences of interest in the Salis Lab RBS calculator to create matrix of quantitative Translation Initiation Rate (TIR) for every represented nucleotide position for which there was an RBS for each sRNA under investigation (rows) in Excel (direct output of RBS Calculator)
3. Pass TIR matrix to [DF_sRNA_candidate_screen.py](https://github.com/mihailom/Dual_function_sRNA_screen/blob/main/DF_sRNA_candidate_screen.py) to compile dual function candidates along with relevant RNA+protein sequence info
