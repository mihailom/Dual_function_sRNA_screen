
#get RNA sequences from a genome .txt file given sRNA features from csv file
#sRNA_coords.csv contains RNA name (column 1), coordinates (start/end columns 2/3), orientation fwd/rev (column 4)
#def get_sequences('sRNA_coords.csv', "ecoligenome_MKM.txt", "sRNA_seqs.csv"):
def get_sequences(csv_file_name, genome_file_name, output_filename):

    sRNAs = []
    starts = []
    ends = []
    orientations = []
    import csv
    with open(csv_file_name,'r') as csvfile:
        read_file = csv.reader(csvfile, delimiter = ",")
        for row in read_file:
            sRNAs.append(row[0])
            orientations.append(row[3])
            starts.append(row[1])
            ends.append(row[2])
    csvfile.close()
    print(starts)
    print(ends)

    #import genome 
    genome = ""
    position = 0
    location = 0
    linecount = 0
    seq = ""
    revcomp = ""
    revseq = ""


    file = open(genome_file_name,"r")
    for raw_line in file:    #for each line
        line = raw_line.rstrip("\r\n")

        if linecount == 0:
            name = line.split(' ',1)[0]# save the name of the genome
            for char in line:
                if position > (int(len(name))+1): 
                    genome = genome + char
                    location +=1
                else:
                    genome = genome
                position+=1
        linecount += 1

        if linecount > 0:
            for char in line:
                genome = genome + char
                location +=1
        linecount += 1
    file.close()

    #get sequence for each sRNA
    sequences = []
    #there are some empty rows for some reason
    for n in range(1,97):
        start = starts[n]
        end = ends[n]
        orientation = orientations[n]

        tcomp = "A"
        ccomp = "G"
        acomp = "T"
        gcomp = "C"

        seq=genome[int(start):int(end)+1]
        
        if orientation == "rev":
            revseq = seq[::-1]
            revcomp = ""
            for nuc in revseq:
                if nuc == "A":
                    revcomp = revcomp + acomp
                if nuc == "T":
                    revcomp= revcomp + tcomp
                if nuc == "C":
                    revcomp = revcomp + ccomp
                if nuc == "G":
                    revcomp= revcomp + gcomp
            sequence= revcomp
        
        else:
            # sequence=genome[int(start):int(end)]
            sequence=seq
        sequences.append(sequence)

    with open(output_filename,"w") as seqs:
        writer = csv.writer(seqs)
        writer.writerows([sequences])
    seqs.close()

