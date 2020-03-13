#!/usr/bin/env python3

#Fixa: Filen som man sätter in  måste sluta på fna!!!!!!!

import sys
usage='This program count the GC content i GC1, GC2 and GC3 for a fasta file, the length of the sequence and the codon usage.'


if len(sys.argv) != 4:
   print('This is a friendly reminder that there needs to be four arguments. The first argument is the name of the software (software_fasta.py) and the second argument needs to be a fasta file, the two final files are your names for the output files.', file=sys.stderr)
   sys.exit()

if not '.fna' in sys.argv[1]:
    print('This is a friendly reminder that the input file needs to be a fasta file and have the suffix .fna', file=sys.stderr)
    sys.exit()

#This is the amino acid codon table.
gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

def get_sequences():
    """This function collect all sequences and store them in a list. The function also work if the
    sequence is split into several lines). """
    string=''
    seq=[]
    with open (sys.argv[1], 'r') as fin:
        for line in fin:
            line=line.rstrip()
            if line:
                if not line.startswith('>'):
                    line=line.upper()
                    string +=(line)
                elif line.startswith('>'):
                    if string:
                        seq.append(string)
                        string = ''
        seq.append(string)
        return seq


def present_gc_content(seq):
    outfile=open (sys.argv[2],'w')
    """This function calculate the GC content in 1:st,2:nd and 3:rd position of the codon."""
    position_1=''
    position_2=''
    position_3=''
    for sequence in seq:
        position_1 +=sequence[::3]
        position_2 +=sequence[1::3]
        position_3 +=sequence[2::3]
    print('GC1\t{}'.format(round(calculate_GC(position_1), 3)))
    print('GC2\t{}'.format(round(calculate_GC(position_2), 3)))
    print('GC3\t{}'.format(round(calculate_GC(position_3), 3)))


def calculate_GC(position):
    """The function calculate and return the GC"""
    g = position.count('G')
    c = position.count('C')
    gc = (g+c)/len(position)
    return gc

def total_GC_content(seq):
    """The function calculate the overal GC-content in percent"""
    all_seq=''
    for sequence in seq:
        all_seq += sequence
    all_GC_content = calculate_GC(all_seq)
    print('GC_Percent\t{}'.format((round(all_GC_content*100, 3))))


def length(seq):
    """This function calculate the length of the sequence"""
    length=0
    number_of_seq=0
    for sequence in seq:
        number_of_seq +=1
        length +=len(sequence)
    length_seq = length/number_of_seq
    print('Length\t{}'.format(length_seq))


def create_codon_dict(seq):
    """The function creates a dictionary with key = codon and value = the number of that codon"""
    gen_code = ['TTT','TTC','TTA','TTG','TCT','TCC','TCA','TCG','TAT','TAC', 'TAA', 'TAG', 'TGA'
    ,'TGT','TGC','TGG','CTT','CTC','CTA','CTG','CCT','CCC','CCA','CCG',
    'CAT','CAC','CAA','CAG','CGT','CGC','CGA','CGG','ATT','ATC','ATA','ATG','ACT',
    'ACC','ACA','ACG','AAT','AAC','AAA','AAG','AGT','AGC','AGA','AGG','GTT','GTC',
    'GTA','GTG','GCT','GCC','GCA','GCG','GAT','GAC','GAA','GAG','GGT','GGC','GGA','GGG']
    stop_codons = ['TAA', 'TGA', 'TAG']
    codon_dict={}
    for codon in gen_code:
        if not codon in stop_codons:
            codon_dict[codon] = 0
    for sequence in seq:
        for i in range(0,len(sequence),3):
            codon = sequence[i:i+3]
            if not len(codon)%3 == 0:
                continue
            else:
                if not codon in stop_codons:
                    codon_dict[codon] +=1
    return codon_dict


def sort_by_value(codon_dict):
    """The function sort the dictionary by value"""
    outfile=open (sys.argv[2],'w')
    sortbyvalue=dict(sorted(codon_dict.items(), key=lambda t: t[1]))
    for key, value in sortbyvalue.items():
        print(key,value, file=outfile)
    return sortbyvalue


def most_less_common_codon(codon_dict_sort):
    """The function select the most common and less codons and the number of them"""
    most_common_element= list(codon_dict_sort.keys())[len(codon_dict_sort)-1]
    number_most_common_element=list(codon_dict_sort.values())[len(codon_dict_sort)-1]
    less_common_element=list(codon_dict_sort.keys())[0]
    number_less_common_element=list(codon_dict_sort.values())[0]
    return (most_common_element,number_most_common_element,less_common_element,number_less_common_element)


print("Statistic")
def main():
    seq=get_sequences()
    length(seq)
    total_GC_content(seq)
    present_gc_content(seq)
    codon_dict = create_codon_dict(seq)
    codon_dict_sort = sort_by_value(codon_dict)

    #Below there is some loops to create the codon frequency.
    outfile=open (sys.argv[3],'w')
    aa_dict={}
    for codon in codon_dict_sort:
        aa = gencode[codon]
        if not aa in aa_dict:
             aa_dict[aa] = [codon]
        else:
            aa_dict[aa].append(codon)

    #The number for each codon is count.
    codon_count={}
    for codon in codon_dict_sort:
        aa = gencode[codon]
        if not aa in codon_count:
            codon_count[aa] = [codon_dict_sort[codon]]
        else:
            codon_count[aa].append(codon_dict_sort[codon])

    #The frequency for each codon
    for aa in codon_count:
        sum_of_aa = sum(codon_count[aa])
        for codon in aa_dict[aa]:           #aa_dict[aa] is a list of codons.
            print(aa, codon, round(codon_dict_sort[codon]/sum_of_aa, 3), file=outfile)                #print both the aminoacid and the codon tha represent the aa.
main()
