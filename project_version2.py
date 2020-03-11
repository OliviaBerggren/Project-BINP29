#!/usr/bin/env python3

#Fixa: Filen som man sätter in  måste sluta på fna!!!!!!!

import sys
usage='This program count the GC content i GC1, GC2 and GC3 for a fasta file, the length of the sequence and the codon usage.'


if len(sys.argv) != 2:
   print('This is a friendly reminder that there needs to be two arguments. The second argument needs to be a fasta file', file=sys.stderr)
   sys.exit()

if not '.fasta' in sys.argv[1] or '.fna' in sys.argv[1]:
    print('This is a friendly reminder that the input file needs to be a fasta file and have the substring .fasta or .fna', file=sys.stderr)
    sys.exit()


#The function collect all sequences and store them in a list.
def get_sequences():
    string=''
    seq=[]
    with open(args.infile.name, 'r') as fin:
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

#This function calculate the GC content in 1:st, 2:nd and 3:rd position of the codon.
def present_gc_content(seq):
    position_1=''
    position_2=''
    position_3=''
    for sequence in seq:
        position_1 +=sequence[::3]
        position_2 +=sequence[1::3]
        position_3 +=sequence[2::3]
    #print('GC content for 1:st position in codon is {}'.format(calculate_GC(position_1)))
    #print('GC content for 2:nd position in codon is {}'.format(calculate_GC(position_2)))
    #print('GC content for 3:rd position in codon is {}'.format(calculate_GC(position_3)))
    print('GC1 \t {}'.format(round(calculate_GC(position_1), 3)))
    print('GC2 \t {}'.format(round(calculate_GC(position_2), 3)))
    print('GC3 \t {}'.format(round(calculate_GC(position_3), 3)))


def calculate_GC(position):
    g = position.count('G')
    c = position.count('C')
    gc = (g+c)/len(position)
    return gc

def total_GC_content(seq):
    all_seq=''
    for sequence in seq:
        all_seq += sequence
    all_GC_content = calculate_GC(all_seq)
    print('GC_Percent \t {}'.format((round(all_GC_content*100, 3))))


#This function calculate the length of the sequence.
def length(seq):
    length=0
    number_of_seq=0
    for sequence in seq:
        number_of_seq +=1
        length +=len(sequence)
    length_seq = length/number_of_seq
    print('Length \t {}'.format(length_seq))

#The function creates a dictionary with key = codon and value = the number of that codon.
def create_codon_dict(seq):
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

#The function sort the dictionary by value
def sort_by_value(codon_dict):
    sortbyvalue=dict(sorted(codon_dict.items(), key=lambda t: t[1]))
    for key, value in sortbyvalue.items():
        print (key, value)
    return sortbyvalue


#The function select the most common and less codons and the number of them.
def most_less_common_codon(codon_dict_sort):
    most_common_element= list(codon_dict_sort.keys())[len(codon_dict_sort)-1]
    number_most_common_element=list(codon_dict_sort.values())[len(codon_dict_sort)-1]
    less_common_element=list(codon_dict_sort.keys())[0]
    number_less_common_element=list(codon_dict_sort.values())[0]
    return (most_common_element,number_most_common_element,less_common_element,number_less_common_element)

def present_codons(codon_dict_sort):
    most_common, number_most_common, less_common, number_less_common = most_less_common_codon(codon_dict_sort)
    print('Most common codon:\t{}\t{}'.format(most_common, number_most_common))
    print('Less common codon:\t{}\t{}'.format(less_common,number_less_common))

print("Statistic"+"\t"+"Species")
def main():
    seq=get_sequences()
    length(seq)
    total_GC_content(seq)
    present_gc_content(seq)
    codon_dict = create_codon_dict(seq)
    codon_dict_sort = sort_by_value(codon_dict)
    present_codons(codon_dict_sort)

main()
