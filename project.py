#!/usr/bin/env python3

#Fixa: Filen som man sätter in  måste sluta på fna!!!!!!!

import sys
import argparse

usage='This program count the GC content i GC1, GC2 and GC3 for a fasta file, the longest/shortest sequence, the average lenght of the sequences and how many times each codon occur.'
parser=argparse.ArgumentParser(description=usage)
parser.add_argument('-v', '--version',action='version',
    version='%(prog)s 1.0')
parser.add_argument('infile', type=argparse.FileType('r'))
args=parser.parse_args()


if len(sys.argv) != 2:
   print('There needs to be two arguments. The second argument needs to be a fasta file', file=sys.stderr)
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
    print('GC content for 1:st position in codon is {}'.format(calculate_GC(position_1)))
    print('GC content for 2:nd position in codon is {}'.format(calculate_GC(position_2)))
    print('GC content for 3:rd position in codon is {}'.format(calculate_GC(position_3)))


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
    print('The total_GC content is {}'.format((all_GC_content)))


#This function calculate the average length of the sequences.
# def average_length(seq):
def average_length(seq):
    length=0
    number_of_seq=0
    for sequence in seq:
        number_of_seq +=1
        length +=len(sequence)
    average = length/number_of_seq
    print('The average length for a sequence is {} nucleotides'.format(average))

#The two functions below calculate the longest and shortest sequence.
def maximum_length(seq):
    return (len(max(seq, key=len)), len(min(seq, key=len)))

def present_maxium_length(seq):
    max_length, min_length = maximum_length(seq)
    print('The longest sequence have {} nucleotides'.format(max_length))
    print('The shortest sequence have {} nucleotides'.format(min_length))


#This function create the reverse compliment sequences.!!!!!!!!!!!!!!!!!!!!!!!!!!!!
def reverse_complement(seq):
    reverse_complement_seq=[]
    for sequence in seq:
        sequence=sequence[::-1]
        reverse_complement_seq.append(sequence.translate(str.maketrans('ACGT','TGCA')))
    #print("The reverse compliment sequence are '{}'".format(reverse_complement_seq[-1]))


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
    print (codon_dict)
    return codon_dict

#The function sort the dictionary by value
def sort_by_value(codon_dict):
    sortbyvalue=dict(sorted(codon_dict.items(), key=lambda t: t[1]))
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
    print('The most common codon is {} and it occur {} times'.format(most_common, number_most_common))
    print('The less common codon is {} and it occur {} times'.format(less_common,number_less_common))


def main():
    seq=get_sequences()
    present_gc_content(seq)
    total_GC_content(seq)
    average_length(seq)
    present_maxium_length(seq)
    reverse_complement(seq)
    codon_dict = create_codon_dict(seq)
    codon_dict_sort = sort_by_value(codon_dict)
    print(codon_dict_sort)
    present_codons(codon_dict_sort)

main()
