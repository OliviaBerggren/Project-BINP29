Project: BINP29
Author: Olivia Berggren
Date: 9 march 2020

>>>>>>>>>
Aim:
The aim was to create a software that could give information about a sequence in a fasta file. The software should be easy to use and the output
should be easy to understand.


>>>>>>> Python script usage
The script calculate the overall GC-content, the GC-content codon position 1,2 and 3, codon usage, codon frequency
and the length of the sequence.

The script needs to be run with a fasta file and have two outputs files. There should be four arguments:
1. name of the scrips (software_fasta1.py)
2. fasta file
3. name for the output-file for the codon usage (the user decide the name)
4. name for the output-file for the codon frequency (the user decide the name)
The information about GC will be print to standard output.

The fasta file must have the suffix .fna.

The software can be tested with the data files, BRAC2_human.fna, BRAC2_mouse.fna, BRAC2_elephant.fna.
Download them from Github resository:
https://github.com/OliviaBerggren/Project-BINP29

Example:
python3 software_fasta1.py BRAC2_human.fna codon_usage_human_BRAC2 codon_frequency_human_BRAC2
