"""
This script will compute the Burrows–Wheeler transform in a memmory efficient way
and uses it to map the reads from a fasta formated file on a reference genome. 
This script map the reads in sens and antisens.
This script support multiple substitions, insertions and deletions, using argparse options.

Script usage:
-------------
python Burrow.py [-h] -fasta FASTA -reads READS [-sub SUBSTITUTION]
python Burrow.py -fasta Hu-1.fasta -reads READSsars_cov_2_1e6.fasta -sub 3 -ins 3 -del 3
"""

__author__ = "Vander Meersche Yann"
__copyright__ = "Universite de Paris"
__credits__ = ["Yann Vander Meersche"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Vander Meersche Yann"
__email__ = "yann-vm@hotmail.fr"
__status__ = "Developpement"

#Modules =======================================================================
import os
import sys
import argparse

import numpy as np
import matplotlib.pyplot as plt
#===============================================================================


def args():
    #Declaration of expexted arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-fasta", "--fasta", help="Input fasta sequences file.", type=str, required=True)
    parser.add_argument("-reads", "--reads", help="Input reads file.", type=str, required=True)
    parser.add_argument("-sub", "--substitution", help="Number of accepted substitution (default = 0).", type=int, required=False, default=0)
    parser.add_argument("-ins", "--insertion", help="Number of accepted insertion (default = 0).", type=int, required=False, default=0)
    parser.add_argument("-del", "--deletion", help="Number of accepted deletion (default = 0).", type=int, required=False, default=0)
    args = parser.parse_args()

    return args.fasta, args.reads, args.substitution, args.insertion, args.deletion


def parse_fasta(fasta_file):
    """
    Parse a fasta formated file into a dictionnary {genome_id: sequence}.
    """
    seq_dict = {}

    with open(fasta_file, "r") as f_in:
        #seq_id = ""
        for line in f_in:
            if line.startswith(">"):
                seq_id = line[1:].strip()
                seq_dict[seq_id] = ""
            else:
                seq_dict[seq_id] += line.strip()

    return seq_dict


def parse_reads(reads_file):
    """
    Parse a fasta formated file into a dictionnary {read_id: sequence}.
    """
    reads_dict = {}

    with open(reads_file, "r") as f_in:
        for line in f_in:
            if line.startswith(">"):
                read_id = line[1:].strip()
                read_id = read_id.split(".")[0].split("-")[0].split("+")[0]
                read_id = read_id[:-1]
            else:
                reads_dict[read_id] = line.strip()

    return reads_dict


'''
def compute_bwt(seq):
    """
    Basic BWT search
    """

    seq = "$" + seq
    list_perm = []

    for i in range(len(seq)):
        list_perm.append(seq)
        seq = seq[-1] + seq[:-1]

    bwt = ""
    for seq in sorted(list_perm):
        bwt += seq[-1]

    return(bwt)
'''

# Optimised BWT ################################################################
def sa(seq):
    """
    Create the suffix array (seq, rank) of the sequence.
    As only half of the matrix is stored and sorted, this implementation is
    memory efficient and faster that the previous implementation.
    Works with every sequence (not just ATGC alphabet).
    """
    half_mat = []
    for i in range(len(seq)):
        half_mat.append((seq[i:], i))
    half_mat = sorted(half_mat)

    seq_index = []
    for tpl in half_mat:
        seq_index.append(tpl[1])

    return seq_index


def compute_bwt(seq):
    """
    Compute the Burrows-Wheeler transform with a suffix array.
    Works with every sequence (not just ATGC alphabet).
    """
    seq = "$" + seq

    bwt = []
    for seq_index in sa(seq):
        bwt.append(seq[seq_index-1])

    return ''.join(bwt)
################################################################################

"""
def generate_C(bwt):
    alphabet="$ACGT"
    c = {}
    for n in alphabet:
        c[n] = 0
    for n in bwt:
        c[n] += 1

    s = 0
    for n in alphabet:
        tmp = c[n]
        c[n] = s
        s += tmp

    return(c)

def generate_num(bwt): 
    num = []
    c = {}

    for base in bwt:
        if base not in c:
            c[base] = 0
        num.append(c[base])
        c[base] += 1

    return(num)
"""


def generate_num_C(bwt):
    """
    Compute the num and C in one step. Slightly more efficient that the generate_C
    and generate_num functions.
    Works with every sequence (not just ATGC alphabet).
    """
    num = []
    c = {}

    for base in bwt:
        if base not in c:
            c[base] = 0
        num.append(c[base])
        c[base] += 1

    s = 0
    for n in sorted(c):
        tmp = c[n]
        c[n] = s
        s += tmp

    return(num, c)


def gener_seq(bwt, c, num):
    """
    Generate seq
    """
    seq = ""
    i = 0

    while bwt[i] != "$":
        seq = bwt[i] + seq
        i = c[bwt[i]] + num[i]

    return(seq)


def gener_pos(bwt, c, num):
    """
    Generate pos
    """
    pos = [-1] * len(bwt)
    seq = ""
    x = len(bwt) - 2
    p = 0

    while bwt[p] != "$":
        pos[p] = x
        p = c[bwt[p]] + num[p]
        x -= 1

    return(pos)


def create_occ(bwt):
    """
    Create occ
    """
    dict_occ = {}
    occ = []

    for n in bwt:
        if n not in dict_occ:
            dict_occ[n] = 0

    for i in range(len(bwt)):
        occ.append({})
        dict_occ[bwt[i]] += 1
        for n in dict_occ:
            occ[i][n] = dict_occ[n]

    return(occ)


'''
def pos_read(read, c, occ):
    """
    Basic pos_read
    """
    i = len(read) - 1
    b = c[read[i]]
    e = c[read[i]] + occ[-1][read[i]] - 1
    
    while i >= 1 and b <= e:
        i = i - 1
        b = c[read[i]] + occ[b-1][read[i]]
        e = c[read[i]] + occ[e][read[i]] - 1
    return(b, e)
'''

# Using substitution, insertion et délétions ###################################
def pos_read_sub(read, c, occ, user_nb_sub):
    """
    Recherche le read en prenant en compte de 0 (pos_read classique) à n substitutions
    """
    nb_subs = 0

    i = len(read) - 1
    b = c[read[i]]
    e = c[read[i]] + occ[-1][read[i]] - 1
    
    while i >= 1 and b <= e:
        b_save = b
        e_save = e

        i = i - 1
        b = c[read[i]] + occ[b-1][read[i]]
        e = c[read[i]] + occ[e][read[i]] - 1
        
        #S'il y a une "erreur" dans le read (b > e), on test les autres bases.
        if b > e and nb_subs < user_nb_sub:
            for base in ["A", "T", "G", "C"]:
                if base != read[i]:
                    b = c[base] + occ[b_save-1][base]
                    e = c[base] + occ[e_save][base] - 1
                    if b <= e:
                        nb_subs += 1
                        break
    return(b, e)


def pos_read_ins(read, c, occ, user_nb_ins):
    """
    Recherche le read en prenant en compte de 0 (pos_read classique) à n insertions
    """
    nb_ins = 0

    i = len(read) - 1
    b = c[read[i]]
    e = c[read[i]] + occ[-1][read[i]] - 1

    while i >= 1 and b <= e:
        b_save = b
        e_save = e

        i = i - 1
        b = c[read[i]] + occ[b-1][read[i]]
        e = c[read[i]] + occ[e][read[i]] - 1
        
        #Du temps qu'il y a des "erreurs" dans le read (b > e), on saute une position.
        while b > e and nb_ins < user_nb_ins:
            i = i - 1
            b = c[read[i]] + occ[b_save-1][read[i]]
            e = c[read[i]] + occ[e_save][read[i]] - 1
            nb_ins += 1

    return(b, e)


def pos_read_del(read, c, occ, user_nb_del):
    """
    Recherche le read en prenant en compte de 0 (pos_read classique) à n délétion
    """
    nb_del = 0

    i = len(read) - 1
    b = c[read[i]]
    e = c[read[i]] + occ[-1][read[i]] - 1

    while i >= 1 and b <= e:
        b_save = b
        e_save = e

        i = i - 1
        b = c[read[i]] + occ[b-1][read[i]]
        e = c[read[i]] + occ[e][read[i]] - 1

        #S'il y a une "erreur" dans le read (b > e), on test les autres bases.
        if b > e and nb_del < user_nb_del:
            for base in ["A", "T", "G", "C"]:
                if base != read[i]:
                    b = c[base] + occ[b_save-1][base]
                    e = c[base] + occ[e_save][base] - 1
                    if b <= e:
                        nb_del += 1    #Mais on ne change pas i pour la prochaine iteration
                        i = i + 1
                        break
    return(b, e)
################################################################################


def map_read(b, e, pos):
    """
    Map the read on the reference genome
    """
    if b == e:
        return(pos[b] + 1)    #Position on the reference genome
    elif e - b > 1:
        return(-2)    #Multiple position mapping reads (not supported)
    else:
        return(-1)    #Not mapped reads


def rev_complement(seq):
    """
    Compute the reverse complement of a DNA sequence
    """
    complement = {"A" : "T", "T" : "A", "G" : "C", "C" : "G"}
    seq = list(seq[::-1])

    for i, base in enumerate(seq):
        seq[i] = complement[seq[i]]

    return "".join(seq)


def stats(map_pos_list):
    """
    Calculate some informative statistics
    """
    nb_map = sum(map_pos_list > -1)
    pourc_map = nb_map / len(map_pos_list) * 100
    print(f"Number of mapped reads: {nb_map} ({pourc_map:.3f}%)")

    nb_unmap = sum(map_pos_list == -1)
    pourc_unmap = nb_unmap / len(map_pos_list) * 100
    print(f"Number of unmapped reads: {nb_unmap} ({pourc_unmap:.3f})%")

    nb_mult = sum(map_pos_list == -2)
    pourc_mult = nb_mult / len(map_pos_list) * 100
    print(f"Number of multiple positions reads: {nb_mult} ({pourc_mult:.3f})%")


def main():
    fasta_file, reads_file, user_nb_sub, user_nb_ins, user_nb_del = args()

    fasta_dict = parse_fasta(fasta_file)
    reads_dict = parse_reads(reads_file)

    for key in fasta_dict:
        seq = fasta_dict[key]

        #Normal
        bwt = compute_bwt(seq)
        num, C = generate_num_C(bwt)
        seqn = gener_seq(bwt, C, num)
        posn = gener_pos(bwt, C, num)
        occ = create_occ(bwt)
        
        #Reverse
        rev_seq = rev_complement(seq)
        rev_bwt = compute_bwt(rev_seq)
        rev_num, rev_C = generate_num_C(rev_bwt)
        rev_seqn = gener_seq(rev_bwt, rev_C, rev_num)
        rev_posn = gener_pos(rev_bwt, rev_C, rev_num)
        rev_occ = create_occ(rev_bwt)
        
        map_pos_list = []

        #Map the reads in normal and reverse sens.
        #Look for substitution, then insertions, and finally 
        for read_key in reads_dict:

            #Substitution Sens
            b, e = pos_read_sub(reads_dict[read_key], C, occ, user_nb_sub)
            map_pos = map_read(b, e, posn)

            #Substitution Reverse
            if map_pos == -1 or map_pos == -2:
                b, e = pos_read_sub(reads_dict[read_key], rev_C, rev_occ, user_nb_sub) 
                map_pos = map_read(b, e, rev_posn)

            #Insertion Sens
            if user_nb_ins > 0:
                if map_pos == -1 or map_pos == -2:
                    b, e = pos_read_ins(reads_dict[read_key], C, occ, user_nb_ins)
                    map_pos = map_read(b, e, posn)
                    
                #Insertion Reverse
                if map_pos == -1 or map_pos == -2:
                    b, e = pos_read_ins(reads_dict[read_key], rev_C, rev_occ, user_nb_ins) 
                    map_pos = map_read(b, e, rev_posn)
                    
            #Deletion Sens
            if user_nb_del > 0:
                if map_pos == -1 or map_pos == -2:
                    b, e = pos_read_del(reads_dict[read_key], C, occ, user_nb_del)
                    map_pos = map_read(b, e, posn)

                #Deletion Reverse
                if map_pos == -1 or map_pos == -2:
                    b, e = pos_read_del(reads_dict[read_key], rev_C, rev_occ, user_nb_del)
                    map_pos = map_read(b, e, rev_posn)         
                        
            map_pos_list.append(map_pos)


        map_pos_list = np.array(map_pos_list)
        
        #Print statistics in the terminal
        stats(map_pos_list)


        #Display the histogram
        binwidth = 100
        mapped_pos = map_pos_list[map_pos_list > -1]
        plt.hist(mapped_pos, bins=np.arange(min(mapped_pos), max(mapped_pos) + binwidth, binwidth))
        plt.ylim(0, 4000)
        plt.savefig(f"plot_{user_nb_sub}_subs_{user_nb_ins}_ins_{user_nb_del}_del")



if __name__ == '__main__':
    main()
