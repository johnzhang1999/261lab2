# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 17:32:48 2018

@author: jdkan
"""
import numpy
import time
import alignment
import copy
import random
import matplotlib
# Your task is to *accurately* predict the primer melting points using machine 
# learning based on the sequence of the primer.

# Load the primers and their melting points.
def GetPhylum(seq_id):
    # input: 16s sequence ID - returns phylum
    # Kingdom, Phylum, Class, Order, Family, Genus, Species   
    p = k.split(';')[1][3:]
    return p




def Load16SFastA(path, fraction = 1.0):
    # from a file, read in all sequences and store them in a dictionary
    # sequences can be randomly ignored for testing by adjusting the fraction
    random.seed(11)
    
    infile = open(path, 'r')
    sequences_16s = {}
    c = 0
    my_seq = ""
    for line in infile:
        if ">" in line:
            my_id = line[1:-1]
            if random.random() < fraction:
                sequences_16s[my_id] = ""
            
            
        else:
            if my_id in sequences_16s:
                sequences_16s[my_id] += line[:-1]
    
       
    return sequences_16s



def ConvertLibaryToKmerSets(library, K=2):
    import hashlib
    # use the hashlib function to generate unique values to store for each k-mer
    
    #m = hashlib.md5()
    #m.update('TEST'.encode('utf-8'))
    #m.digest() -> returns the hashed value of the updated input
    new_lib = {}
    c = 0
    for k in library.keys():
        new_lib[k] = set()
        
        # add your code here to build the k-mer set
        
    return new_lib

def JaccardIndex(s1, s2):
    numerator = float(len(s1.intersection(s2)))
    denominator = float(len(s1.union(s2)))
    return numerator/denominator

def KmerMatch(sequence_kmer_set, library_kmer_set):
    best_score = 0.0
    best_match = None
    
    #add your code here to find the best kmer match
    
    return best_score, best_match


def AlignmentMatch(sequence, library):
    best_score = -10000000000
    best_match = None
    
    #add your code here to find the best match using alignment
    
    return best_score, best_match
if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
   

   location_dict = {}
   location_dict['C1'] = ("Control", 1)
   location_dict['C2'] = ("Control", 2)
   location_dict['C3'] = ("Control", 3)
   
   location_dict['S8'] = ("Point-Mon", 1)
   location_dict['S9'] = ("Point-Mon", 2)
   location_dict['S10'] = ("Point-Mon", 3)
   
   location_dict['S2'] = ("Point-Allegheny", 1)
   location_dict['S3'] = ("Point-Allegheny", 2)
   location_dict['S4'] = ("Point-Allegheny", 3)
   
   location_dict['S1'] = ("Sharpsburg", 1)
   location_dict['S15'] = ("Sharpsburg", 2)
   location_dict['S16'] = ("Sharpsburg", 3)
   location_dict['S17'] = ("Sharpsburg", 4)
   
   location_dict['S12'] = ("Braddock", 1)
   location_dict['S13'] = ("Braddock", 2)
   location_dict['S14'] = ("Braddock", 3)
   
   location_dict['S5'] = ("Neville Island", 1)
   location_dict['S6'] = ("Neville Island", 2)
   location_dict['S7'] = ("Neville Island", 3)
   
   fn = "bacterial_16s_genes.fa"
   sequences_16s = Load16SFastA(fn, fraction = 1.0)
   fn = "Fall2018CleanReads.fa"
   sample_sequences = Load16SFastA(fn, fraction = 1.0)
   
   print ("Loaded %d 16s sequences." % len(sequences_16s))
   print ("Loaded %d sample sequences." % len(sample_sequences))
   
   kmer_16s_sequences = ConvertLibaryToKmerSets(sequences_16s, K=6)
   kmer_sample_sequences = ConvertLibaryToKmerSets(sample_sequences, K=6)
   
