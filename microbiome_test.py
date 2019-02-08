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
            if my_id == 'RS_GCF_000158815.1~NZ_GG657738.1-#2 d__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Corynebacteriales;f__Micromonosporaceae;g__Micromonospora;s__ 272 6798591':
                sequences_16s[my_id] = ""
            if my_id == 'C1_0':
                sequences_16s[my_id] = ""
            if random.random() < fraction:
                sequences_16s[my_id] = ""
              
        else:
            if my_id in sequences_16s:
                sequences_16s[my_id] += line[:-1]
    
       
    return sequences_16s



def ConvertLibaryToKmerSets(library, K=2):
    import hashlib
    # use the hashlib function to generate unique values to store for each k-mer
    
    m = hashlib.md5()
    #m.update('TEST'.encode('utf-8'))
    #m.digest() -> returns the hashed value of the updated input
    new_lib = {}
    c = 0
    for k in library.keys():
        new_lib[k] = []
        # add your code here to build the k-mer set

        i_step = 0
        for kmer_i in range (len(library[k])//K):
          m.update((library[k][i_step:i_step + K]).encode('utf-8'))
          new_lib[k].append(m.digest())
          i_step += 1
        
    return new_lib


def JaccardIndex(s1, s2):
    print("s1, s2", len(s1), len(s2))
    s1_set = set()
    s2_set = set()

    if (len(s1) == len(s2)):
      s1_set = set(s1)
      s2_set = set(s2)

      numerator = float(len(s1_set.intersection(s2_set)))
      denominator = float(len(s1_set.union(s2_set)))

      return numerator/denominator
    
    best_JI = -1

    if (len(s1) > len(s2)):
      start_offset = 0
      while(start_offset + len(s2) < len(s1)):
        this_JI = JaccardIndex(s2, s1[start_offset:start_offset + len(s2)])
        start_offset += 1
        if (this_JI > best_JI): best_JI = this_JI

    elif (len(s1) < len(s2)):
      start_offset = 0
      while(start_offset + len(s1) < len(s2)):
        this_JI = JaccardIndex(s1, s2[start_offset:start_offset + len(s1)])
        start_offset += 1
        if (this_JI > best_JI): best_JI = this_JI

    return best_JI

def KmerMatch(sequence_kmer_set, library_kmer_set):
    best_score = 0.0
    best_match = None
    
    #add your code here to find the best kmer match
    
    return best_score, best_match


def AlignmentMatch(sequence, library):
    best_score = -10000000000
    best_match = None
    
    #add your code here to find the best match using alignment
    for item in library:
      (score, optloc, A) = alignment.local_align(item, sequence)
      if score > best_score:
        best_score = score
        best_match = item
    
    return best_score, best_match

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
   

   location_dict = {}
   location_dict['C1'] = ("Control", 1)
   location_dict['C2'] = ("Control", 2)
   location_dict['C3'] = ("Control", 3)
   
   location_dict['S8'] = ("Point-Mon", 1)
   location_dict['S9'] = ("Point-Mon", 2)
   location_dict['S11'] = ("Point-Mon", 3)
   
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
   sequences_16s = Load16SFastA(fn, fraction = 0.05)
   
   fn = "Fall2018CleanReads.fa"
   sample_sequences = Load16SFastA(fn, fraction = 0.05)


   
   print ("Loaded %d 16s sequences." % len(sequences_16s))
   print ("Loaded %d sample sequences." % len(sample_sequences))
   
   kmer_16s_sequences = ConvertLibaryToKmerSets(sequences_16s, K=6)
   kmer_sample_sequences = ConvertLibaryToKmerSets(sample_sequences, K=6)
#    k1 = 'RS_GCF_000158815.1~NZ_GG657738.1-#2 d__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Corynebacteriales;f__Micromonosporaceae;g__Micromonospora;s__ 272 6798591'
#    k2 = 'C1_0'
#    print ()

   for a16skey in kmer_16s_sequences.keys():
       for asamplekey in kmer_sample_sequences.keys():
           print(JaccardIndex(kmer_sample_sequences[asamplekey],kmer_16s_sequences[a16skey]))

