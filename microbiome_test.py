# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 17:32:48 2018

@author: jdkan

Group Members: Anupam Pokharel, John Zhang
"""
import numpy
import time
import alignment
import copy
import random
import pickle
import matplotlib.pyplot as plt
from ast import literal_eval as make_tuple
from collections import Counter

# https://codeinterview.io/NLNZMXPANV

# Your task is to *accurately* predict the primer melting points using machine
# learning based on the sequence of the primer.

# Load the primers and their melting points.


def GetPhylum(seq_id):
    # input: 16s sequence ID - returns phylum
    # Kingdom, Phylum, Class, Order, Family, Genus, Species
    p = seq_id.split(';')[1][3:]
    return p


def Load16SFastA(path, fraction=1.0):
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

    m = hashlib.md5()
    # m.update('TEST'.encode('utf-8'))
    # m.digest() -> returns the hashed value of the updated input
    new_lib = {}
    c = 0
    for k in library.keys():
        new_lib[k] = []
        # add your code here to build the k-mer set

        i_step = 0
        for kmer_i in range(len(library[k])-K+1):
            # m.update((library[k][i_step:i_step + K]).encode('utf-8'))
            # new_lib[k].append(m.digest())
            new_lib[k].append(library[k][i_step:i_step + K])
            i_step += 1

    return new_lib


def JaccardIndex(s1, s2):
    # print("s1, s2", len(s1), len(s2))
    s1 = set(s1)
    s2 = set(s2)

    ji_balance = len(s1)/len(s2)

    numerator = float(len(s1.intersection(s2)))
    denominator = float(len(s1.union(s2)))

    return ji_balance * numerator/denominator

def KmerMatchWrapper(all_sequences, kmer_lib, keys=None):
    result = {}
    for key, seq in all_sequences.items(): 
        if key in keys:
            score, match_key = KmerMatch(seq, kmer_lib)
            result[key] = (score, match_key)
    return result

def KmerMatch(sequence_kmer_set, library_kmer_set):
    best_score = 0.0
    best_match = None
    for key,val in library_kmer_set.items():
        this_score = JaccardIndex(val, sequence_kmer_set)

        if (this_score > best_score):
            best_score = this_score
            best_match = key

    return best_score, best_match

def AlignmentMatch(sequence, library):
    best_score = -10000000000
    best_match = None
    counter = 0

    # add your code here to find the best match using alignment
    for key in library.keys():
        # print(library[key])
        (score, optloc, A) = alignment.local_align(library[key], sequence)
        print('checking library item:', counter)
        if score > best_score:
            best_score = score
            best_match = key
        counter += 1

    return best_score, best_match

def generateTruth(lib, sample):
    align_score_file = open("./align_scores.txt", "w")
    counter = 0
    itemList = sample.items()
    sampleLen = len(itemList)
    for key, value in itemList:
        # if counter % 100 == 0:
        print('aligning samples:', counter)
        this_align_score = AlignmentMatch(value, lib)
        align_score_file.write(key + '\t' + str(this_align_score) + '\n')
        counter += 1
    align_score_file.close()

def processTruth(sample):
	f = open('./truth_percents.txt', 'w')
	with open('./align_scores.txt', 'r') as file:
		for line in file:
			s = line.strip().split('\t')
			samp = s[0]
			score = make_tuple(s[1])[0]
			target_score = len(sample[samp]) * 10
			percent = score/target_score
			f.write(samp + '\t' + str(percent) + '\t' + make_tuple(s[1])[1] + '\n')
	f.close()

def convertTruthsToDict (filepath): #alignment scores dictionary
	result = {}
	with open(filepath, 'r') as file:
		for line in file:
			s = line.strip().split('\t')
			key = s[0]
			val = s[1]
			result[key] = val

	return result

def processTruthAbove90():
	f = open('./truth_percents_good.txt', 'w')
	with open('./truth_percents.txt', 'r') as file:
		for line in file:
			s = line.strip().split('\t')
			score = s[1]
			if (float(score) >= .9):
				f.write(str(score) + '\n')

def ExtractKeys(filepath):
	result = []
	with open(filepath, 'r') as file:
		for line in file:
			s = line.strip().split('\t')
			key = s[0]
			result.append(key)

	return result

def assessment(align_dict, kmer_dict, kmer_thresh): #if the align score and kmer score correspond
    result = {}
    correct = 0
    total = len(align_dict.keys())
    assert (total == len(kmer_dict.keys()))
    for key in align_dict.keys():
        align_val = float(align_dict[key])
        kmer_val = float(kmer_dict[key][0])
        cond = ((align_val > .9 and kmer_val > kmer_thresh) 
                            or (align_val < .9 and kmer_val < kmer_thresh))
        if (cond): 
            result[key] = True
            correct += 1
        else: result[key] = False

    return result, (correct/total)

def runAssessment (sequences_16s, sample_sequences, K):
	print("\n\n\n", K)
	print("Loaded %d 16s sequences." % len(sequences_16s))
	print("Loaded %d sample sequences." % len(sample_sequences))

	kmer_16s_sequences = ConvertLibaryToKmerSets(sequences_16s, K)
	kmer_sample_sequences = ConvertLibaryToKmerSets(sample_sequences, K)

    # uncomment the below two lines to generate truth (WARNING: runs extremely slow)
	# generateTruth(sequences_16s, sample_sequences)
	# processTruth(sample_sequences)
	filepath = './truth_percents.txt'
	keys = ExtractKeys(filepath)

	align_socre_dict = convertTruthsToDict(filepath)
	kmer_scores_dict = KmerMatchWrapper(kmer_sample_sequences, kmer_16s_sequences, keys)

	# print(align_socre_dict)
	# print('###############')
	# print(kmer_scores_dict)
	
	for i in range (0, 20):
		thresh = i/20
		final_matches, success_ratio = assessment(align_socre_dict, kmer_scores_dict, thresh)
		print (thresh, success_ratio)

def calcMixFrac(lib, samples, K, thresh):
    # lib: dictionary of k-mer sets
    # samples: dictionary of k-mer sets from all places
    # returns a dict (keys: C1, ...) of dict (keys: phylums)

    print("Loaded %d 16s sequences." % len(lib))
    print("Loaded %d sample sequences." % len(samples))

    kmer_lib = ConvertLibaryToKmerSets(lib, K)
    kmer_samples = ConvertLibaryToKmerSets(samples, K)
    l = len(kmer_samples)

    unmatched = {}
    result = {}
    counter = 0


    for loc_i, seq in kmer_samples.items():
        if (counter % 10 == 0): print('progress:', counter/l)
        loc = loc_i.split('_')[0]
        if not loc in result: 
            result[loc] = {}
        score, match = KmerMatch(seq, kmer_lib)
        phylum = GetPhylum(match)

        if score > thresh:
            if not phylum in result[loc]:
                result[loc][phylum] = 0
            result[loc][phylum] += 1
        else:
            unmatched[loc_i] = seq
        counter += 1

    return result, unmatched
        
def mergeLoc(raw, loc_dict):
    result = {}
    for loc_i, phylum_count_dict in raw.items():
        loc = loc_dict[loc_i][0]
        if not loc in result:
            result[loc] = Counter({})
        result[loc] += Counter(phylum_count_dict)
    return result

def normalize(raw):
    arr = numpy.array(raw)
    sum = numpy.sum(arr)
    return numpy.round(arr/sum*100)

def write(path, obj):
    pickle_out = open(path,'wb')
    pickle.dump(obj, pickle_out)
    pickle_out.close()
    print('written:', path)

def load(path):
    pickle_in = open(path,'rb')
    obj = pickle.load(pickle_in)
    return obj

def generatePlot(loc, phylum_count_dict, path):
    labels = list(phylum_count_dict.keys())
    counts = list(phylum_count_dict.values())
    sizes = normalize(counts)

    plt.title(loc)
    plt.pie(sizes, labels=labels, autopct='%1.1f%%',
        startangle=90)
    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.savefig(path)
    plt.close()


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
    sequences_16s = Load16SFastA(fn, fraction=0.05)

    fn = "Fall2018CleanReads.fa"
    sample_sequences = Load16SFastA(fn, fraction=0.01)

    # uncomment the below two lines to run K and K-mer threshold assessments.
    # for K in range (2, 11, 1):
    # 	runAssessment(sequences_16s, sample_sequences, K)

    result, unmatched = calcMixFrac(sequences_16s, sample_sequences, 8, 0.6)

    write('matched_result.pickle', result)
    write('unmatched.pickle', unmatched)

    merged = mergeLoc(result, location_dict)

    write('matched_merged.pickle', merged)

    for loc, phylum_count_dict in merged.items():
        path = loc + '.png'
        generatePlot(loc, phylum_count_dict, path)