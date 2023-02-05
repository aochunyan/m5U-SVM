#!/usr/bin/env python
# _*_coding:utf-8_*_
# @Time : 2022.10.20
# @Author : aochunyan
# @Email : acy196707@163.com
# @IDE : PyCharm
# @File : m5USVM.py

from collections import Counter
import itertools
import pickle
import sys
from feature_scripts import read_fasta_sequences
from Bio import SeqIO
import pandas as pd
import numpy as np
from gensim import utils
from gensim.models.word2vec import Word2Vec
######################################################################################################

myDiIndex = {
    'AA': 0, 'AC': 1, 'AG': 2, 'AT': 3,
    'CA': 4, 'CC': 5, 'CG': 6, 'CT': 7,
    'GA': 8, 'GC': 9, 'GG': 10, 'GT': 11,
    'TA': 12, 'TC': 13, 'TG': 14, 'TT': 15
}
myTriIndex = {
    'AAA': 0, 'AAC': 1, 'AAG': 2, 'AAT': 3,
    'ACA': 4, 'ACC': 5, 'ACG': 6, 'ACT': 7,
    'AGA': 8, 'AGC': 9, 'AGG': 10, 'AGT': 11,
    'ATA': 12, 'ATC': 13, 'ATG': 14, 'ATT': 15,
    'CAA': 16, 'CAC': 17, 'CAG': 18, 'CAT': 19,
    'CCA': 20, 'CCC': 21, 'CCG': 22, 'CCT': 23,
    'CGA': 24, 'CGC': 25, 'CGG': 26, 'CGT': 27,
    'CTA': 28, 'CTC': 29, 'CTG': 30, 'CTT': 31,
    'GAA': 32, 'GAC': 33, 'GAG': 34, 'GAT': 35,
    'GCA': 36, 'GCC': 37, 'GCG': 38, 'GCT': 39,
    'GGA': 40, 'GGC': 41, 'GGG': 42, 'GGT': 43,
    'GTA': 44, 'GTC': 45, 'GTG': 46, 'GTT': 47,
    'TAA': 48, 'TAC': 49, 'TAG': 50, 'TAT': 51,
    'TCA': 52, 'TCC': 53, 'TCG': 54, 'TCT': 55,
    'TGA': 56, 'TGC': 57, 'TGG': 58, 'TGT': 59,
    'TTA': 60, 'TTC': 61, 'TTG': 62, 'TTT': 63
}
baseSymbol = 'ACGT'

myDictDefault = {
    'PseDNC': {
               'RNA': ['Rise (RNA)', 'Roll (RNA)', 'Shift (RNA)', 'Slide (RNA)', 'Tilt (RNA)', 'Twist (RNA)']},

}
myDataFile = {

    'PseDNC': {'DNA': 'didnaPhyche.data', 'RNA': 'dirnaPhyche.data'},

}
AA = 'ACGT'
######################################################################################################
def get_min_sequence_length(seqs):
    minLen = 10000
    for i in seqs:
        if minLen > len(i[1]):
            minLen = len(i[1])
    return minLen

def get_cksnap(seqs, gap):
    encodings = []
    aaPairs = []
    for aa1 in AA:
        for aa2 in AA:
            aaPairs.append(aa1 + aa2)
    for i in seqs:
        sequence = i
        code = [ ]
        for g in range(gap + 1):
            myDict = {}
            for pair in aaPairs:
                myDict[pair] = 0
            sum = 0
            for index1 in range(len(sequence)):
                index2 = index1 + g + 1
                if index1 < len(sequence) and index2 < len(sequence) and sequence[index1] in AA and sequence[
                    index2] in AA:
                    myDict[sequence[index1] + sequence[index2]] = myDict[sequence[index1] + sequence[index2]] + 1
                    sum = sum + 1
            for pair in aaPairs:
                code.append(myDict[pair] / sum)
        encodings.append(code)
    return encodings

#####################################################################################################
def get_enac(seqs, window):
    encodings = []
    for sequence in seqs:
        code = []
        for j in range(len(sequence)):
            if j < len(sequence) and j + window <= len(sequence):
                count = Counter(sequence[j:j + window])
                for key in count:
                    count[key] = count[key] / len(sequence[j:j + window])
                for aa in AA:
                    code.append(count[aa])
        encodings.append(code)
    return encodings
###########################################################################################################
def kmerArray(sequence, k):
    kmer = []
    for i in range(len(sequence) - k + 1):
        kmer.append(sequence[i:i + k])
    return kmer

def get_kmer(seqs, k):
    encoding = []
    header = []
    NA = 'ACGT'

    for kmer in itertools.product(NA, repeat=k):
        header.append(''.join(kmer))
    for i in seqs:
        sequence = i
        kmers = kmerArray(sequence, k)
        count = Counter()
        count.update(kmers)
        for key in count:
            count[key] = count[key] / len(kmers)
        code = []
        for j in range(len(header)):
            if header[j] in count:
                code.append(count[header[j]])
            else:
                code.append(0)
        encoding.append(code)
    return encoding

#####################################################################################PAAC
def check_Pse_arguments(method,type,weight,lamadaValue):
    myNum = 0

    myIndex = []
    myProperty = {}
    dataFile = ''
    data_path = "./data/"

    if myNum ==0:
        myIndex = myDictDefault[method][type]
        dataFile = myDataFile[method][type]
    if dataFile != '':
        with open(data_path + dataFile, 'rb') as f:
            myProperty = pickle.load(f)

    if len(myIndex) == 0 or len(myProperty) == 0:
        print('Error: arguments is incorrect.')
        sys.exit(1)

    return myIndex, myProperty, lamadaValue, weight

def get_kmer_frequency(sequence, kmer):
    myFrequency = {}
    for pep in [''.join(i) for i in list(itertools.product(baseSymbol, repeat=kmer))]:
        myFrequency[pep] = 0
    for i in range(len(sequence) - kmer + 1):
        myFrequency[sequence[i: i + kmer]] = myFrequency[sequence[i: i + kmer]] + 1
    for key in myFrequency:
        myFrequency[key] = myFrequency[key] / (len(sequence) - kmer + 1)
    return myFrequency

def correlationFunction(pepA, pepB, myIndex, myPropertyName, myPropertyValue):
    CC = 0
    for p in myPropertyName:
        CC = CC + (float(myPropertyValue[p][myIndex[pepA]]) - float(myPropertyValue[p][myIndex[pepB]])) ** 2
    return CC / len(myPropertyName)

def get_theta_array(myIndex, myPropertyName, myPropertyValue, lamadaValue, sequence, kmer):
    thetaArray = []
    for tmpLamada in range(lamadaValue):
        theta = 0
        for i in range(len(sequence) - tmpLamada - kmer):
            theta = theta + correlationFunction(sequence[i:i + kmer],
                                                sequence[i + tmpLamada + 1: i + tmpLamada + 1 + kmer], myIndex,
                                                myPropertyName, myPropertyValue)
        thetaArray.append(theta / (len(sequence) - tmpLamada - kmer))
    return thetaArray



def get_psednc(seqs, myPropertyName, myPropertyValue, lamadaValue, weight):
    encodings = []
    myIndex = myDiIndex

    for i in seqs:
        sequence = i
        code = []
        dipeptideFrequency = get_kmer_frequency(sequence, 2)
        thetaArray = get_theta_array(myIndex, myPropertyName, myPropertyValue, lamadaValue, sequence, 2)
        for pair in sorted(myIndex.keys()):
            code.append(dipeptideFrequency[pair] / (1 + weight * sum(thetaArray)))
        for k in range(17, 16 + lamadaValue + 1):
            code.append((weight * thetaArray[k - 17]) / (1 + weight * sum(thetaArray)))
        encodings.append(code)
    return encodings

#####################################################################################PAAC

def kmer_seq(sequence, k):
    list_seq=list(sequence)
    return (list("".join(list_seq[i:i+k]) for i in range(len(sequence)-k+1)))

def get_w2v(sequences,model):
    word_features=np.zeros((len(sequences), 100))
    for i in range(len(sequences)):
        array=np.array([model.wv[w] for w in sequences[i] if w in model.wv ])
        idx = np.argwhere(np.all(array[..., :] == 0, axis=1))
        array = np.delete(array, idx, axis=0)
        S= pd.Series(array.mean(axis=0))
        word_features[i]=S.values
    return word_features

def w2v_kmer_corpus(seq_file):
    with open('./data/w2v_kmer.txt',"w") as fil1:
            for seq_record in SeqIO.parse(seq_file, "fasta"):
                seq_id=seq_record.id
                seq=seq_record.seq
                sub_list=kmer_seq(seq,7)
                fil1.write(str(seq_id)+",")
                fil1.writelines(str(x)+"," for x in sub_list)
                fil1.write("\n")
#####################################################################################PAAC

def get_features(file, type_=None):
    sequences, names = read_fasta_sequences.read_RNA_sequences(file)
    encoding1 = get_cksnap(sequences, 3)
    encoding2 = get_enac(sequences, 3)
    encoding3 = get_kmer(sequences, 4)
    my_property_name, my_property_value, lamada_value, weight = check_Pse_arguments( "PseDNC", "RNA", 0.1, 5)
    encoding4 = get_psednc(sequences, my_property_name, my_property_value, 5, 0.1)

    w2v_kmer_corpus(file)

    with utils.smart_open('./data/w2v_kmer.txt', 'r', encoding='utf-8-sig', ) as infile:
        its_list = list(infile)

    list_seq = []
    for x in range(len(its_list)):
        list_seq.append(its_list[x].split(",")[1:])

    if type_ == "Full_transcript":
        model = Word2Vec.load("./model/F_w2v.model")
    else:
        model = Word2Vec.load("./model/M_w2v.model")
    encoding5 = get_w2v(np.array(list_seq), model)


    encoding1 = pd.DataFrame(encoding1)
    encoding2 = pd.DataFrame(encoding2)
    encoding3 = pd.DataFrame(encoding3)
    encoding4 = pd.DataFrame(encoding4)
    encoding5 = pd.DataFrame(encoding5)

    encoding =pd.concat([encoding1,encoding2,encoding3,encoding4,encoding5],axis=1,ignore_index=True)
    return encoding




