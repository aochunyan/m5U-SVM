#!/usr/bin/env python
# _*_coding:utf-8_*_
# @Time : 2022.10.20
# @Author : aochunyan
# @Email : acy196707@163.com
# @IDE : PyCharm
# @File : m5USVM.py

import sys, os, platform
import pandas as pd
import numpy as np
import joblib
from feature_scripts import feature_generate
from feature_scripts import read_fasta_sequences
from feature_scripts import feature_selection
import argparse
import time
import warnings
warnings.filterwarnings('ignore')

if not os.path.exists("results"):
    os.makedirs("results")

pPath = os.path.split(os.path.realpath("__file__"))[0]
sys.path.append(pPath)
data_path = os.path.abspath(pPath) + r'\data' if platform.system() == 'Windows' else os.path.abspath(pPath)  + r'/data'

feature_index_F = pd.read_csv(r"./data/F_features_index.csv")
feature_index_M = pd.read_csv(r"./data/M_features_index.csv")

#######################################################################################

def  m5Usvm_pred(FastaDATA,Modes):

    print("Sequence checking......")
    seq, name = read_fasta_sequences.read_RNA_sequences(FastaDATA)

    if Modes=="Full_transcript":
        print("feature extraction......")
        Data = feature_generate.get_features(FastaDATA,type_=Modes)
        print("feature selection......")
        x = feature_selection.select_features(Data,feature_index_F)
        feature = pd.DataFrame(x)
        scaler = joblib.load("./model/Fscaler.pkl")
        feature = scaler.transform(feature)

        svmm5u = joblib.load("./model/F_SVM_model.pkl")
    elif Modes=="Mature_mRNA":
        print("feature extraction......")
        Data = feature_generate.get_features(FastaDATA, type_=Modes)
        print("feature selection......")
        x = feature_selection.select_features(Data, feature_index_M)
        feature = pd.DataFrame(x)
        scaler = joblib.load("./model/Mscaler.pkl")
        feature = scaler.transform(feature)

        svmm5u = joblib.load("./model/M_SVM_model.pkl")

    print("Predicting......")
    y_pred=svmm5u.predict(feature)
    y_pred_prob=svmm5u.predict_proba(feature)
    results=pd.DataFrame(np.zeros([len(y_pred),4]),columns=["Index","Sequences","Prediction","Confidence"])

    for i in range(len(y_pred)):
        if y_pred[i]==1:
            y_prob=str(round(y_pred_prob[i][1]*100,2))+"%"
            results.iloc[i,:]=[round(i+1),seq[i],"m5U",y_prob]
        if y_pred[i]==0:
            y_prob=str(round(y_pred_prob[i][0]*100,2))+"%"
            results.iloc[i,:]=[round(i+1),seq[i],"Non-m5U",y_prob]
    os.chdir("Results")
    results.to_csv(args.o, index=False)
    print("job finished!")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="it's usage tip.",description="m5U-SVM: Identification of RNA 5-methyluridine modification sites based on multi-view features of physicochemical features and distributed representation")
    parser.add_argument("-i", required=True, default=None, help="input fasta file")
    parser.add_argument("-m", required=True, default=None, help="sequecne mode Full_transcript or Mature_mRNA")
    parser.add_argument("-o", default="Results.csv", help="output a CSV results file")

    args = parser.parse_args()

    time_start = time.time()


    m5Usvm_pred(args.i,args.m)

    time_end = time.time()
    print('Total time cost', time_end - time_start, 'seconds')
