#!/usr/bin/env python
# _*_coding:utf-8_*_
# @Time : 2022.10.20
# @Author : aochunyan
# @Email : acy196707@163.com
# @IDE : PyCharm
# @File : m5USVM.py

import re

def read_RNA_sequences(file):
    with open(file) as f:
        records = f.read()
    records = records.split('>')[1:]
    fasta_sequences = []
    sequence_name = []
    for fasta in records:
        array = fasta.split('\n')
        header, sequence = array[0].split()[0], re.sub('[^ACGTU-]', '-', ''.join(array[1:]).upper())
        header_array = header.split('|')
        name = header_array[0]
        sequence = re.sub('U', 'T', sequence)
        fasta_sequences.append(sequence)
        sequence_name.append(str(name))
    return fasta_sequences, sequence_name





