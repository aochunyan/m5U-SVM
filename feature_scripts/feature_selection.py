#!/usr/bin/env python
# _*_coding:utf-8_*_
# @Time : 2022.10.20
# @Author : aochunyan
# @Email : acy196707@163.com
# @IDE : PyCharm
# @File : m5USVM.py


import numpy as np
import pandas as pd

def select_features(allfeatures, feature_index):
    print('Feature selection...')
    new_features = []
    orignal_data = pd.DataFrame(allfeatures)
    for i in list(feature_index):
        new_features.append(orignal_data[int(i)])
    features = np.array(new_features).T

    return features








