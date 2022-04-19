from PDB import PDB
from settings import PDB_Entry,RingData,MAXNUM,RAD
from get_shift import get_EF_shift,get_rc_shift,get_ring_shift,get_shift_from_file
from dataTransform import getThreeAA,getSimilarity,get_H_info
from ANNPyTorch import ANN_prediction,data_Normalization

import matplotlib.pyplot as plt 
import seaborn as sns

import torch
from torch import nn
from torch.autograd import Variable

import numpy as np
import pandas as pd

import os
import re
import math

atom_color = {'C':'pink','CA':'blue','CB':'cyan','N':'green','HN':'purple','HA':'orange'}

def main(atomName:str,OutlierName:dict):
    '''
    （未完成）Outlier值的特征分析，三氨基酸的模式等
    '''
    data = pd.read_csv('Working/Sparta+/data/train_set_'+atomName+'.csv',sep=',',header=0,dtype=np.float32)
    for key in data.columns:
        data.fillna({key:data[key].median()},inplace=True)
    min_threshold = data['shift1'].quantile(0.005)
    max_threshold = data['shift1'].quantile(0.995)
    OutlierData = data.loc[(data['shift1'] < min_threshold) & (data['shift1'] > max_threshold)]
    a = OutlierData['KeyName'].value
    
    for name in a:
        try:
            OutlierName[name]+=1
        except:
            OutlierName[name]=1
    


    #x = data.loc[:,(data.columns!='KeyName') | (data.columns!='shift')].values
    
def data_analysis(atomName):
    '''
    生成交叉数据，方便后续进行绘图
    '''
    data = pd.read_csv('Working/Sparta+/data/pred_set_'+atomName+'.csv',sep=',',header=0)
    data = data.loc[:,data.columns!='name']
    data = data.loc[:,data.columns!='shift1']
    data = data.astype(np.float32)

    for key in data.columns:
        data.fillna({key:data[key].median()},inplace=True)
    shift_name = 'shift2'
    
    shift_list1 = data.loc[:,data.columns==shift_name].values.tolist()
    data1 = (data.loc[:,data.columns!=shift_name]-data.loc[:,data.columns!=shift_name].min())/\
        (data.loc[:,data.columns!=shift_name].max()-data.loc[:,data.columns!=shift_name].min())
    x1 = data1.loc[:,data1.columns!=shift_name].values
    x1 = torch.from_numpy(x1)

    data = data_Normalization(data,shift_name)
    shift_list3 = data.loc[:,data.columns==shift_name].values.tolist()
    data2 = (data.loc[:,data.columns!=shift_name]-data.loc[:,data.columns!=shift_name].min())/\
        (data.loc[:,data.columns!=shift_name].max()-data.loc[:,data.columns!=shift_name].min())
    x2 = data2.loc[:,data2.columns!=shift_name].values
    x2 = torch.from_numpy(x2)

    #模型读取与预测
    ANN_Model = torch.load('Working/Sparta+/results/ANN_'+atomName+'_models.pth.tar')
    outputs = ANN_Model(x1)
    y_pred = outputs.detach().numpy().tolist()
    shift_list2 = y_pred.copy()
    outputs = ANN_Model(x2)
    y_pred = outputs.detach().numpy().tolist()
    shift_list4 = y_pred.copy()
    with open('Working/Sparta+/data/pred_set_shift_'+atomName+'.csv','w') as write_object:
        write_object.write('name,shift,label\n')
        for i in range(len(shift_list1)):
            write_object.write('{},{},{}\n'.format('experiment',shift_list1[i][0],'100%'))
        for i in range(len(shift_list2)):
            write_object.write('{},{},{}\n'.format('pred',shift_list2[i][0],'100%'))
        for i in range(len(shift_list3)):
            write_object.write('{},{},{}\n'.format('experiment',shift_list3[i][0],'99%'))
        for i in range(len(shift_list4)):
            write_object.write('{},{},{}\n'.format('pred',shift_list4[i][0],'99%'))

def plot_violin(atomName,y_name,data_path:str=''):
    '''
    绘制小提琴图
    '''
    df = pd.read_csv('Working/Sparta+/data/pred_set_shift_'+atomName+'.csv',sep=',',header=0)
    plt.figure(figsize=(16,12))
    sns.violinplot(x='label',y=y_name,hue='name',width=0.3,color=atom_color[atomName],data=df)
    plt.title('shift distribution between y with pred_y ('+atomName+')')
    plt.xlabel('samples don\'t drop out')
    plt.ylabel('shift value(ppm)')
    plt.legend()
    plt.savefig('Working/Sparta+/results/violin/figure_'+atomName+'_1.png')
    plt.close()


if __name__ == '__main__':
    atomNames = ['C','CA','CB','N','HN','HA']
    OutlierName = {}
    for atomName in atomNames:
        #main(atomName,OutlierName)
        #data_analysis(atomName)
        plot_violin(atomName,'shift')
    #print(OutlierName)














