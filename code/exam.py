import math
import numpy as np
import torch
from torch import nn
from torch.autograd import Variable
import pandas as pd
from ANNPyTorch import ANN_prediction
#import pdblib

#exam1 = pdb.Pdb('/home/zhangjs/Working/Sparta+/data/shiftx2-trainset-June2011/PDB-training/R001_1QRXA.pdb')
#b = exam1.get_names()
#print(b)
with open('Working/Sparta+/results/RMSD.txt','w') as write_object:
    atomName = ['N','CA','C','HN','HA']
    for atom in atomName:
        data = pd.read_csv('Working/Sparta+/data/pred_set_'+atom+'.csv',sep=',',header=0,dtype=np.float32)
        for key in data.columns:
            data.fillna({key:data[key].median()},inplace=True)

        y = data.loc[:,data.columns=='shift'].values
        data = (data.loc[:,data.columns!='shift']-data.loc[:,data.columns!='shift'].min())/\
            (data.loc[:,data.columns!='shift'].max()-data.loc[:,data.columns!='shift'].min())

        x = data.loc[:,data.columns!='shift'].values
        y = torch.from_numpy(y)
        x = torch.from_numpy(x)
        ANN_Model = torch.load('Working/Sparta+/results/ANN_'+atom+'_models.pth.tar')
        outputs = ANN_Model(x)
        y = y.numpy().tolist()
        y_pred = outputs.detach().numpy().tolist()
        plot_list = [[],[],[]]
        for i,j in zip(y,y_pred):
            if i[0]>80:continue
            plot_list[0].append(i[0])
            plot_list[1].append(j[0])
            plot_list[2].append(i[0]-j[0])
        a = [x**2 for x in plot_list[2]]
        mse = np.sum(a)/len(plot_list[2])
        rmse1 = mse**0.5

        data = pd.read_csv('Working/Sparta+/data/train_set_'+atom+'.csv',sep=',',header=0,dtype=np.float32)
        for key in data.columns:
            data.fillna({key:data[key].median()},inplace=True)

        y = data.loc[:,data.columns=='shift'].values
        data = (data.loc[:,data.columns!='shift']-data.loc[:,data.columns!='shift'].min())/\
            (data.loc[:,data.columns!='shift'].max()-data.loc[:,data.columns!='shift'].min())

        x = data.loc[:,data.columns!='shift'].values
        y = torch.from_numpy(y)
        x = torch.from_numpy(x)
        ANN_Model = torch.load('Working/Sparta+/results/ANN_'+atom+'_models.pth.tar')
        outputs = ANN_Model(x)
        y = y.numpy().tolist()
        y_pred = outputs.detach().numpy().tolist()
        plot_list = [[],[],[]]
        for i,j in zip(y,y_pred):
            if i[0]>80:continue
            plot_list[0].append(i[0])
            plot_list[1].append(j[0])
            plot_list[2].append(i[0]-j[0])
        a = [x**2 for x in plot_list[2]]
        mse = np.sum(a)/len(plot_list[2])
        rmse2 = mse**0.5
        write_object.write('{},{}\n'.format(rmse1,rmse2))



