import numpy as np
import torch
from torch import nn
from torch.autograd import Variable
import pandas as pd

import math
import matplotlib.pyplot as plt 
import seaborn as sns

from ANNPyTorch import ANN_prediction,data_Normalization
from PDB import PDB
from settings import PDB_Entry,RingData,MAXNUM,RAD
from get_shift import get_EF_shift,get_rc_shift,get_ring_shift,get_shift_from_file
from dataTransform import getThreeAA,getSimilarity,get_H_info

import sys
import re

#绘图用颜色
atom_color = {'C':'pink','CA':'blue','CB':'cyan','N':'green','HN':'purple','HA':'orange'}

def validation(atomNames):
    '''
    用测试集以及训练集检测结果的好坏
    '''
    for atomName in atomNames:
        #数据输入与预处理
        data = pd.read_csv('Working/Sparta+/data/train_set_'+atomName+'.csv',sep=',',header=0)
        data = data.loc[:,data.columns!='name']
        data = data.astype(np.float32)
        for key in data.columns:
            data.fillna({key:data[key].median()},inplace=True)
        shift_name = 'shift1'
        data = data_Normalization(data,shift_name)
        y = data.loc[:,data.columns==shift_name].values
        data = data.loc[:,data.columns!='shift2']
        data = (data.loc[:,data.columns!=shift_name]-data.loc[:,data.columns!=shift_name].min())/\
            (data.loc[:,data.columns!=shift_name].max()-data.loc[:,data.columns!=shift_name].min())
        #drop_list = ['HB_1O_Energy','HB_2HN_Energy','HB_2HA_Energy','HB_2O_Energy','HB_3HN_Energy']
        #for feature_name in drop_list:
        #    data = data.loc[:,data.columns!=feature_name]
        x = data.loc[:,data.columns!=shift_name].values
        y = torch.from_numpy(y)
        x = torch.from_numpy(x)

        #模型读取与预测
        ANN_Model = torch.load('Working/Sparta+/results/ANN_'+atomName+'_models.pth.tar')
        outputs = ANN_Model(x)

        #验证集效果检测
        ErrorLoss = nn.MSELoss()
        dist1 = ErrorLoss(outputs,y)
        ErrorLoss = nn.L1Loss()
        dist2 = ErrorLoss(outputs,y)
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
        rmse = mse**0.5
        print(dist1.data,dist2.data,rmse)

        plt.figure(figsize=(8,8))
        plt.plot([min(plot_list[0]),max(plot_list[0])],\
            [min(plot_list[0]),max(plot_list[0])],\
            'r')
        plt.scatter(plot_list[0],plot_list[1],s=2,alpha=0.2,cmap='rainbow',color=atom_color[atomName])
        plt.title('correlation between y and y_pred('+atomName+')')
        plt.xlabel('shift_value_experiment y(ppm)')
        plt.ylabel('shift_value_pred y\'(ppm)')
        plt.savefig('Working/Sparta+/results/figure_'+atomName+'3.png')
        plt.close()

        plt.figure(figsize=(8,8))
        plt.hist(plot_list[2],bins=100,density=True,color=atom_color[atomName])
        plt.title('difference between y and y_pred('+atomName+')')
        plt.xlabel('shift_value experiment-pred (ppm)')
        plt.ylabel('percentage (%)')
        plt.savefig('Working/Sparta+/results/figure_'+atomName+'4.png')

def get_predict_data(pdbFilePath,writeFilePath,atomName):
    '''
    预测，针对包含H化的PDB结构进行化学位移预测（未完成）
    '''
    """
    MISSING_ATOm = {}
    with open(pdbFilePath) as read_object:
        flag = 0
        for line in read_object:
            if line.startswith('REMARK') and re.search('M RES CSSEQI  ATOMS',line):
                flag = 1
                continue
            if flag:
                info = line.strip().split(' ')
                while '' in info:
                    info.remove('')
                chain = info[3]
                resName = info[2]
                resID = info[4]
                pass##
    """


    pdb_object = PDB('/home/zhangjs/Toolkits/SPARTA+/SPARTA+/test/gb3.pdb')
    #pdb_object = PDB(pdbFilePath)

    chain_id = 0
    redict = pdb_object.residueList[chain_id]
    ThreeReList = getThreeAA(redict)
    tripeptides = []
    X = []
    shift = []

    for itA in ThreeReList:
        flag = 0
        if itA[0]!=pdb_object.r1 and itA[2]!=pdb_object.rN:
            for nu in range(itA[1]-2,itA[1]+3,1):
                if nu not in redict.keys():
                    flag = 1
                    break
        if flag:
            continue
        

        aa1 = pdb_object.getOneAAName(redict[itA[0]])
        aa2 = pdb_object.getOneAAName(redict[itA[1]])
        aa3 = pdb_object.getOneAAName(redict[itA[2]])
        ##获取相似性矩阵
        similarity_matrix=getSimilarity(aa1,aa2,aa3)

        ##获取Phi,Psi,Chi1,Chi2
        Angle_matrix = []
        for i in itA:
            Phi = pdb_object.getPhi(chain_id,i) if i!=pdb_object.r1 else 9999.000
            Psi = pdb_object.getPsi(chain_id,i) if i!=max(redict.keys()) else 9999.000
            Chi1 = pdb_object.getChi1(chain_id,i)
            Chi2 = pdb_object.getChi2(chain_id,i)
            Angle_matrix.append([Phi,Psi,Chi1,Chi2])
        phi_psi_matrix = []
        chi_matrix = []
        for i in range(len(Angle_matrix)):
            for j in range(len(Angle_matrix[0])):
                if j<=1:
                    phi_psi_matrix.append(math.cos(Angle_matrix[i][j]/RAD) if Angle_matrix[i][j]!=9999.0 else 0)
                    phi_psi_matrix.append(math.sin(Angle_matrix[i][j]/RAD) if Angle_matrix[i][j]!=9999.0 else 0)
                else:
                    if Angle_matrix[i][j] == 9999.0:
                        chi_matrix.append([0,0,0])
                    else:
                        chi_matrix.append([math.cos(Angle_matrix[i][j]/RAD),math.sin(Angle_matrix[i][j]/RAD),1])
        phi_psi_matrix = phi_psi_matrix[2:-2]


        ##获取HN_S2值
        HN_S2 = [0,0,0]
        for nu in range(3):
            try:HN_S2[nu] = pdb_object.HN_S2[itA[nu]]
            except:HN_S2[nu] = np.nan
        ##获取氢键的距离，角度，能量等值
        H_bond_info = []
        H_bond_info.append(get_H_info(pdb_object,chain_id,itA[0],'O'))
        H_bond_info.append(get_H_info(pdb_object,chain_id,itA[1],'HN'))
        H_bond_info.append(get_H_info(pdb_object,chain_id,itA[1],'HA'))
        H_bond_info.append(get_H_info(pdb_object,chain_id,itA[1],'O'))
        H_bond_info.append(get_H_info(pdb_object,chain_id,itA[2],'HN'))

        ##校正shift值
        RCprev,RC,RCadj,RCnext = get_rc_shift([redict[itA[0]],redict[itA[1]],redict[itA[2]]],atomName)
        shift_ring = get_ring_shift(pdb_object,chain_id,itA[1],atomName)
        shift_EF = get_EF_shift(pdb_object.ElectricField,itA[1],atomName)
        if shift_ring == None:
            shift_ring = 0
        if shift_EF == None:
            shift_EF = 0
        temp = []
        for mono_simi in similarity_matrix:
            temp.extend(mono_simi)
        temp.extend(phi_psi_matrix)
        for chi_num in chi_matrix:
            temp.extend(chi_num) 
        temp.extend(HN_S2)
        for mono_info in H_bond_info:
            temp.extend(mono_info)
        shift_delta.append(RC+RCadj+RCprev+RCnext,shift_ring,shift_EF,0)
        if atomName=='HN' and redict[itA[1]]=='PRO':
            continue
        elif atomName=='CB' and redict[itA[1]]=='GLY':
            continue
        X.append(temp.copy())
        shift.append(shift_delta.copy())
        tripeptides.append(aa1+aa2+aa3)
       
    
    ANN_Model = torch.load('Working/Sparta+/results/ANN_'+atomName+'_models.pth.tar')
    X = torch.from_numpy(X)
    y = ANN_Model(X)
    y = y.numpy()
    for i,num in enumerate(y):
        shift[i][-1] = num + shift[i][0] + shift[i][1] + shift[i][2]
    with open(writeFilePath,'w') as write_object:
        for i in range(len(shift)):
            write_object.write('{},{},{},{},{}\n'.format(tripeptides[i],shift[i][0],shift[i][1],shift[i][2],shift[i][3]))

if __name__=='__main__':
    atomNames = ['C','CA','CB','N','HN','HA']
    validation(atomNames)
    #get_predict_data()
