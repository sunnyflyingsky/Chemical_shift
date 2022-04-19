from PDB import PDB
from settings import PDB_Entry,RingData,MAXNUM,RAD
from get_shift import get_EF_shift,get_rc_shift,get_ring_shift,get_shift_from_file
#from ANNPyTorch import ANN_prediction
import numpy as np
import os
import re
import math

def get_data_set(pdbFilePath:str,shiftFilePath:str,writeFilePath:str,atomName:str):
    '''
    获取数据集并存储为数据格式
    '''
    with open(pdbFilePath) as read_object:
        text = read_object.read()
        if re.search('MISSING ATOM',text):
            return None
    pdb_object = PDB(pdbFilePath)
    #pdb_object = PDB('/home/zhangjs/Toolkits/SPARTA+/SPARTA+/test/gb3.pdb')
    shift_value = get_shift_from_file(shiftFilePath)

    chain_id = 0
    redict = pdb_object.residueList[chain_id]
    ThreeReList = getThreeAA(redict)

    for itA in ThreeReList:
        flag = 0
        if itA[0]!=pdb_object.r1 and itA[2]!=pdb_object.rN:
            for nu in range(itA[1]-2,itA[1]+3,1):
                if nu not in redict.keys():
                    flag = 1
                    break
        if flag:
            continue
        

        ##shift
        if itA[1] in shift_value.keys() and atomName in shift_value[itA[1]].keys():
            shift = shift_value[itA[1]][atomName]
        else:
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
        #feature_x.append(temp.copy())
        #pred_y.append(shift-RC-RCadj-RCprev-RCnext-shift_ring-shift_EF)
        if atomName=='HN' and redict[itA[1]]=='PRO':
            continue
        elif atomName=='CB' and redict[itA[1]]=='GLY':
            continue
        with open(writeFilePath,'a') as write_object:
            write_object.write('{},'.format(aa1+aa2+aa3))
            for num in temp.copy():
                write_object.write('{},'.format(str(num)))
            write_object.write('{},{}\n'.format(shift-RC-RCadj-RCprev-RCnext-shift_ring-shift_EF,shift-RC-RCadj-RCprev-RCnext))

    return 1

def getThreeAA(seq:dict):
    '''
    三氨基酸滑窗
    '''
    relist = []

    Temp = []
    for reNum,reName in seq.items():
        if len(Temp)<3:
            Temp.append(reNum)
        else:
            relist.append(Temp.copy())
            del Temp[0]
            Temp.append(reNum)
    relist.append(Temp.copy())
        
    return relist

def getSimilarity(aa1,aa2,aa3):
    '''
    读取氨基酸相似性
    '''
    homology_matrix = [[],[],[]]
    i = 0
    with open('Working/Sparta+/parameters/homology.tab') as read_object:
        for line in read_object:
            if i < 7 :
                i+=1
            elif not line.startswith(' '):
                num_info = line.strip().split('  ')
                for j,num in enumerate(num_info):
                    if j==0:
                        if str(num)==aa1:
                            flag = 0
                        elif str(num)==aa2:
                            flag = 1
                        elif str(num)==aa3:
                            flag = 2
                        else:
                            break
                    if j==0 or j==5 or j == 10:
                        continue
                    else:
                        homology_matrix[flag].append(int(num))
    return homology_matrix

def get_H_info(pdb_object:PDB,chainID:int,resID:int,atomName:str):
    '''
    '''
    H_info_list = []
    if atomName == 'HA' and atomName not in pdb_object.Conformers_Atom[chainID][resID].keys():
        atomName = 'HA1'
    if atomName not in pdb_object.HBDistList[chainID][resID].keys():
        return [0,0,0,0,0]
    dist = pdb_object.HBDistList[chainID][resID][atomName]
    energy = pdb_object.HBEnergyList[chainID][resID][atomName]
    try:
        donorAngle = pdb_object.HB_DHO_AngleList[chainID][resID][atomName]
        acceptorAngle = pdb_object.HB_HOA_AngleList[chainID][resID][atomName]
    except:
        return [0,0,0,0,0]
    H_info_list = [dist,math.cos(donorAngle/RAD),math.cos(acceptorAngle/RAD),energy,1].copy()  #,energy
    return H_info_list

def Angle_correct(angle:float):
    '''
    针对特殊角度的校正，将这些角度归一化为固定的离散值
    '''
    if angle<=60:
        pass
    elif angle>60 and angle<=120:
        pass
    elif angle>120 and angle<=180:
        pass
    return angle

def pred_set_(atomName:str):
    '''
    生成验证数据集
    '''
    #输出文件路径
    writeFilePath = '/home/zhangjs/Working/Sparta+/data/pred_set_' + atomName +'.csv'

    feature_title = []
    mono = 'S'
    for i in range(60):
        feature_title.append(mono+str(i))
    with open(writeFilePath,'w') as write_object:
        write_object.write('name,'+\
            ','.join(feature_title)+\
            ',Psi_1_cos,Psi_1_sin,Phi_2_cos,Phi_2_sin,Psi_2_cos,Psi_2_sin,Phi_3_cos,Phi_3_sin'+\
            ',Chi1_1_cos,Chi1_1_sin,chi1_1_Exist,Chi2_1_cos,Chi2_1_sin,chi2_1_Exist'+\
            ',Chi1_2_cos,Chi1_2_sin,chi1_2_Exist,Chi2_2_cos,Chi2_2_sin,chi2_2_Exist'+\
            ',Chi1_3_cos,Chi1_3_sin,chi1_3_Exist,Chi2_3_cos,Chi2_3_sin,chi2_3_Exist'+\
            ',HN_S2_1,HN_S2_2,HN_S2_3'+\
            ',HB_1O_dist,HB_1O_DonorAngle,HB_1O_AcceptorAngle,HB_1O_Energy,HB_1O_Exist'+\
            ',HB_2HN_dist,HB_2HN_DonorAngle,HB_2HN_AcceptorAngle,HB_2HN_Energy,HB_2HN_Exist'+\
            ',HB_2HA_dist,HB_2HA_DonorAngle,HB_2HA_AcceptorAngle,HB_2HA_Energy,HB_2HA_Exist'+\
            ',HB_2O_dist,HB_2O_DonorAngle,HB_2O_AcceptorAngle,HB_2O_Energy,HB_2O_Exist'+\
            ',HB_3HN_dist,HB_3HN_DonorAngle,HB_3HN_AcceptorAngle,HB_3HN_Energy,HB_3HN_Exist'+\
            ',shift1,shift2\n')
            #,HB_1O_Energy,HB_2HN_Energy,HB_2HA_Energy,HB_2O_Energy,HB_3HN_Energy

    ##testset
    folderPath = 'Working/Sparta+/data/shiftx2-testset-June2011' #testset trainset
    pdb_dir_list = os.listdir(folderPath + '/PDB-testset-addHydrogens') #testset training
    shift_dir_list = os.listdir(folderPath + '/CS-corrected-testset-addPDBresno')

    dir_list_dict = {}

    #包装对应每个pdb文件和shift文件
    for fileName1,fileName2 in zip(pdb_dir_list,shift_dir_list):
        try:
            dir_list_dict[fileName1[:4]][0] = fileName1
        except:
            dir_list_dict[fileName1[:4]] = [fileName1,'']
        try:
            dir_list_dict[fileName2[:4]][1] = fileName2
        except:
            dir_list_dict[fileName2[:4]] = ['',fileName2]
    for value in dir_list_dict.values():
        ##for debug
        #if value[0]!='R202_2TRXA.pdbH':
        #    continue
        fileName1,fileName2 = value[0],value[1]
        
        pdbFilePath = folderPath+'/PDB-testset-addHydrogens/' + fileName1
        shiftFilePath = folderPath+ '/CS-corrected-testset-addPDBresno/'+fileName2
        res = get_data_set(pdbFilePath,shiftFilePath,writeFilePath,atomName)
        if res != None:
            print('one step finished!')
        else:print(fileName1)

def train_set_(atomName:str):
    '''
    生成测试数据集
    '''
    #输出文件路径
    writeFilePath = '/home/zhangjs/Working/Sparta+/data/train_set_' + atomName +'.csv'

    feature_title = []
    mono = 'S'
    for i in range(60):
        feature_title.append(mono+str(i))
    with open(writeFilePath,'w') as write_object:
        write_object.write('name,'+\
            ','.join(feature_title)+\
            ',Psi_1_cos,Psi_1_sin,Phi_2_cos,Phi_2_sin,Psi_2_cos,Psi_2_sin,Phi_3_cos,Phi_3_sin'+\
            ',Chi1_1_cos,Chi1_1_sin,chi1_1_Exist,Chi2_1_cos,Chi2_1_sin,chi2_1_Exist'+\
            ',Chi1_2_cos,Chi1_2_sin,chi1_2_Exist,Chi2_2_cos,Chi2_2_sin,chi2_2_Exist'+\
            ',Chi1_3_cos,Chi1_3_sin,chi1_3_Exist,Chi2_3_cos,Chi2_3_sin,chi2_3_Exist'+\
            ',HN_S2_1,HN_S2_2,HN_S2_3'+\
            ',HB_1O_dist,HB_1O_DonorAngle,HB_1O_AcceptorAngle,HB_1O_Energy,HB_1O_Exist'+\
            ',HB_2HN_dist,HB_2HN_DonorAngle,HB_2HN_AcceptorAngle,HB_2HN_Energy,HB_2HN_Exist'+\
            ',HB_2HA_dist,HB_2HA_DonorAngle,HB_2HA_AcceptorAngle,HB_2HA_Energy,HB_2HA_Exist'+\
            ',HB_2O_dist,HB_2O_DonorAngle,HB_2O_AcceptorAngle,HB_2O_Energy,HB_2O_Exist'+\
            ',HB_3HN_dist,HB_3HN_DonorAngle,HB_3HN_AcceptorAngle,HB_3HN_Energy,HB_3HN_Exist'+\
            ',shift1,shift2\n')
            #,HB_1O_Energy,HB_2HN_Energy,HB_2HA_Energy,HB_2O_Energy,HB_3HN_Energy

    ##trainset
    folderPath = 'Working/Sparta+/data/shiftx2-trainset-June2011' #testset trainset
    pdb_dir_list = os.listdir(folderPath + '/PDB-training-addHydrogens') #testset training
    shift_dir_list = os.listdir(folderPath + '/CS-corrected-training-addPDBresno')

    dir_list_dict = {}

    #包装对应每个pdb文件和shift文件
    for fileName1,fileName2 in zip(pdb_dir_list,shift_dir_list):
        try:
            dir_list_dict[fileName1[:4]][0] = fileName1
        except:
            dir_list_dict[fileName1[:4]] = [fileName1,'']
        try:
            dir_list_dict[fileName2[:4]][1] = fileName2
        except:
            dir_list_dict[fileName2[:4]] = ['',fileName2]
    for value in dir_list_dict.values():
        ##for debug
        #if value[0]!='R202_2TRXA.pdbH':
        #    continue
        fileName1,fileName2 = value[0],value[1]
        pdbFilePath = folderPath+'/PDB-training-addHydrogens/' + fileName1
        shiftFilePath = folderPath+ '/CS-corrected-training-addPDBresno/'+fileName2
       
        res = get_data_set(pdbFilePath,shiftFilePath,writeFilePath,atomName)
        if res != None:
            print('one step finished!')
        else:print(fileName1)

if __name__ == '__main__':
    atomNameList = ['CB','CA','HN', 'N', 'HA', 'C']
    for atomName in atomNameList:
        pred_set_(atomName)
        train_set_(atomName)
    


    




