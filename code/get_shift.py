import numpy as np
import re
import math
from PDB import PDB

def get_shift_from_file(filepath:str):
    '''
    获取对应原子的shift值作为预测的y值
    '''
    shift_value = {}
    flag = 0
    with open(filepath) as read_object:
        for line in read_object:
            if re.search('_Atom_shift_assign_ID',line):
                flag = 1
                continue
            if flag:
                if re.search('stop_',line):
                    flag = 0
                    break
                info = line.strip().split(' ')
                if len(info)<8:
                    continue
                while '' in info:
                    info.remove('')
                key = int(info[1].strip())
                if key not in shift_value.keys():
                    shift_value[key] = {}
                atomname =  info[4].strip()
                if atomname == 'H':atomname = 'HN'
                shift_value[key][atomname] = float(info[6].strip())
    return shift_value

def get_rc_shift(resName:list,atomName:str):
    '''
    获取rc值
    '''
    #if 'MSE' or 'FTR' in resName:
    #    if resName[0]=='MSE' or resName[0]=='FTR':rcprev = 0
    #    if resName[1]=='MSE' or resName[1]=='FTR':rc,rcadj = 0,0
    #    if resName[2]=='MSE' or resName[2]=='FTR':rcnext = 0

    link_atom = {'HA':2,'CA':3,'CB':4,'C':5,'N':6,'HN':7}
    filepath = '/home/zhangjs/Working/Sparta+/parameters' 
    filename = ['rcprev.tab','randcoil.tab','rcadj.tab','rcnext.tab']
    for i,name in enumerate(filename):
        with open(filepath+'/'+name) as read_object:
            for line in read_object:
                if line.startswith('FORMAT'):
                    break
            read_object.readline()
            for line in read_object:
                info = line.strip().split(' ')
                while '' in info:
                    info.remove('')
                if i==0 and resName[0]==info[0].strip():
                    rcprev = info[link_atom[atomName]]
                elif i==1 and resName[1]==info[0].strip():
                    rc = info[link_atom[atomName]]
                elif i == 2 and resName[1]==info[0].strip():
                    rcadj = info[link_atom[atomName]]
                elif i==3 and resName[2]==info[0].strip():
                    rcnext = info[link_atom[atomName]]
    return float(rcprev),float(rc),float(rcadj),float(rcnext)

def get_ring_shift(pdb_object:PDB,chain_id:int=0,resID:int=1,atomName:str=''):
    '''
    计算环平面电子磁场效应
    '''
    if atomName == 'H' or atomName[:2] == 'HN':
        atomName = 'HN'
        if pdb_object.residueList[chain_id][resID] == 'PRO':
            return None
    elif atomName == 'HA' and pdb_object.residueList[chain_id][resID] == 'GLY':
        #ring_shift = pdb_object.getOrbitalShift(chain_id,resID,'HA3')
        atomName = 'HA2'
    try:
        ring_shift = pdb_object.getOrbitalShift(chain_id,resID,atomName)
        return ring_shift
    except:
        return None

def get_EF_shift(EF_shift:dict,resID:int,atomname:str):
    '''
    针对HA和HN原子计算EF_shift效应
    '''
    if atomname == 'HN' or atomname == 'HA':
        if resID not in EF_shift.keys():
            return None
        elif atomname in EF_shift[resID].keys():
            return EF_shift[resID][atomname]
    return None






