from settings import PDB_Entry,RingData,MAXNUM,RAD
import math
import re



class PDB(object):
    '''
    用于加载PDB结构，这是一个自定义的PDB结构解析器，如果你有更好的选择，比如调包，可以尝试更可靠的方式
    '''
    def __init__(self,filePath:str):
        self.ATOM_info = []  #简单序列排列的原子列表
        self.sequenceList = {}  #氨基酸链
        self.Conformers = []  ##chain_id->atomnum->ATOM_PDB
        self.Conformers_Atom = []  ##chain_id->resnum->atomname->ATOM_PDB
        self.residueList = []  ##chain_id->resnum->resname
        self.EMPTY = PDB_Entry()  #储存空的原子结构
        itnum = 0
        flag = 0
        #RemarkNum = None
        #MissingAtomResid = {}
        with open(filePath) as read_object:
            for line in read_object:
                """
                if line.startswith('REMARK') and re.search('M RES CSSEQI  ATOMS',line):
                    info = line.strip().split(' ')
                    RemarkNum = info[1]
                    continue
                if RemarkNum and line.startswith('REMARK '+RemarkNum):
                    info = line.strip().split(' ')
                    while '' in info:
                        info.remove('')
                    chain = info[3]
                    resID = int(info[4])
                    if chain not in MissingAtomResid.keys():
                        MissingAtomResid[chain] = []
                    MissingAtomResid[chain].append(resID)
                """
                if line.startswith('ATOM'):
                #   flag = 1
                #if flag and line.startswith('TER'):
                #    flag = 0
                #if flag and not line.startswith('ANISOU'):
                    itnum+=1
                    num = int(self.getField(line,2))
                    name = self.getField(line,3)
                    if name == 'H' or name=='1H' or name=='2H' or name=='3H':
                        name = name + 'N'
                    if name[0].isdigit():
                        name = name[1:] + name[0]
                    residue = self.getField(line,4)
                    chain_id = self.getField(line,5)
                    seq = int(self.getField(line,6))
                    X = float(self.getField(line,7))
                    Y = float(self.getField(line,8))
                    Z = float(self.getField(line,9))
                    #occupancy = getField(line.strip(),10)
                    b_factor = float(self.getField(line.strip(),11))
                    pdb_temp = PDB_Entry(atomNum=num,resNum=seq,atomName=name,resName=residue,\
                        chainName=chain_id,X=X,Y=Y,Z=Z,B_Factor=b_factor,Coord=[X,Y,Z])
                    self.ATOM_info.append(pdb_temp)
                    #if seq in MissingAtomResid[chain_id]:
                    #    continue
                    if chain_id in self.sequenceList.keys():
                        self.sequenceList[chain_id]+=self.getOneAAName(residue)
                        self.Conformers[-1][itnum] = PDB_Entry(atomNum=itnum,resNum=seq,atomName=name,resName=residue,\
                            chainName=chain_id,X=X,Y=Y,Z=Z,B_Factor=b_factor,Coord=[X,Y,Z])
                        if seq not in self.Conformers_Atom[-1].keys():
                            self.Conformers_Atom[-1][seq] = {}
                        self.Conformers_Atom[-1][seq][name] = PDB_Entry(atomNum=itnum,resNum=seq,atomName=name,resName=residue,\
                            chainName=chain_id,X=X,Y=Y,Z=Z,B_Factor=b_factor,Coord=[X,Y,Z])
                        self.residueList[-1][seq] = residue
                    else:
                        self.sequenceList[chain_id]=self.getOneAAName(residue)
                        self.Conformers.append({})
                        self.Conformers_Atom.append({})
                        self.Conformers[-1][itnum] = PDB_Entry(atomNum=itnum,resNum=seq,atomName=name,resName=residue,\
                            chainName=chain_id,X=X,Y=Y,Z=Z,B_Factor=b_factor,Coord=[X,Y,Z])
                        if seq not in self.Conformers_Atom[-1].keys():
                            self.Conformers_Atom[-1][seq] = {}
                        self.Conformers_Atom[-1][seq][name] = PDB_Entry(atomNum=itnum,resNum=seq,atomName=name,resName=residue,\
                            chainName=chain_id,X=X,Y=Y,Z=Z,B_Factor=b_factor,Coord=[X,Y,Z])
                        self.residueList.append({seq:residue})
                elif line.startswith(''):
                    pass
        #判断二硫键的形成
        for Cid in range(len(self.residueList)):
            for it in self.residueList[Cid].keys():
                if self.residueList[Cid][it] == 'CYS':
                    if self.isSSBonded(Cid,it):
                        self.residueList[Cid][it] = 'cys'
        self.acceptorList = []  ##chain_id->atomnum->atomnum
        self.donorList = []   ##chain_id->atomnum->atomnum
        self.HBDistList = []  ##chain_id->resid->atomnum->dist
        self.HBEnergyList = []  ##chain_id->resid->atomnum->energy
        self.HB_DHO_AngleList = []  ##chain_id->resid->atomnum->angle
        self.HB_HOA_AngleList = []  ##chain_id->resid->atomnum->angle
        for i in range(len(self.Conformers)):
            self.acceptorList.append({})
            self.donorList.append({})
            self.HBDistList.append({})
            self.HBEnergyList.append({})
            self.HB_DHO_AngleList.append({})
            self.HB_HOA_AngleList.append({})
            for resNum_key,it_value in self.Conformers_Atom[i].items():
                self.HB_DHO_AngleList[i][resNum_key] = {}
                self.HB_HOA_AngleList[i][resNum_key] = {}
                self.HBEnergyList[i][resNum_key] = {}
                self.HBDistList[i][resNum_key] = {}
                for atomName in it_value.keys():
                    if atomName[:2] == 'HN':
                        atomName = 'HN'
                    self.HBDistList[i][resNum_key][atomName] = 0
                    self.HBEnergyList[i][resNum_key][atomName] = 0
            #初始化氢键距离，氢键拐角
            self.initHBond(i)
        
        ##计算HN_S2和EF效应
        self.r1 = min(list(self.residueList[0].keys()))
        self.rN = max(list(self.residueList[0].keys()))
        self.HN_S2 = {}
        self.ElectricField = {}
        self.calc_HN_S2()
        self.calc_ElectricField()

        self.Rings = []
        self.RingNo = 0
        self.initOrbitalShift()

    #####
    # 载入PDB参数
    #####

    def getField(self,text:str,index:int):
        '''
        从PDB的原子field信息中获取对应的具体信息的方法
        '''
        if index==1:return text[0:6].strip()  ##"ATOM"
        elif index==2:return text[6:11].strip()  ##Atom number
        elif index==3:return text[11:16].strip()  ##Atom name
        elif index==4:return text[17:21].strip()  ##Residue name
        elif index==5:return text[21:22].strip()  ##Chain ID
        elif index==6:return text[22:26].strip()  ##Residue seq
        elif index==7:return text[30:38].strip()  ##X
        elif index==8:return text[38:46].strip()  ##Y
        elif index==9:return text[46:54].strip()  ##Z
        elif index==10:return text[54:60].strip()  ##Occupancy
        elif index==11:return text[60:66].strip()  ##B-factor
        else:
            return None
    
    def getEntry(self,conformerID:int,rNum:int,aName:str=''):
        '''
        获取对应的原子
        '''
        if aName=='':
            return self.Conformers[conformerID][rNum]
        else:
            return self.Conformers_Atom[conformerID][rNum][aName]
    
    #####
    # 计算氨基酸的键转角,
    #####
    def getBondAngle(self, A: list, B: list, C: list):
        '''
        得到键键之间的夹角
        获得的是角度值
        '''
        v,v1,v2 = B.copy(),A.copy(),C.copy()
        self.Vec3Sub(v1,v) #向量BA
        self.Vec3Sub(v2,v) #向量BC
        self.Vec3Norm(v1)
        self.Vec3Norm(v2)
        c = self.Vec3Scalar(v1,v2)  #内积
        self.Vec3Cross(v1,v2)  #外积
        s = self.Vec3Abs(v1)

        return math.atan2(s,c)*RAD
    
    def getBondAngle_2(self, a: PDB_Entry, b: PDB_Entry, c:PDB_Entry):
        '''
        函数重载，获取原子之间的键角
        '''
        return self.getBondAngle(a.Coord,b.Coord,c.Coord)
    
    def getDihedralAngle(self,a: PDB_Entry,b: PDB_Entry,c:PDB_Entry,d:PDB_Entry):
        '''
        求平面夹角
        返回的是角度值
        '''
        cb,n1,n2 = [0,0,0],[0,0,0],[0,0,0]
        co = 0
        if a.atomName == '' or  b.atomName == '' or c.atomName == '' or  d.atomName == '':
            return MAXNUM
        #bc向量
        cb[0],cb[1],cb[2] = c.X-b.X, c.Y-b.Y,c.Z-b.Z
        ##平面abc的法向量
        n1[0] = (b.Y - a.Y) * cb[2] + (a.Z - b.Z) * cb[1]
        n1[1] = (b.Z - a.Z) * cb[0] + (a.X - b.X) * cb[2]
        n1[2] = (b.X - a.X) * cb[1] + (a.Y - b.Y) * cb[0]
        ##平面bcd的法向量
        n2[0] = cb[1] * (d.Z - c.Z) + cb[2] * (c.Y - d.Y)
        n2[1] = cb[2] * (d.X - c.X) + cb[0] * (c.Z - d.Z)
        n2[2] = cb[0] * (d.Y - c.Y) + cb[1] * (c.X - d.X)
        TEMP,TEMP1,TEMP2,TEMP3,TEMP4,TEMP5 = \
            n1[0],n1[1],n1[2],n2[0],n2[1],n2[2]
        #计算两个平面的夹角，即abc平面与bcd平面的夹角的cos值
        #cos(a-b-c-d) = |n1·n2|/||*||
        co = float((n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2])/\
            math.sqrt((TEMP**2+TEMP1**2+TEMP2**2)*(TEMP3**2+TEMP4**2+TEMP5**2)))
        return RAD * (math.acos(co) * \
            self.sgn(float((n1[1]*n2[2]-n1[2]*n2[1])*cb[0] + \
            (n1[2]*n2[0]-n1[0]*n2[2])*cb[1] + \
            (n1[0]*n2[1]-n1[1]*n2[0])*cb[2])))
        ##这个正负取得是什么值，两个法向量与BC垂直，则

    def getPhi(self,conformerID:int,resNum:int):
        '''
        '''
        return self.getDihedralAngle(self.getEntry(conformerID,resNum-1,'C'),\
            self.getEntry(conformerID,resNum,'N'),\
            self.getEntry(conformerID,resNum,'CA'),\
            self.getEntry(conformerID,resNum,'C'))
      
    def getPsi(self,conformerID:int,resNum:int):
        '''
        '''
        return self.getDihedralAngle(self.getEntry(conformerID,resNum,'N'),\
            self.getEntry(conformerID,resNum,'CA'),\
            self.getEntry(conformerID,resNum,'C'),\
            self.getEntry(conformerID,resNum+1,'N'))

    def getOmega(self,conformerID:int,resNum:int):
        '''
        '''
        return self.getDihedralAngle(self.getEntry(conformerID,resNum-1,'CA'),\
            self.getEntry(conformerID,resNum-1,'C'),\
            self.getEntry(conformerID,resNum,'N'),\
            self.getEntry(conformerID,resNum,'CA'))
    
    def getChi1(self,conformerID:int,resNum:int):
        '''
        '''
        rName = self.residueList[conformerID][resNum]
        if rName == 'GLY' or rName == 'ALA':
            return MAXNUM
        elif rName == 'VAL':
            R = self.getEntry(conformerID,resNum,'CG2')
        elif rName == 'ILE':
            R = self.getEntry(conformerID,resNum,'CG1')
        elif rName == 'THR':
            R = self.getEntry(conformerID,resNum,'OG1')
        elif rName == 'SER':
            R = self.getEntry(conformerID,resNum,'OG')
        elif rName == 'CYS' or rName == 'cys':
            R = self.getEntry(conformerID,resNum,'SG')
        else:
            R = self.getEntry(conformerID,resNum,'CG')

        return self.getDihedralAngle(self.getEntry(conformerID,resNum,'N'),\
            self.getEntry(conformerID,resNum,'CA'),\
            self.getEntry(conformerID,resNum,'CB'),R)

    def getChi2(self,conformerID:int,resNum:int):
        '''
        '''
        rName = self.residueList[conformerID][resNum]
        if rName=='GLY' or rName=='ALA':
            return MAXNUM

        if rName == 'VAL':
            R3 = self.getEntry(conformerID,resNum,'CG1')
            try:R4 = self.getEntry(conformerID,resNum,'HG11')
            except:return MAXNUM
        elif rName == 'THR':
            R3 = self.getEntry(conformerID,resNum,'CG2')
            try:R4 = self.getEntry(conformerID,resNum,'HG21')
            except:return MAXNUM
        elif rName == 'SER':
            R3 = self.getEntry(conformerID,resNum,'OG')
            try:R4 = self.getEntry(conformerID,resNum,'HG')
            except:return MAXNUM
        elif rName == 'ILE':
            R3 = self.getEntry(conformerID,resNum,'CG1')
            try:R4 = self.getEntry(conformerID,resNum,'CD1')
            except:return MAXNUM
        elif (rName == 'CYS' or rName == 'cys') and 'HG' in self.Conformers_Atom[conformerID][resNum]:
            R3 = self.getEntry(conformerID,resNum,'SG')
            R4 = self.getEntry(conformerID,resNum,'HG')
        elif (rName == 'CYS' or rName == 'cys') and 'HG' not in self.Conformers_Atom[conformerID][resNum]:
            return MAXNUM
        else:
            R3 = self.getEntry(conformerID,resNum,'CG')
            if rName=='ASN' or rName == 'ASP':
                try:
                    R4 = self.getEntry(conformerID,resNum,'OD1')
                except:return MAXNUM
            elif rName=='HIS':
                R4 = self.getEntry(conformerID,resNum,'ND1')
            elif rName=='MET':
                R4 = self.getEntry(conformerID,resNum,'SD')
            elif rName=='LEU' or rName=='PHE' or rName=='TYR' or rName=='TRP' or rName== 'FTR':
                R4 = self.getEntry(conformerID,resNum,'CD1')
            elif rName=='MSE':
                R4 = self.getEntry(conformerID,resNum,'SE')
            else:
                R4 = self.getEntry(conformerID,resNum,'CD')

        return self.getDihedralAngle(self.getEntry(conformerID,resNum,'CA'),\
            self.getEntry(conformerID,resNum,'CB'),\
            R3,R4)

    def getDist(self,A:list,B:list):
        '''
        计算距离
        '''
        v = [0,0,0]
        v = B.copy()
        self.Vec3Sub(v,A)
        return self.Vec3Abs(v)
    
    def isSSBonded(self,conformerID:int,resNum:int):
        '''
        判断SG是否形成了二硫键
        '''
        if self.residueList[conformerID][resNum]!='CYS' and self.residueList[conformerID][resNum] !='cys':
            return False
        CYS_SG = self.getEntry(conformerID,resNum,'SG')
        for i in range(len(self.residueList[conformerID][resNum])):
            if (self.residueList[conformerID][resNum][i]=='CYS' or self.residueList[conformerID][resNum][i]=='cys') and abs(i-resNum)>=4:
                if self.getDist(CYS_SG.Coord,self.Conformers_Atom[0][i]['SG'].Coord)<=2.5:
                    return True
        return False

    #####
    # 获取氢键的信息
    #####
    def initHBond(self,chainNum:int=0):
        '''
        '''
        conf = self.Conformers[chainNum]
        for i in conf.keys():
            a = self.isAcceptor(conf[i],chainNum)
            if a.atomName!='':
                self.acceptorList[chainNum][i] = a.atomNum
            else:
                temp = self.isDonor(conf[i],chainNum)
                if temp.atomName=='': continue
                self.donorList[chainNum][i] = temp.atomNum
        for itA in self.acceptorList[chainNum].keys():
            A = conf[itA]
            A_Heavy = conf[self.acceptorList[chainNum][itA]]
            for itD in self.donorList[chainNum].keys():
                D = conf[itD]
                D_Heavy = conf[self.donorList[chainNum][itD]]

                if abs(A.resNum - D.resNum) < 2: continue

                D_ON = self.getDist(A.Coord,D_Heavy.Coord) 
                D_CH = self.getDist(A_Heavy.Coord,D.Coord)
                D_OH = self.getDist(A.Coord,D.Coord)
                D_CN = self.getDist(A_Heavy.Coord,D_Heavy.Coord)

                if D_OH > D_CN : continue

                HBond_E = 332.0*0.42*0.20
                HBond_E *= (1.0/D_ON + 1.0/D_CH - 1.0/D_OH - 1.0/D_CN)

                if HBond_E < -0.5:
                    aN = 'HN' if A.atomName[:2] == 'HN' else A.atomName
                    dN = 'HN' if D.atomName[:2] == 'HN' else D.atomName
                    if self.HBDistList[chainNum][A.resNum][aN] == 0 or \
                        self.HBDistList[chainNum][A.resNum][aN] > D_OH or \
                        self.HBEnergyList[chainNum][A.resNum][aN] > HBond_E:
                        self.HBDistList[chainNum][A.resNum][aN] = D_OH
                        self.HB_DHO_AngleList[chainNum][A.resNum][aN] = self.getBondAngle_2(D_Heavy,D,A)
                        self.HB_HOA_AngleList[chainNum][A.resNum][aN] = self.getBondAngle_2(D,A,A_Heavy)
                        self.HBEnergyList[chainNum][A.resNum][aN] = HBond_E
                    if self.HBDistList[chainNum][D.resNum][dN] == 0 or \
                        self.HBDistList[chainNum][D.resNum][dN] > D_OH or \
                        self.HBEnergyList[chainNum][D.resNum][dN] > HBond_E:
                        self.HBDistList[chainNum][D.resNum][dN] = D_OH
                        self.HB_DHO_AngleList[chainNum][D.resNum][dN] = self.getBondAngle_2(D_Heavy,D,A)
                        self.HB_HOA_AngleList[chainNum][D.resNum][dN] = self.getBondAngle_2(D,A,A_Heavy)
                        self.HBEnergyList[chainNum][D.resNum][dN] = HBond_E
        
    def getHBondDist(self,resNum:int,atomName:str):
        '''
        获取氢键的长度
        '''
        return self.HBDistList[resNum][atomName]

    def isAcceptor(self,A:PDB_Entry,id:int=0):
        '''
        判断原子是否是H的acceptor
        '''
        if A.atomName=='O' or A.atomName=='OT1' or A.atomName=='OT2':
            return self.Conformers_Atom[id][A.resNum]['C']
        if A.resName == 'ASN' and A.atomName == 'OD1':
            return self.Conformers_Atom[id][A.resNum]['CG']
        if A.resName == 'ASP' and (A.atomName == 'OD1' or A.atomName == 'OD2'):
            return self.Conformers_Atom[id][A.resNum]['CG']
        if (A.resName == 'CYS' or A.resName == 'cys') and A.atomName == 'SG':
            return self.Conformers_Atom[id][A.resNum]['CB']
        if A.resName == 'GLN' and A.atomName == 'OE1':
            return self.Conformers_Atom[id][A.resNum]['CD']
        if A.resName == 'GLU' and (A.atomName=='OE1' or A.atomName=='OE2'):
            return self.Conformers_Atom[id][A.resNum]['CD']
        if A.resName == 'HIS' and (A.atomName == 'ND1' or A.atomName=='NE2'):
            return self.Conformers_Atom[id][A.resNum]['CE1']
        if A.resName == 'MET' and A.atomName == 'SD':
            return self.Conformers_Atom[id][A.resNum]['CG']
        if A.resName == 'SER' and A.atomName == 'OG':
            return self.Conformers_Atom[id][A.resNum]['CB']
        if A.resName == 'THR' and A.atomName == 'OG1':
            return self.Conformers_Atom[id][A.resNum]['CB']
        if A.resName == 'TYR' and A.atomName == 'OH':
            return self.Conformers_Atom[id][A.resNum]['CZ']
        
        return self.EMPTY

    def isDonor(self,D:PDB_Entry,id:int=0):
        '''
        判断原子是否是H的donor
        '''
        if re.search('HN',D.atomName) != None:
            return self.Conformers_Atom[id][D.resNum]['N']
        if re.search('HA',D.atomName) != None:
            return self.Conformers_Atom[id][D.resNum]['CA']
        if D.resName == 'ARG':
            if re.search('HE',D.atomName) != None:
                return self.Conformers_Atom[id][D.resNum]['NE']
            if re.search('HH1',D.atomName) != None:
                return self.Conformers_Atom[id][D.resNum]['NH1']
            if re.search('HH2',D.atomName) != None:
                return self.Conformers_Atom[id][D.resNum]['NH2']
        if D.resName == 'ASP' and re.search('HD2',D.atomName) != None:
            return self.Conformers_Atom[id][D.resNum]['ND2']
        if (D.resName == 'CYS' or D.resName == 'cys') and re.search('HG',D.atomName) != None:
            return self.Conformers_Atom[id][D.resNum]['SG']
        if D.resName == 'GLU' and re.search('HE2',D.atomName) != None:
            return self.Conformers_Atom[id][D.resNum]['NE2']
        if D.resName == 'HIS':
            if re.search('HD1',D.atomName) != None:
                return self.Conformers_Atom[id][D.resNum]['ND1']
            if re.search('HE2',D.atomName) != None:
                return self.Conformers_Atom[id][D.resNum]['NE2']
        if D.resName == 'LYS' and (D.atomName=='HZ1' or D.atomName=='HZ2' or D.atomName=='HZ3'):
            return self.Conformers_Atom[id][D.resNum]['NZ']
        if D.resName == 'SER' and D.atomName == 'HG':
            return self.Conformers_Atom[id][D.resNum]['OG']
        if D.resName == 'THR' and D.atomName == 'HG1':
            return self.Conformers_Atom[id][D.resNum]['OG1']
        if D.resName == 'TRP' and D.atomName == 'HE1':
            return self.Conformers_Atom[id][D.resNum]['NE1']
        if D.resName == 'TYR' and D.atomName == 'HH':
            return self.Conformers_Atom[id][D.resNum]['OH']
        
        return self.EMPTY

    #####
    # 获取环平面电流产生的化学位移shifts_ring
    #####
    def initOrbitalShift(self,):
        '''
        初始化环的各个原子,以及环平面
        '''
        conf = {}
        PF6Atoms = ["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]
        W6Atoms = ["CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3"]
        H5Atoms = ["CG", "ND1", "CE1", "NE2", "CD2"]
        W5Atoms = ["CG", "CD1", "NE1", "CE2", "CD2"]

        self.RingNo = 0
        for i in self.residueList[0].keys():
            resID = i
            resName = self.residueList[0][i]
            if resName=='PHE' or resName=='TYR' or resName=='TRP': ## 计算六元环
                self.Rings.append(RingData(resID=resID,resName=resName,atomNo=6,coordA=[[],[],[],[],[],[]]))
                if resName=='PHE':self.Rings[self.RingNo].ringFact = 1.46
                if resName=='TYR':self.Rings[self.RingNo].ringFact = 1.24
                if resName=='TRP':self.Rings[self.RingNo].ringFact = 1.24
                NO_MISSED_RING_ATOMS = True
                for j in range(6):
                    atom = PDB_Entry()
                    if resName=='PHE' or resName=='TYR':
                        atom = self.getEntry(0,resID,PF6Atoms[j])
                        if atom.atomName == PF6Atoms[j]:
                            self.Rings[self.RingNo].coordA[j] = atom.Coord.copy()
                        else:
                            ## missed ring atom in coordinates file
                            NO_MISSED_RING_ATOMS = False
                    elif resName=='TRP':
                        atom = self.getEntry(0,resID,W6Atoms[j])
                        if atom.atomName == W6Atoms[j]:
                            self.Rings[self.RingNo].coordA[j] = atom.Coord.copy()
                        else:
                            ## missed ring atom in coordinates file
                            NO_MISSED_RING_ATOMS = False
                if NO_MISSED_RING_ATOMS:self.RingNo+=1
            if resName=='HIS' or resName == 'TRP':   ## 计算五元环
                self.Rings.append(RingData(resID=resID,resName=resName,atomNo=5,coordA=[[],[],[],[],[]].copy()))
                if resName=='HIS':self.Rings[self.RingNo].ringFact = 1.35
                if resName=='TRP':self.Rings[self.RingNo].ringFact = 1.32
                NO_MISSED_RING_ATOMS = True
                for j in range(5):
                    atom = PDB_Entry()
                    if resName=='HIS':
                        atom = self.getEntry(0,resID,H5Atoms[j])
                        if atom.atomName == H5Atoms[j]:
                            self.Rings[self.RingNo].coordA[j] = atom.Coord.copy()
                        else:
                            NO_MISSED_RING_ATOMS = False
                    elif resName == 'TRP':
                        atom = self.getEntry(0,resID,W5Atoms[j])
                        if atom.atomName == W5Atoms[j]:
                            self.Rings[self.RingNo].coordA[j] = atom.Coord.copy()
                        else:
                            NO_MISSED_RING_ATOMS = False
                if NO_MISSED_RING_ATOMS:self.RingNo+=1
        for k in range(self.RingNo):
            self.calcPlane(self.Rings[k])

    def calcPlane(self,ringP:RingData):
        '''
        计算环平面的各项参数
        '''
        v1,v2 = [0,0,0],[0,0,0]
        self.Vec3Zero(ringP.center)
        for atomI in range(ringP.atomNo):
            self.Vec3Add(ringP.center,ringP.coordA[atomI])
        self.Vec3Scale(ringP.center,float(1.0/ringP.atomNo))   ##根据环周围的原子坐标，计算环中心点的坐标（均值）
        #if ringP.resID == 184:
        #    print(ringP.resName)
        self.Vec3Zero(ringP.norm)
        ringP.rad = 0.0
        v1 = ringP.coordA[ringP.atomNo-1].copy()
        self.Vec3Sub(v1,ringP.center)
        for atomI in range(ringP.atomNo):
            v2 = ringP.coordA[atomI]
            self.Vec3Sub(v2,ringP.center)
            ringP.rad += self.Vec3Abs(v2)
            self.Vec3Cross(v1,v2)
            self.Vec3Add(ringP.norm,v1)
            v1 = v2.copy()
        self.Vec3Norm(ringP.norm)
        ringP.rad /= ringP.atomNo
        
        # transM代表了什么参数?
        ringP.transM[0][2] = ringP.norm[0]
        ringP.transM[1][2] = ringP.norm[1]
        ringP.transM[2][2] = ringP.norm[2]

        ## 此处也并不知道在计算什么
        if ringP.norm[0] > 0.5 or ringP.norm[0] < -0.5:
            v1[0] = -ringP.norm[1]
            v1[1] = ringP.norm[0]
            v1[2] = 0.0
        else:
            v1[0] = 0.0
            v1[1] = -ringP.norm[2]
            v1[2] = ringP.norm[1]
        self.Vec3Norm(v1)
        ringP.transM[0][1] = v1[0]
        ringP.transM[1][1] = v1[1]
        ringP.transM[2][1] = v1[2]

        self.Vec3Cross(v1,ringP.norm)
        ringP.transM[0][0] = v1[0]
        ringP.transM[1][0] = v1[1]
        ringP.transM[2][0] = v1[2]

    def getOrbitalShift(self,conformerID:int,resNum:int,aName:str):
        '''
        获取ring_shift
        '''
        v = [0,0,0]
        v = self.getEntry(conformerID,resNum,aName).Coord.copy()
        ringShiftSum = 0.0
        for i in range(self.RingNo):
            vt,v1,v2 = [0,0,0],[0,0,0],[0,0,0]
            if self.Rings[i].resID == resNum and aName!='HN':
                continue
            atomNo = self.Rings[i].atomNo
            g = 0.0
            for atomI in range(atomNo):
                v1 = self.Rings[i].coordA[atomI].copy()
                v2 = self.Rings[i].coordA[(atomI+1) % atomNo].copy()
                r1 = self.Vec3DiffAbs(v1,v)
                r2 = self.Vec3DiffAbs(v2,v)
                vt = v.copy()

                self.Vec3Sub(vt,self.Rings[i].center)
                self.Vec3Sub(v1,self.Rings[i].center)
                self.Vec3Sub(v2,self.Rings[i].center)

                self.Mat3VecMult(vt,self.Rings[i].transM)
                self.Mat3VecMult(v1,self.Rings[i].transM)
                self.Mat3VecMult(v2,self.Rings[i].transM)

                self.Vec3Sub(v1,vt)
                self.Vec3Sub(v2,vt)

                area = v1[0]*v2[1]-v1[1]*v2[0]
                g += area*(1.0/(r1**3)+1.0/(r2**3))
            g *= 0.5
            b = float(5.4548e-6)
            ringShiftSum +=self.Rings[i].ringFact*b*g
        
        return ringShiftSum*(float(-1.0e6))

    #####
    # 计算HN_S2以及EF_shift
    #####
    def calc_HN_S2(self,):
        '''
        计算HN_S2
        '''
        conf = self.Conformers[0]

        for itA in self.residueList[0].keys():
            resID = itA
            resName = self.residueList[0][itA]

            try:
                O_prev = self.Conformers_Atom[0][resID-1]['O']
            except:O_prev = self.EMPTY
            try:
                H = self.Conformers_Atom[0][resID]['H']
            except:H = self.EMPTY
            try:
                HN = self.Conformers_Atom[0][resID]['HN']
            except:HN = self.EMPTY

            if O_prev.atomName == '':
                continue
            if H.atomName == '' and HN.atomName == '' and resName!='PRO': 
                self.HN_S2[resID] = 0
                continue    ##如果不是O/HN原子，就跳过
            if H.atomName == '': H = HN
            if resName == 'PRO': H = self.Conformers_Atom[0][resID]['N']

            S2 = 0
            for i in conf.keys():
                b = conf[i]
                resID2 = b.resNum
                if resID2==resID or resID2 == resID-1 or b.atomName[0] == 'H': continue

                D_OK = self.getDist(O_prev.Coord,b.Coord)
                D_HK = self.getDist(H.Coord,b.Coord)

                S2 += math.exp(-1.0*D_OK)
                S2 += 0.8 * math.exp(-1.0*D_HK)
            S2 *= 2.656
            S2 = (math.exp(S2)-math.exp(-S2))/(math.exp(S2)+math.exp(-S2)) - 0.1
            self.HN_S2[resID] = S2
        #if self.r1+1 not in self.HN_S2.keys():
        #    self.HN_S2[self.r1+1] = 0
        if self.HN_S2[self.r1+2]>0 and self.HN_S2[self.r1+1]>0:
            self.HN_S2[self.r1] = self.HN_S2[self.r1+1]-float(abs(self.HN_S2[self.r1+1]-self.HN_S2[self.r1+2]))
        if self.HN_S2[self.rN-2]>0 and self.HN_S2[self.rN-1]>0:
            self.HN_S2[self.rN] = self.HN_S2[self.rN-1]-float(abs(self.HN_S2[self.rN-1]-self.HN_S2[self.rN-2]))

    def calc_ElectricField(self,):
        '''
        计算Xef，电子场（力场）作用下的值
        什么值，具体定义还没有弄清楚
        '''
        conf = self.Conformers[0]

        targetAtomList = {'HN':'N','HA':'CA'}
        Qlist = {'C':-0.9612,'O':1.39374,'N':0.7209}

        #self.ElectricField = {}
        for itA in self.residueList[0].keys():
            resID = itA
            resName = self.residueList[0][itA]
            for itT in targetAtomList.keys():
                try:
                    target = self.Conformers_Atom[0][resID][itT]
                except:
                    target = self.EMPTY
                partner = self.Conformers_Atom[0][resID][targetAtomList[itT]]
                if target.atomName == '': continue
                for i in conf.keys():
                    b = conf[i]
                    resID2 = b.resNum
                    if float(abs(resID2-resID)) <= 1 : continue
                    atomName2 = b.atomName
                    if atomName2 == 'O' and target.atomName == 'HN': continue
                    if atomName2== 'O' or atomName2[0:2] == 'OD' or\
                        atomName2[0:2] == 'OE' or atomName2 == 'C' or atomName2 == 'N':
                        c = math.cos(self.getBondAngle_2(partner,target,b)*math.pi/180)
                        dist = self.getDist(target.Coord,b.Coord)
                        if dist>3.0 : continue
                        if resID not in self.ElectricField.keys():
                            self.ElectricField[resID] = {}
                        if target.atomName not in self.ElectricField[resID].keys():
                            self.ElectricField[resID][target.atomName] = 0
                        self.ElectricField[resID][target.atomName] += Qlist[atomName2[0:1]]*c/(dist*dist)

    #####
    # 氨基酸名字的转换
    #####
    def getThreeAAName(self,a:str):
        '''
        将氨基酸的单字符简称转化成为三字符简称
        '''
        if a == 'A': return "ALA"
        elif a == 'C': return "CYS"
        elif a == 'c': return "cys"
        elif a == 'D': return "ASP"
        elif a == 'E': return "GLU"
        elif a == 'F': return "PHE"
        elif a == 'G': return "GLY"
        elif a == 'H': return "HIS"
        elif a == 'I': return "ILE"
        elif a == 'K': return "LYS"
        elif a == 'L': return "LEU"
        elif a == 'M': return "MET"
        elif a == 'N': return "ASN"
        elif a == 'P': return "PRO"
        elif a == 'Q': return "GLN"
        elif a == 'R': return "ARG"
        elif a == 'S': return "SER"
        elif a == 'T': return "THR"
        elif a == 'V': return "VAL"
        elif a == 'W': return "TRP"
        elif a == 'Y': return "TYR"
        
        return "???"

    def getOneAAName(self,a:str):
        '''
        将氨基酸的三字符简称转化成为单字符简称
        '''
        if a == "ALA": return "A"
        elif a == "CYS": return "C"
        elif a == "cys": return "c"
        elif a == "ASP": return "D"
        elif a == "GLU": return "E"
        elif a == "PHE": return "F"
        elif a == "GLY": return "G"
        elif a == "HIS" or a == "HIS+": return "H"
        elif a == "ILE": return "I"
        elif a == "LYS" or a == "LYS+": return "K"
        elif a == "LEU": return "L"
        elif a == "MET": return "M"
        elif a == "ASN": return "N"
        elif a == "PRO": return "P"
        elif a == "GLN": return "Q"
        elif a == "ARG" or a == "ARG+": return "R"
        elif a == "SER": return "S"
        elif a == "THR": return "T"
        elif a == "VAL": return "V"
        elif a == "TRP": return "W"
        elif a == "TYR": return "Y"

        return "?"
    
    #####
    # 某一些运算符
    #####
    def sgn(self,x):
        '''
        正负指示器
        '''
        return ((x >= 0) - (x < 0))
    
    def Vec3Zero(self,v:list):
        '''
        向量归零
        '''
        for i in range(len(v)):
            v[i] = float(0.0)
    
    def Vec3Abs(self, v: list):
        '''
        取模
        '''
        return math.sqrt(v[0]**2+v[1]**2+v[2]**2)
    
    def Vec3DiffAbs(self,v1:list, v2:list):
        '''
        向量终点之间的距离
        '''
        return math.sqrt((v1[0]-v2[0])**2 + \
            (v1[1]-v2[1])**2 + \
            (v1[2]-v2[2])**2 )
 
    def Vec3Norm(self, v: list):
        '''
        归一化
        '''
        a = self.Vec3Abs(v)
        if a != 0:
            for i in range(len(v)):
                v[i] /= a

    def Vec3Scale(self,v:list,s:float):
        '''
        向量伸缩
        '''
        for i in range(len(v)):
            v[i] *= s

    def Vec3Add(self,v1:list,v2:list):
        '''
        向量加法
        '''
        for i in range(len(v1)):
            v1[i]+=v2[i]

    def Vec3Sub(self,v1:list,v2:list):
        '''
        向量减法
        v1 - v2
        '''
        for i in range(len(v1)):
            v1[i]-=v2[i]

    def Vec3Scalar(self,v1:list,v2:list):
        '''
        向量点乘（内积）
        '''
        return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

    def Vec3Cross(self,v1:list,v2:list):
        '''
        向量叉乘（外积）
        '''
        vRes = [0,0,0]
        vRes[0] = v1[1]*v2[2] - v1[2]*v2[1]
        vRes[1] = v1[2]*v2[0] - v1[0]*v2[2]
        vRes[2] = v1[0]*v2[1] - v1[1]*v2[0]
        v1 = vRes.copy()

    def Mat3VecMult(self,v:list,m: list):
        '''
        向量和矩阵(3x3)的外积
        '''
        vRes = [0,0,0]
        vRes[0] = v[0]*m[0][0] + v[1]*m[1][0] + v[2]*m[2][0]
        vRes[1] = v[0]*m[0][1] + v[1]*m[1][1] + v[2]*m[2][1]
        vRes[2] = v[0]*m[0][2] + v[1]*m[1][2] + v[2]*m[2][2]
        v = vRes.copy()

    def Vec3ScaleAdd(self,v1:list,s:float,v2:list):
        '''
        向量伸缩加法
        '''
        for i in range(len(v1)):
            v1[i] += s*v2[i]

    def isNum(self,a):
        try:
            t = float(a)
            return True
        except:
            return False

