import numpy as np
import math
import re
import time
from settings import RingData,PDB_Entry,MAXNUM,RAD
#
class PDB(object):
    '''
    '''

    def __init__(self,fileName:str,):
        '''
        '''
        self.PDBfileName = fileName
        self.loadPDB(self.PDBfileName)
        self.Rings = [RingData(),]  #芳香环的列表
        self.EMPTY = PDB_Entry()  #分子
        self.RingNo = 0
        self.Conformers = {}  # Mols, conformer ID + atom number 检索的原子列表
        self.ATOMS = {}     # conformer ID + residue number + atom name 检索的原子列表
        self.r1,self.rN = 9999,-9999
        self.resideList,self.resideListOne = [],[]   # 三氨基酸和单氨基酸的residue列表,

        self.acceptorList = {}   # 原子index检索的H键acceptor
        self.donorList = {}   # 原子标记或者关联的重原子检索的H键donor

        self.HBDistList = {}  # residue number+ atom name 检索的氢键距离
        self.HBEnergyList = {}  # residue number+ atom name 检索的氢键能量
        self.HB_DHO_AngleList = {}  # residue number+ atom name 检索的氢键夹角(Donator-H-O)
        self.HB_HOA_AngleList = {}   # residue number+ atom name 检索的氢键夹角(H-O-Acceptor_base)
        
        self.SpherePointNo = 0
        self.Star_SpherePoints = [0,0,0]
        self.SurfPrec = 3
        self.NeighborList = {}   # Neighboring atom list, indexed by atom number of target atom, and vector list of neighbor atoms
        self.ResSurfaceFullList,self.ResSurfacePartList = {},{}
        # Residue Surface list, indexed by residue number, and surface area for this residue
        self.AtomSurfaceFullList,self.AtomSurfacePartList = {},{}
        # Atom Surface list, indexed by residue number, atom name, and surface area for this residue
        self.VDW_RAD = {} # VDW radius for different type of atoms
        self.HN_S2 = {}
        self.ElectricField = {}  # indexed by resID and atomName
        

    #####
    # 加载PDB的坐标
    #####
    def getField(self,text:str,index:int):
        '''
        从PDB的原子field信息中获取对应的具体信息的方法
        '''
        if index==1:return text[0::6]  ##"ATOM"
        elif index==2:return text[6::5]  ##Atom number
        elif index==3:return text[11::5]  ##Atom name
        elif index==4:return text[17::4]  ##Residue name
        elif index==5:return text[21::1]  ##Chain ID
        elif index==6:return text[22::4]  ##Residue seq
        elif index==7:return text[30::8]  ##X
        elif index==8:return text[38::8]  ##Y
        elif index==9:return text[46::8]  ##Z
        elif index==10:return text[54::6]  ##Occupancy
        elif index==11:return text[60::6]  ##B-factor
        else:
            return None

    def loadPDB_Entry(self,text:str,entry:PDB_Entry):
        '''
        从PDB文件中提取信息
        '''
        entry.atomNum = int(self.getField(text,2))
        atomName = str(self.simplifyWhiteSpace(self.getField(text,3)))
        entry.atomName = atomName
        if atomName == 'H':
            entry.atomName = 'HN'
        if atomName[0]>='1' and atomName[0]<='3':
            entry.atomName = atomName[1::4] + atomName[0]
        entry.resName =  self.simplifyWhiteSpace(self.getField(str,4)[0::3])
        entry.chainName = self.getField(str,5)
        entry.resNum = int(self.getField(str,6))
        entry.X = float(self.getField(str,7))
        entry.Y = float(self.getField(str,8))
        entry.Z = float(self.getField(str,9))
        entry.B_Factor = float(self.getField(str,11))

        entry.Coord[0] = entry.X
        entry.Coord[1] = entry.Y
        entry.Coord[2] = entry.Z

    
    def loadPDB(self,fileName:str):
        '''
        '''



    #####
    # 获取原子
    #####
    def getEntry(self,conformerID:int,rNum:int,aName:str=''):
        '''
        '''
        if aName=='':
            return self.Conformers[conformerID][rNum]
        else:
            return self.ATOMS[conformerID][rNum][aName]


    def isSSBonded(self,conformerID:int,resNum:int):
        '''
        '''
        if self.resideListOne[resNum]!='C' and self.resideListOne[resNum] !='c':
            return False
        CYS_SG = self.getEntry(1,resNum,'SG')
        for i in range(len(self.resideListOne)):
            if (self.resideListOne[i]=='C' or self.resideListOne[i]=='c') and abs(i-resNum)>=4:
                if self.getDist(CYS_SG.Coord,self.ATOMS[1][i]['SG'].Coord)<=2.5:
                    return True
        return False
    
    #####
    # 获取角度和距离
    #####
    def getBondAngle(self, A: list, B: list, C: list):
        '''
        得到键键之间的夹角
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
        函数重载，
        '''
        return self.getBondAngle(a.Coord,b.Coord,c.Coord)
    
    def getDihedralAngle(self,a: PDB_Entry,b: PDB_Entry,c:PDB_Entry,d:PDB_Entry):
        '''
        求平面夹角
        '''
        cb,n1,n2 = [0,0,0],[0,0,0],[0,0,0]
        co = 0
        TEMP,TEMP1,TEMP2,TEMP3,TEMP4,TEMP5= 0,0,0,0,0
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
        TEMP,TEMP1,TEMP2,TEMP3,TEMP4,TEMP4,TEMP5 = \
            n1[0],n1[1],n1[2],n2[0],n2[1],n2[2]
        #计算两个平面的夹角，即abc平面与bcd平面的夹角的cos值
        #cos(a-b-c-d) = |n1·n2|/||*||
        co = float((n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2])/\
            math.sqrt((TEMP1**2+TEMP2**2+TEMP3**2)*(TEMP3**2+TEMP4**2+TEMP5**2)))
        return RAD * (math.acos(co) * \
            self.sgn(float((n1[1]*n2[2]-n1[2]*n2[1])*cb[0] + \
            (n1[2]*n2[0]-n1[0*n2[2]])*cb[1] + \
            (n1[0]*n2[1]-n1[1*n2[0]])*cb[2])))
        ##这个正负取得是什么值，两个法向量与BC垂直，则

    #####
    # 计算五个个不清楚的参数
    #####
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
        rName = self.resideList[resNum]
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
        rName = self.resideList[resNum]
        if rName.upper()=='GLY' or rName.upper()=='ALA':
            return MAXNUM

        if rName == 'VAL':
            R3 = self.getEntry(conformerID,resNum,'CG1')
            R4 = self.getEntry(conformerID,resNum,'HG11')
        elif rName == 'THR':
            R3 = self.getEntry(conformerID,resNum,'CG2')
            R4 = self.getEntry(conformerID,resNum,'HG21')
        elif rName == 'SER':
            R3 = self.getEntry(conformerID,resNum,'OG')
            R4 = self.getEntry(conformerID,resNum,'HG')
        elif rName == 'ILE':
            R3 = self.getEntry(conformerID,resNum,'CG1')
            R4 = self.getEntry(conformerID,resNum,'CD1')
        elif rName == 'CYS' or rName == 'cys':
            R3 = self.getEntry(conformerID,resNum,'SG')
            R4 = self.getEntry(conformerID,resNum,'HG')
        else:
            R3 = self.getEntry(conformerID,resNum,'CG')
        
        if rName=='ASN' or rName == 'ASP':
            R4 = self.getEntry(conformerID,resNum,'OD1')
        elif rName=='HIS':
            R4 = self.getEntry(conformerID,resNum,'ND1')
        elif rName=='MET':
            R4 = self.getEntry(conformerID,resNum,'SD')
        elif rName=='LEU' or rName=='PHE' or rName=='TYR' or rName=='TRP':
            R4 = self.getEntry(conformerID,resNum,'CD1')
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
    # 获取环当前的化学位移shifts
    #####
    def initOrbitalShift(self,):
        '''
        初始化环的各个原子及其位移
        '''
        conf = {}
        PF6Atoms = ["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]
        W6Atoms = ["CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3"]
        H5Atoms = ["CG", "ND1", "CE1", "NE2", "CD2"]
        W5Atoms = ["CG", "CD1", "NE1", "CE2", "CD2"]

        self.RingNo = 0
        for i in range(len(self.resideList)):
            resID = i
            resName = self.resideList[i]
            if resName=='PHE' or resName=='TYR' or resName=='TRP': ## 计算六元环
                self.Rings[self.RingNo].resNum = resID
                self.Rings[self.RingNo].resName = resName
                self.Rings[self.RingNo].atomNo = 6
                if resName=='PHE':self.Rings[self.RingNo].ringFact = 1.46
                if resName=='TYR':self.Rings[self.RingNo].ringFact = 1.24
                if resName=='TRP':self.Rings[self.RingNo].ringFact = 1.24
                NO_MISSED_RING_ATOMS = True
                for j in range(6):
                    atom = PDB_Entry()
                    if resName=='PHE' or resName=='TYR':
                        atom = self.getEntry(1,resID,PF6Atoms[i])
                        if atom.atomName == PF6Atoms[i]:
                            self.Rings[self.RingNo].coordA[i] = atom.Coord.copy()
                        else:
                            ## missed ring atom in coordinates file
                            NO_MISSED_RING_ATOMS = False
                    elif resName=='TRP':
                        atom = self.getEntry(1,resID,W6Atoms[i])
                        if atom.atomName == W6Atoms[i]:
                            self.Rings[self.RingNo].coordA[i] = atom.Coord.copy()
                        else:
                            ## missed ring atom in coordinates file
                            NO_MISSED_RING_ATOMS = False
                if NO_MISSED_RING_ATOMS:self.RingNo+=1
            if resName=='HIS' or resName == 'TRP':   ## 计算五元环
                self.Rings[self.RingNo].resID = resID
                self.Rings[self.RingNo].resName = resName
                self.Rings[self.RingNo].atomNo = 5
                if resName=='HIS':self.Rings[self.RingNo].ringFact = 1.35
                if resName=='TRP':self.Rings[self.RingNo].ringFact = 1.32
                NO_MISSED_RING_ATOMS = True
                for i in range(5):
                    atom = PDB_Entry()
                    if resName=='HIS':
                        atom = self.getEntry(1,resID,H5Atoms[i])
                        if atom.atomName == H5Atoms[i]:
                            self.Rings[self.RingNo].coordA[i] = atom.Coord.copy()
                        else:
                            NO_MISSED_RING_ATOMS = False
                    elif resName == 'TRP':
                        atom = self.getEntry(1,resID,W5Atoms[i])
                        if atom.atomName == W5Atoms[i]:
                            self.Rings[self.RingNo].coordA[i] = atom.Coord.copy()
                        else:
                            NO_MISSED_RING_ATOMS = False
                if NO_MISSED_RING_ATOMS:self.RingNo+=1
        for i in range(self.RingNo):
            self.calcPlane(self.Rings[i])

    def calcPlane(self,ringP:RingData):
        '''
        计算环平面的各项参数
        '''
        v1,v2 = [0,0,0],[0,0,0]
        self.Vec3Zero(ringP.center)
        for atomI in range(ringP.atomNo):
            self.Vec3Add(ringP.center,ringP.coordA[atomI])
        self.Vec3Scale(ringP.center,float(1.0/ringP.atomNo))   ##根据环周围的原子坐标，计算环中心点的坐标（均值）

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
    # 获取氢键的信息
    #####
    def initHBond(self,DIST:float,ANGLE:float):
        '''
        '''
        conf = self.Conformers[1]
        self.acceptorList = {}
        self.donorList = {}
        self.HBDistList = {}
        for i in range(len(conf)):
            a = self.isAcceptor(conf[i])
            if a.atomName!='':
                self.acceptorList[i] = a.atomNum
            else:
                temp = self.isDonor(conf[i])
                if temp.atomName=='': continue
                self.donorList[i] = temp.atomNum
        for itA in range(len(self.acceptorList)):
            A = conf[itA]
            A_Heavy = conf[self.acceptorList[itA]]
            for itD in range(len(self.donorList)):
                D = conf[itD]
                D_Heavy = conf[self.donorList[itD]]

                if abs(A.resNum - D.resNum) < 2: continue

                D_ON = self.getDist(A,D_Heavy) 
                D_CH = self.getDist(A_Heavy,D)
                D_OH = self.getDist(A,D)
                D_CN = self.getDist(A_Heavy,D_Heavy)

                if D_OH > D_CN : continue

                HBond_E = 332.0*0.42*0.20
                HBond_E *= (1.0/D_ON + 1.0/D_CH - 1.0/D_OH - 1.0/D_CN)

                if HBond_E < -0.5:
                    if self.HBDistList[A.resNum][A.resName] == 0 or \
                        self.HBDistList[A.resNum][A.atomName] > D_OH or \
                        self.HBEnergyList[A.resNum][A.atomName] > HBond_E:
                        self.HBDistList[A.resNum][A.atomName] = D_OH
                        self.HB_DHO_AngleList[A.resNum][A.atomName] = self.getBondAngle_2(D_Heavy,D,A)
                        self.HB_HOA_AngleList[A.resNum][A.atomName] = self.getBondAngle_2(D,A,A_Heavy)
                        self.HBEnergyList[A.resNum][A.atomName] = HBond_E
                    if self.HBDistList[D.resNum][D.atomName] == 0 or \
                        self.HBDistList[D.resNum][D.atomName] > D_OH or \
                        self.HBEnergyList[D.resNum][D.atomName] > HBond_E:
                        self.HBDistList[D.resNum][D.atomName] = D_OH
                        self.HB_DHO_AngleList[D.resNum][D.atomName] = self.getBondAngle_2(D_Heavy,D,A)
                        self.HB_HOA_AngleList[D.resNum][D.atomName] = self.getBondAngle_2(D,A,A_Heavy)
                        self.HBEnergyList[D.resNum][D.atomName] = HBond_E
        
    def getHBondDist(self,resNum:int,atomName:str):
        '''
        获取氢键的长度
        '''
        return self.HBDistList[resNum][atomName]

    def isAcceptor(self,A:PDB_Entry):
        '''

        '''
        if A.atomName=='O' or A.atomName=='OT1' or A.atomName=='OT2':
            return self.ATOMS[1][A.resNum]['C']
        if A.resName == 'ASN' and A.atomName == 'OD1':
            return self.ATOMS[1][A.resNum]['CG']
        if A.resName == 'ASP' and (A.atomName == 'OD1' or A.atomName == 'OD2'):
            return self.ATOMS[1][A.resNum]['CG']
        if (A.resName == 'CYS' or A.resName == 'cys') and A.atomName == 'SG':
            return self.ATOMS[1][A.resNum]['CB']
        if A.resName == 'GLN' and A.atomName == 'OE1':
            return self.ATOMS[1][A.resNum]['CD']
        if A.resName == 'GLU' and (A.atomName=='OE1' or A.atomName=='OE2'):
            return self.ATOMS[1][A.resNum]['CD']
        if A.resName == 'HIS' and (A.atomName == 'ND1' or A.atomName=='NE2'):
            return self.ATOMS[1][A.resNum]['CE1']
        if A.resName == 'MET' and A.atomName == 'SD':
            return self.ATOMS[1][A.resNum]['CG']
        if A.resName == 'SER' and A.atomName == 'OG':
            return self.ATOMS[1][A.resNum]['CB']
        if A.resName == 'THR' and A.atomName == 'OG1':
            return self.ATOMS[1][A.resNum]['CB']
        if A.resName == 'TYR' and A.atomName == 'OH':
            return self.ATOMS[1][A.resNum]['CZ']
        
        return self.EMPTY
 
    def isDonor(self,D:PDB_Entry):
        '''
        将H替换成H的donor？
        '''
        if D.atomName == 'HN':
            return self.ATOMS[1][D.resNum]['N']
        if re.search('HA',D.atomName) != None:
            return self.ATOMS[1][D.resNum]['CA']
        if D.resName == 'ARG':
            if D.atomName == 'HE':
                return self.ATOMS[1][D.resNum]['NE']
            if D.atomName == 'HH11' or D.atomName == 'HH12':
                return self.ATOMS[1][D.resNum]['NH1']
            if D.atomName == 'HH21' or D.atomName == 'HH22':
                return self.ATOMS[1][D.resNum]['NH2']
        if D.resName == 'ASP' and (D.atomName=='HD21' or D.atomName=='HD22'):
            return self.ATOMS[1][D.resNum]['ND2']
        if (D.resName == 'CYS' or D.resName == 'cys') and D.atomName == 'HG':
            return self.ATOMS[1][D.resNum]['SG']
        if D.resName == 'GLU' and (D.atomName=='HE21' or D.atomName=='HE22'):
            return self.ATOMS[1][D.resNum]['NE2']
        if D.resName == 'HIS':
            if D.atomName == 'HD1':
                return self.ATOMS[1][D.resNum]['ND1']
            if D.atomName == 'HE2':
                return self.ATOMS[1][D.resNum]['NE2']
        if D.resName == 'LYS' and (D.atomName=='HZ1' or D.atomName=='HZ2' or D.atomName=='HZ3'):
            return self.ATOMS[1][D.resNum]['NZ']
        if D.resName == 'SER' and D.atomName == 'HG':
            return self.ATOMS[1][D.resNum]['OG']
        if D.resName == 'THR' and D.atomName == 'HG1':
            return self.ATOMS[1][D.resNum]['OG1']
        if D.resName == 'TRP' and D.atomName == 'HE1':
            return self.ATOMS[1][D.resNum]['NE1']
        if D.resName == 'TYR' and D.atomName == 'HH':
            return self.ATOMS[1][D.resNum]['OH']
        
        return self.EMPTY

    #####
    # 某一些运算符
    #####
    def sgn(self,x):
        '''
        正负指示器
        '''
        return ((x >= 0) - (x < 0))

    """
    def arccos_(self,x:float):
        '''
        求arccos角度大小
        math.acos(x)
        '''
        if x>1: x=1
        elif x<-1: x=-1
        if abs(x)>=0.5:
            if x>0:
                return math.atan(math.sqrt(1.0-x*x)/x)
            else:
                return math.pi+math.atan(math.sqrt(1.0-x*x)/x)
        else:
            return 0.5*math.pi - math.atan(x/math.sqrt(1.0-x*x))
        
    """
    
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
    
    #####
    # 函数调用
    #####
    def simplifyWhiteSpace(self,text:str):
        '''
        简化空格，去掉头尾的空格以及中间的多余空格
        '''
        return re.sub(' +', ' ', text).strip()

        """
        if text=='':return text
        result = ''
        n = len(text)
        i,l = 0,0
        while True:
            while l<=n-1 and self.isSpace(text[l:l+1]):
                l+=1
            while l<=n-1 and not self.isSpace(text[l:l+1]):
                i+=1
                result+=text[l:l+1]
                l+=1
            if l<=n-1:
                i+=1
                result+=' '
            else:
                break
        if i>0 and result[i-1]==' ':
            i-=1
        result = result[:i+1]
        return result
        """

    def isSpace(self,c:str):
        '''
        判断字符非空
        '''
        if ord(c) >= 9 and ord(c) <= 13:
            return True
        return False
    
    def initSurface(self,rad_sol:float):
        '''
        初始化原子
        原子化半径
        化合物表面
        '''
        SolventRad = float(1.4)
        self.VDW_RAD['H'],self.VDW_RAD['HN'] = 1.2,1.2
        
        self.VDW_RAD['N'],self.VDW_RAD['ND1'],self.VDW_RAD['ND2'],self.VDW_RAD['ND2'] = 1.55,1.55,1.55,1.55,1.55
        self.VDW_RAD['NE'],self.VDW_RAD['NE1'],self.VDW_RAD['NE2'] = 1.55,1.55,1.55
        self.VDW_RAD['NH1'],self.VDW_RAD['NH2'],self.VDW_RAD['NZ'] = 1.55,1.55,1.55
        
        self.VDW_RAD['C'] = 2.3  ##1.1的C-H键长加上1.2的H原子半径
        
        self.VDW_RAD['CO'],self.VDW_RAD['CA'],self.VDW_RAD['CB'] = 1.7,1.7,1.7
        self.VDW_RAD['CG'],self.VDW_RAD['CG1'],self.VDW_RAD['CG2'] = 1.7,1.7,1.7
        self.VDW_RAD['CD'],self.VDW_RAD['CD1'],self.VDW_RAD['CD2'] = 1.7,1.7,1.7
        self.VDW_RAD['CE'],self.VDW_RAD['CE1'],self.VDW_RAD['CE2'] = 1.7,1.7,1.7
        self.VDW_RAD['CZ'],self.VDW_RAD['CZ2'],self.VDW_RAD['CZ3'],self.VDW_RAD['CH'] = 1.7,1.7,1.7,1.7

        self.VDW_RAD['O'],self.VDW_RAD['OD1'],self.VDW_RAD['OD2'] = 1.52,1.52,1.52
        self.VDW_RAD['OG'],self.VDW_RAD['OG1'] = 1.52,1.52
        self.VDW_RAD['OE1'],self.VDW_RAD['OE2'],self.VDW_RAD['OH'] = 1.52,1.52,1.52

        self.VDW_RAD['SG'],self.VDW_RAD['SD'] = 1.8,1.8

        self.SurfPrec = 3
        self.SpherePointNo = 0

        ## 计算表面位点 surface points
        start_T = time.time()
        self.SphereCalcPoints()
        finish_T = time.time()
        print('SphereCalcPoints running time: {} seconds.'.format(float(finish_T-start_T)))
        
        ## 计算邻近原子
        start_T = time.time()
        self.findNeighors(rad_sol)
        finish_T = time.time()
        print(''.format())

        ## 计算表面
        start_T = time.time()
        self.calcSurface(rad_sol)
        finish_T = time.time()
        print(''.format())

    def calcSurface(self,rad_sol:float):
        '''
        计算XX表面
        '''
        conf = self.Conformers[1]
        STD_AREA = {'ALA':(124,28),'CYS':(94,25),'HIS':(201,26),'MET':(215,28),'THR':(152,25),\
            'ARG':(244,27),'GLU':(187,23),'ILE':(194,23),'PHE':(221,26),'TRP':(265,28),\
            'ASN':(161,28),'GLN':(190,26),'LEU':(198,23),'PRO':(150,22),'TYR':(236,25),\
            'ASP':(154,25),'GLY':(89,30),'LYS':(214,26),'SER':(126,26),'VAL':(169,24),'cys':(94,25),}
        for i in range(len(conf)):
            b = conf[i]
            if b.atomName != 'HN' and b.atomName!='N' and b.atomName!='CA' \
                and b.atomName!='CB' and b.atomName!='C' and b.atomName!='O':
                continue
            #if b.atomName[0] == 'H': continue   ## for all heavy atoms
            #if b.atomName != 'O': continue     ## for Carbonyl O only

            resID = b.resNum
            Neighnors = self.NeighborList[b.atomNum]
            fullNo = 0
            partNo = 0
            cent,nCent,x,dx = [0,0,0],[0,0,0,],[0,0,0],[0,0,0]

            b_atomName = b.atomName
            if b.atomName[0]=='C' and b.atomName!= 'C': b_atomName = 'C'  ## for all non-C C
            if b.atomName == 'C':b_atomName='CO' ## for all C=O C
            b_rad = float(self.VDW_RAD[b_atomName]+rad_sol)
            cent = b.Coord.copy()
            for pointI in range(self.SpherePointNo):
                x = cent.copy()
                self.Vec3ScaleAdd(x,b_rad,self.Star_SpherePoints[pointI])
                fullInside,partInside = False,False
                for atomI in range(len(Neighnors)):
                    nAtomP = conf[Neighnors[atomI]]
                    nCent = nAtomP.Coord.copy()
                    dx = x.copy()
                    self.Vec3Sub(dx,nCent)
                    
                    if nAtomP.resNum < self.r1 or nAtomP>self.rN: continue
                    nAtomP_atomName = nAtomP.atomName
                    if nAtomP.atomName[0]=='C':nAtomP_atomName='C'
                    if nAtomP.atomName == 'C': nAtomP_atomName='CO'
                    nRad = self.VDW_RAD[nAtomP_atomName]+rad_sol
                    if dx[0]**2+dx[1]**2+dx[2]**2 > nRad**2:
                        continue
                    fullInside = True
                    if resID == nAtomP.resNum:
                        partInside = True
                        break

                if not fullInside:
                    fullNo+=1
                if not partInside:
                    partNo+=1
            self.AtomSurfaceFullList[resID][b.atomName] = float(fullNo*b_rad*b_rad*4.0*math.pi/self.SpherePointNo)
            self.AtomSurfacePartList[resID][b.atomName] = float(partNo*b_rad*b_rad*4.0*math.pi/self.SpherePointNo)

            self.ResSurfaceFullList[resID] += float(fullNo*b_rad*b_rad*4.0*math.pi/self.SpherePointNo)
            self.ResSurfacePartList[resID] += float(partNo*b_rad*b_rad*4.0*math.pi/self.SpherePointNo)
        
        for itS in range(len(self.ResSurfaceFullList)):
            resName = self.resideList[itS]

    def findNeighors(self,rad_sol:float):
        '''
        find the neighoring atoms, distance < rad_a + rad_b + rad_sol\n
        for atoms with given type, for all non-H atoms
        '''
        conf = self.Conformers[1]
        for itA in range(len(conf)):
            a = conf[itA]
            if a.atomName != 'HN' and a.atomName!='N' and a.atomName!='CA' \
                and a.atomName!='CB' and a.atomName!='C' and a.atomName!='O':
                continue
            a_atomName = a.atomName
            if a.atomName[0] =='C': a_atomName = 'C'
            if a.atomName == 'C': a_atomName = 'CO'
            rad_a = self.VDW_RAD[a_atomName]
            for i in range(len(conf)):
                b = conf[i]
                if b.atomName[0]=='H' and b.atomName!='HN' and b.atomName!='H':
                    continue
                b_atomName = b.atomName
                if b.atomName[0] == 'C':b_atomName='C'
                if b.atomName == 'C':b_atomName='CO'
                rad_b = self.VDW_RAD[b_atomName]
                if self.getDist(a,b)<rad_a+rad_b+rad_sol:
                    self.NeighborList[a.atomNum].append(b.atomNum)

    def calcTriangles(self,\
        x0:float,y0:float,z0:float,\
        x1:float,y1:float,z1:float,\
        x2:float,y2:float,z2:float,\
        rowStartA:list,rowNo:int,quad:int,\
        row0:int,ind0:int,ind1:int,\
        row2:int,ind2:int,\
        pointA:list):
        '''
        计算三原子夹角
        这是个迭代函数，但是我并没有搞清楚它的迭代运算过程
        '''
        if row0 + 1 == row2 or row2 + 1 == row0:
            row0Size,row2Size = 0,0
            if row0 == -rowNo or row0 == rowNo:
                row0Size = 1
            elif row0<0:
                row0Size = rowNo + row0
            else:
                row0Size = rowNo - row0

            if row2 == -rowNo or row2 == rowNo:
                row2Size = 1
            elif row2 < 0 :
                row2Size = rowNo + row2
            else:
                row2Size = rowNo - row2
            
            if ind0 < (quad+1)*row0Size:
                ind0 += rowStartA[rowNo+row0]
                pointA[ind0][0] = float(x0)
                pointA[ind0][1] = float(y0)
                pointA[ind0][2] = float(z0)

            if ind1 < (quad+1)*row0Size:
                ind1 += rowStartA[rowNo+row0]
                pointA[ind1][0] = float(x1)
                pointA[ind1][1] = float(y1)
                pointA[ind1][2] = float(z1)

            if ind2 < (quad+1)*row0Size:
                ind2 += rowStartA[rowNo+row0]
                pointA[ind2][0] = float(x2)
                pointA[ind2][1] = float(y2)
                pointA[ind2][2] = float(z2)
        else:
            Matrix_xyz = [[0,0,0],[0,0,0],[0,0,0]]
            Matrix_xyz[0][0] = x0 + x1
            Matrix_xyz[0][1] = y0 + y1
            Matrix_xyz[0][2] = z0 + z1
            a = math.sqrt(sum(i**2 for i in Matrix_xyz[0]))
            for iter_I in range(len(Matrix_xyz[0])):
                Matrix_xyz[0][iter_I] /= a
            Matrix_xyz[1][0] = x1 + x2
            Matrix_xyz[1][1] = y1 + y2
            Matrix_xyz[1][2] = z1 + z2
            a = math.sqrt(sum(i**2 for i in Matrix_xyz[1]))
            for iter_I in range(len(Matrix_xyz[1])):
                Matrix_xyz[1][iter_I] /= a
            Matrix_xyz[2][0] = x2 + x0
            Matrix_xyz[2][1] = y2 + y0
            Matrix_xyz[2][2] = z2 + z0
            a = math.sqrt(sum(i**2 for i in Matrix_xyz[2]))
            for iter_I in range(len(Matrix_xyz[2])):
                Matrix_xyz[2][iter_I] /= a
            
            rowMid = (row0+row2)/2
            indMid01 = (ind0+ind1)/2
            indMid12 = (ind1+ind2)/2
            indMid02 = (ind0+ind2)/2

            self.calcTriangles(\
                x0,y0,z0,\
                Matrix_xyz[0][0],Matrix_xyz[0][1],Matrix_xyz[0][2],\
                Matrix_xyz[2][0],Matrix_xyz[2][1],Matrix_xyz[2][2],\
                rowStartA,rowNo,quad,\
                row0,ind0,indMid01,\
                rowMid,indMid02,\
                pointA)
            self.calcTriangles(\
                Matrix_xyz[0][0],Matrix_xyz[0][1],Matrix_xyz[0][2],\
                x1,y1,z1,\
                Matrix_xyz[1][0],Matrix_xyz[1][1],Matrix_xyz[1][2],\
                rowStartA,rowNo,quad,\
                row0,indMid01,ind1,\
                rowMid,indMid12,\
                pointA)
            self.calcTriangles(\
                Matrix_xyz[2][0],Matrix_xyz[2][1],Matrix_xyz[2][2],\
                Matrix_xyz[1][0],Matrix_xyz[1][1],Matrix_xyz[1][2],\
                Matrix_xyz[0][0],Matrix_xyz[0][1],Matrix_xyz[0][2],\
                rowStartA,rowNo,quad,\
                rowMid,indMid02,indMid12,\
                row0,indMid01,\
                pointA)
            self.calcTriangles(\
                Matrix_xyz[2][0],Matrix_xyz[2][1],Matrix_xyz[2][2],\
                Matrix_xyz[1][0],Matrix_xyz[1][1],Matrix_xyz[1][2],\
                x2,y2,z2,\
                rowStartA,rowNo,quad,\
                rowMid,indMid02,indMid12,\
                row2,ind2,\
                pointA)
        
    def SphereCalcPoints(self,):
        '''
        计算表面原子位点
        得到Star_SpherePoints和SpherePointNo
        '''
        pointA = [0,0,0]
        #rowStartA = 4  ## int *rowStartA
        pointNo,rowNo,rowSize = 0,0,0

        rowNo = 1 << self.SurfPrec
        l = int((2*rowNo+1)*4)
        rowStartA = [0]*l
        rowStartA[0] = 0
        rowStartA[1] = 1
        rowSize = 4

        for i in range(2,rowNo+1):
            rowStartA[i] = rowStartA[i-1]+rowSize
            rowSize+=4
        for i in range(rowNo+1,2*rowNo+1):
            rowStartA[i] = rowStartA[i-1]+rowSize
            rowSize-=4
        
        pointNo = 4*rowNo*rowNo+2
        l = pointNo*4*len(pointA)   ###### 注意sizeof(Vec3)
        pointA = [[0,0,0]]*l
        self.calcTriangles(\
            1.0, 0.0, 0.0,\
		    0.0, 1.0, 0.0,\
		    0.0, 0.0, 1.0,\
		    rowStartA, rowNo, 0,\
		    0, 0, rowNo,\
		    rowNo, 0,\
		    pointA)
        self.calcTriangles(
            0.0, 1.0, 0.0,
            -1.0, 0.0, 0.0,
            0.0, 0.0, 1.0,
            rowStartA, rowNo, 1,
            0, rowNo, 2 * rowNo,
            rowNo, 0,
            pointA)
        self.calcTriangles(
            -1.0, 0.0, 0.0,
            0.0, -1.0, 0.0,
            0.0, 0.0, 1.0,
            rowStartA, rowNo, 2,
            0, 2 * rowNo, 3 * rowNo,
            rowNo, 0,
            pointA)
        self.calcTriangles(
            0.0, -1.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 0.0, 1.0,
            rowStartA, rowNo, 3,
            0, 3 * rowNo, 4 * rowNo,
            rowNo, 0,
            pointA)
        self.calcTriangles(
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, -1.0,
            rowStartA, rowNo, 0,
            0, 0, rowNo,
            - rowNo, 0,
            pointA)
        self.calcTriangles(
            0.0, 1.0, 0.0,
            -1.0, 0.0, 0.0,
            0.0, 0.0, -1.0,
            rowStartA, rowNo, 1,
            0, rowNo, 2 * rowNo,
            - rowNo, 0,
            pointA)
        self.calcTriangles(
            -1.0, 0.0, 0.0,
            0.0, -1.0, 0.0,
            0.0, 0.0, -1.0,
            rowStartA, rowNo, 2,
            0, 2 * rowNo, 3 * rowNo,
            - rowNo, 0,
            pointA)
        self.calcTriangles(
            0.0, -1.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 0.0, -1.0,
            rowStartA, rowNo, 3,
            0, 3 * rowNo, 4 * rowNo,
            - rowNo, 0,
            pointA)
        
        del rowStartA[:]
        self.Star_SpherePoints = pointA
        self.SpherePointNo = pointNo

    def calc_HN_S2(self,):
        '''
        计算HN_S2
        '''
        conf = self.Conformers[1]

        for itA in range(len(self.resideList)):
            resID = itA
            resName = self.resideList[itA]

            O_prev = self.ATOMS[1][resID-1]['O']
            H = self.ATOMS[1][resID]['H']
            HN = self.ATOMS[1][resID]['HN']

            if O_prev.atomName == '' or (H.atomName == '' and HN.atomName == '' and resName!='PRO'): continue    ##如果不是O/HN原子，就跳过
            if H.atomName == '': H = HN
            if resName == 'PRO': H = self.ATOMS[1][resID]['N']

            S2 = 0
            for i in range(len(conf)):
                b = conf[i]
                resID2 = b.resNum
                if resID2==resID or resID2 == resID-1 or b.atomName[0] == 'H' or (b.atomName[1]=='H' and b.atomName[1] != 'N'): continue

                D_OK = self.getDist(O_prev,b)
                D_HK = self.getDist(H,b)

                S2 += math.exp(-1.0*D_OK)
                S2 += 0.8 * math.exp(-1.0*D_HK)
            S2 *= 2.656
            S2 = (math.exp(S2)-math.exp(-S2))/(math.exp(S2)+math.exp(-S2)) - 0.1
            self.HN_S2[resID] = S2
        if self.HN_S2[self.r1+2]>0 and self.HN_S2[self.r1+1]>0:
            self.HN_S2[self.r1] = self.HN_S2[self.r1+1]-float(abs(self.HN_S2[self.r1+1]-self.HN_S2[self.r1+2]))
        if self.HN_S2[self.rN-2]>0 and self.HN_S2[self.rN-1]>0:
            self.HN_S2[self.rN] = self.HN_S2[self.rN-1]-float(abs(self.HN_S2[self.rN-1]-self.HN_S2[self.rN-2]))

    def calc_ElectricField(self,):
        '''
        计算Xef，电子场（力场）作用下的值
        什么值，具体定义还没有弄清楚
        '''
        conf = self.Conformers[1]
        targetAtomList = {'HN':'N','HA':'CA'}
        Qlist = {'C':-0.9612,'O':1.39374,'N':0.7209}
        self.ElectricField = {}
        for itA in range(len(self.resideList)):
            resID = itA
            resID = self.resideList[itA]
            for itT in range(len(targetAtomList)):
                target = self.ATOMS[1][resID][itT]
                partner = self.ATOMS[1][resID][targetAtomList[itT]]
                if target.atomName == '': continue
                for i in range(len(conf)):
                    b = conf[i]
                    resID2 = b.resNum
                    if float(abs(resID2-resID)) <= 1 : continue
                    atomName2 = b.atomName
                    if atomName2 == 'O' and target.atomName == 'HN': continue
                    if atomName2== 'O' or atomName2[0::2] == 'OD' or\
                        atomName2[0::2] == 'OE' or atomName2 == 'C' or atomName2 == 'N':
                        c = math.cos(self.getBondAngle_2(partner,target,b)*math.pi/180)
                        dist = self.getDist(target,b)
                        if dist>3.0 : continue
                        self.ElectricField[resID][target.atomName] += Qlist[atomName2[0::1]]*c/(dist*dist)

    def collect_HN_S2_and_EF(self,):
        '''
        结合上述两种计算
        '''
        conf = self.Conformers[1]
        self.ElectricField = {}
        targetAtomList = {}
        targetAtomList['HN'] = 'N'
        targetAtomList['HA'] = 'CA'

        Qlist = {}
        Qlist['C'] = -0.9612
        Qlist['O'] = 1.39374
        Qlist['N'] = 0.7209

        for itA in range(len(self.resideList)):
            resID = itA
            resName = self.resideList[itA]

            O_prev = self.ATOMS[1][resID-1]['O']
            H = self.ATOMS[1][resID]['H']
            HN = self.ATOMS[1][resID]['HN']

            if O_prev.atomName == '' or (H.atomName == '' and\
                HN.atomName == '' and resName != 'PRO'):
                continue
            if H.atomName == '': H = HN
            if resName == 'PRO': H = self.ATOMS[1][resID]['N']

            EF_target_HA = self.ATOMS[1][resID]['HA']
            EF_partner_HA = self.ATOMS[1][resID]['CA']
            EF_target_HN = self.ATOMS[1][resID]['HN']
            EF_partner_HN = self.ATOMS[1][resID]['N']

            S2 = 0
            for i in range(len(conf)):
                b = conf[i]
                resID2 = b.resNum
                if not (resID2==resID or resID2 == resID - 1 or b.atomName[0]=='H'\
                    or (b.atomName[1] == 'H' and b.atomName[1] != 'N')):
                    D_OK = self.getDist(O_prev,b)
                    D_HK = self.getDist(H,b)
                    S2 += math.exp(-1.0*D_OK)
                    S2 += 0.8*math.exp(-1.0*D_HK)
                
                if float(abs(resID2-resID))<=1: continue
                atomName2 = b.atomName

                if EF_target_HA.atomName != '':
                    if atomName2 == 'O' or atomName2[0::2] == 'OD' or atomName2[0::2] == 'OE'\
                        or atomName2=='C' or atomName2 == 'N':
                        dist = self.getDist(EF_target_HA,b)
                        if dist<=3.0:
                            self.ElectricField[resID]['HA'] += Qlist[atomName2[0::1]]*math.cos(\
                                self.getBondAngle_2(EF_partner_HA,EF_target_HA,b)*math.pi/180)/(dist*dist)
                if EF_target_HN.atomName != '':
                    if atomName2[0::2] == 'OD' or atomName2[0::2] == 'OE'\
                        or atomName2=='C' or atomName2 == 'N':
                        dist = self.getDist(EF_target_HN,b)
                        if dist<=3.0:
                            self.ElectricField[resID]['HN'] += Qlist[atomName2[0::1]]*math.cos(\
                                self.getBondAngle_2(EF_partner_HN,EF_target_HN,b)*math.pi/180)/(dist*dist)
            S2*=2.656
            S2 = (math.exp(S2)-math.exp(-S2))/(math.exp(S2)+math.exp(-S2))-0.1
            self.HN_S2[resID] = S2
        if self.HN_S2[self.r1+2]>0 and self.HN_S2[self.r1+1]>0:
            self.HN_S2[self.r1] = self.HN_S2[self.r1+1]-float(abs(self.HN_S2[self.r1+1]-self.HN_S2[self.r1+2]))
        if self.HN_S2[self.rN-2]>0 and self.HN_S2[self.rN-1]>0:
            self.HN_S2[self.rN] = self.HN_S2[self.rN-1]-float(abs(self.HN_S2[self.rN-1]-self.HN_S2[self.rN-2]))

    
    
    









