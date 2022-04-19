MAXNUM = 9999
RAD = 57.29577951
###
# Vec3, Mat3,DoubleD2,DoubleD1,
# Mol,Mols,
###
class PDB_Entry(object):
    def __init__(self,atomNum:int=0,resNum:int=0,atomName:str='',resName:str='',\
        chainName:str='',X:float=0,Y:float=0,Z:float=0,B_Factor:float='',Coord:list=[0,0,0].copy()):
        self.atomNum = atomNum
        self.resNum = resNum
        self.atomName = atomName
        self.resName = resName
        self.chainName = chainName
        self.X = X
        self.Y = Y
        self.Z = Z
        self.B_Factor = B_Factor
        self.Coord = Coord.copy()


class RingData(object):
    def __init__(self,resID=0,atomNo=0,resName='',atoms='',coordA=[].copy(),\
        center=[0,0,0].copy(),norm=[0,0,0].copy(),transM=[[0,0,0].copy(),[0,0,0].copy(),[0,0,0].copy()].copy(),rad=0,ringFact=0):
        self.resID = resID
        self.atomNo = atomNo
        self.resName = resName
        self.atoms = atoms
        self.coordA = coordA.copy()   ##原子坐标
        self.center = center.copy()  ##环中心
        self.norm = norm.copy()  ##法向量
        self.transM = [].copy()
        self.transM.append(transM[0].copy())
        self.transM.append(transM[1].copy())
        self.transM.append(transM[2].copy())
        self.rad = rad   ##半径
        self.ringFact = ringFact

