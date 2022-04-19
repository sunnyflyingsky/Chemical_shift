import numpy as np
import math

class ANN(object):
    '''
    '''
    
    def __init__(self,dPATH:str,dNAME_PREFIX:str,N1_nodeI:int=96,N1_nodeH:int=20,N1_nodeO:int=3,\
        N2_nodeI:int=9,N2_nodeH:int=6,N2_nodeO:int=3):
        '''
        '''
        self.buf = []
        self.slash_char = ''
        #####
        #code for input format
		# 0 - all CS and AA data, default
		# 1 - no CS for residue i-1
		# 2 - no CS for residue i
		# 3 - no CS for residue i+1
        #####
        self.input_code = 0  #0-3
        self.DB_PATH,self.DB_NAME_PREFIX = dPATH,dNAME_PREFIX
        self.DB_FNAME_LEVEL = [['','',''],['','','']]
        self.N1_NODE = [N1_nodeI,N1_nodeH,N1_nodeO]
        self.N2_NODE = [N2_nodeI,N2_nodeH,N2_nodeO]

        self.WI_1, self.WI_2, self.WI_3 = [],[],[]
        self.BI_1, self.BI_2, self.BI_3 = [],[],[]
        self.WL1_1, self.WL1_2, self.WL1_3 = [],[],[]
        self.BL1_1, self.BL1_2, self.BL1_3 = [],[],[]
        self.WL2_1, self.WL2_2, self.WL2_3 = [],[],[]
        self.BL2_1, self.BL2_2, self.BL2_3 = [],[],[]

        self.W2I_1, self.W2I_2, self.W2I_3 = [],[],[]
        self.B2I_1, self.B2I_2, self.B2I_3 = [],[],[]
        self.W2L1_1, self.W2L1_2, self.W2L1_3 = [],[],[]
        self.B2L1_1, self.B2L1_2, self.B2L1_3 = [],[],[]
        self.W2L2_1, self.W2L2_2, self.W2L2_3 = [],[],[]
        self.B2L2_1, self.B2L2_2, self.B2L2_3 = [],[],[]

        self.ANN_IN_MTX = []
        self.ANN_IN_MTX_LEVEL1 = []
        self.ANN_IN_MTX_LEVEL2 = []
        self.ANN_OUT_MTX_LEVEL1 = []
        self.it = {}
        self.ANN_OUT_MTX_LEVEL2 = []



    
    def set_input_code(self,c:int):
        '''
        '''
        self.input_code=c
    
    def loadWeights(self,):
        '''
        '''
        wName = ''
        wName = self.DB_PATH+self.slash_char+self.DB_NAME_PREFIX+'.level1.WI.tab'
        self.loadWeightBias3(wName,self.WI_1,self.BI_1,self.WI_2,self.BI_2,self.WI_3,self.BI_3,self.N1_NODE[0],self.N1_NODE[0],self.N1_NODE[0])

        wName = self.DB_PATH+self.slash_char+self.DB_NAME_PREFIX+'.level1.WL1.tab'
        self.loadWeightBias3(wName,self.WL1_1,self.BL1_1,self.WL1_2,self.BL1_2,self.WL1_3,self.BL1_3,self.N1_NODE[1],self.N1_NODE[0],self.N1_NODE[1])

        wName = self.DB_PATH+self.slash_char+self.DB_NAME_PREFIX+'.level1.WL2.tab'
        self.loadWeightBias3(wName,self.WL2_1,self.BL2_1,self.WL2_2,self.BL2_2,self.WL2_3,self.BL2_3,self.N1_NODE[2],self.N1_NODE[1],self.N1_NODE[2])

    def loadWeightBias3(self,fName:str,W1:list,B1:list,W2:list,B2:list,W3:list,B3:list,N_W_row:int,N_W_col:int,N_B:int):
        '''
        加载N_W_row*N_W_col大小的weighting矩阵以及bias(N_B)
        '''
        
    
    #####
    #执行两层ANN计算
    #####
    def runANN(self,inMatrix:list):
        '''
        '''
        pass

    def calcLevel1(self,):
        '''
        计算第一层
        '''
        for i in range(len(self.ANN_IN_MTX_LEVEL1)):
            IL1,IL2,IL3 = [],[],[]
            self.applyANNTransformation(self.ANN_IN_MTX_LEVEL1[i],self.WI_1,self.BI_1,IL1,1)
            self.applyANNTransformation(self.ANN_IN_MTX_LEVEL1[i],self.WI_2,self.BI_2,IL2,1)
            self.applyANNTransformation(self.ANN_IN_MTX_LEVEL1[i],self.WI_3,self.BI_3,IL3,1)

            HL1,HL2,HL3 = [],[],[]
            self.applyANNTransformation(IL1,self.WL1_1,self.BL1_1,HL1,1)
            self.applyANNTransformation(IL2,self.WL1_2,self.BL1_2,HL2,1)
            self.applyANNTransformation(IL3,self.WL1_3,self.BL1_3,HL3,1)

            OL1,OL2,OL3 = [],[],[]
            self.applyANNTransformation(HL1,self.WL2_1,self.BL2_1,OL1,0)
            self.applyANNTransformation(HL2,self.WL2_2,self.BL2_2,OL2,0)
            self.applyANNTransformation(HL3,self.WL2_3,self.BL2_3,OL3,0)

            OUT1 = []
            self.applyVecAverage(OL1,OL2,OL3,OUT1)
            self.ANN_OUT_MTX_LEVEL1[i] = OUT1

    def calcLevel2(self,):
        '''
        计算第二层(隐藏层)
        '''
        for i in range(len(self.ANN_IN_MTX_LEVEL2)):
            IL1,IL2,IL3 = [],[],[]
            self.applyANNTransformation(self.ANN_IN_MTX_LEVEL2[i],self.W2I_1,self.B2I_1,IL1,1)
            self.applyANNTransformation(self.ANN_IN_MTX_LEVEL2[i],self.W2I_2,self.B2I_2,IL2,1)
            self.applyANNTransformation(self.ANN_IN_MTX_LEVEL2[i],self.W2I_3,self.B2I_3,IL3,1)

            HL1,HL2,HL3 = [],[],[]
            self.applyANNTransformation(IL1,self.W2L1_1,self.B2L1_1,HL1,1)
            self.applyANNTransformation(IL2,self.W2L1_2,self.B2L1_2,HL2,1)
            self.applyANNTransformation(IL3,self.W2L1_3,self.B2L1_3,HL3,1)

            OL1,OL2,OL3 = [],[],[]
            self.applyANNTransformation(HL1,self.W2L2_1,self.B2L2_1,OL1,0)
            self.applyANNTransformation(HL2,self.W2L2_2,self.B2L2_2,OL2,0)
            self.applyANNTransformation(HL3,self.W2L2_3,self.B2L2_3,OL3,0)

            OUT2 = []
            self.applyVecAverage(OL1,OL2,OL3,OUT2)
            self.ANN_OUT_MTX_LEVEL2[i] = OUT2

    def runSpartaANN(self,inMatrix:list):
        '''
        开始执行ANN
        '''
        self.ANN_IN_MTX_LEVEL1 = inMatrix
        self.calcLevel1() 

    
    def applyANNTransformation(self,inp:list,w:list,b:list,out:list,code:int):
        '''
        跃迁函数以及层数之间的切换
        '''
        if len(inp) != len(w[0]) or len(w)!=len(b):
            print(' ANN prediction failed with inconsistent data!')
            exit(0)
        for i in range(len(w)):
            sum = 0
            for j in range(len(inp)):
                sum+=inp[j]*w[i][j]
            sum+=b[i]
            if code==1:
                out.append(2.0/(1.0+math.exp(-2.0*sum))-1.0)
            elif code == 0:
                out.append(sum)
    
    def applyVecAverage(self,v1:list,v2:list,v3:list,vout:list):
        '''
        向量加和均值
        '''
        if len(v1)==len(v2) and len(v1)==len(v3):
            for i in range(len(v1)):
                vout.append((v1[i]+v2[i]+v3[i])/3.0)
    
    def applyVecNormalization(self,v:list):
        '''
        向量标准化处理
        '''
        a,b,c = v[0],v[1],v[2]
        if a>1:a=1.0 
        elif a<0:a=0.0
        if b>1:b=1.0 
        elif b<0:b=0.0
        if c>1:c=1.0 
        elif c<0:c=0.0
        sum_ = a+b+c
        a /= sum_
        b /= sum_
        c /= sum_
        v[0],v[1],v[2] = a,b,c

    def getConfidence(self,v:list):
        '''
        max - mid
        '''
        if len(v)!=3:return -1.0
        return 2.0*max(v) - sum(v) + min(v)
    
    def getNumberMissCS(self,v:list):
        '''
        检查给定残基的缺失原子数
        '''
        cnt = 0
        for i in range(1,13,2):
            cnt+=(v[int]==1)
        return cnt

    
    def itoa(self,n:int,buff:str,base:int):
        '''
        输出整数型数的字符串形式？
        '''
        print('{:%d}{}'.format(buff,n))
        return buff
    
    def ftoa(self,n:float,buff:str,f:str,prec:int):
        '''
        输出浮点型数的字符串形式？
        '''
        if not (f.lower()=='f' or f.lower()=='e' or f.lower()=='g') :
            f = 'f'
        formatx = ''
        fs = 0
        print(buff,formatx,n)
        return buff
    
    def getSlashChar(self,):
        '''
        '''








