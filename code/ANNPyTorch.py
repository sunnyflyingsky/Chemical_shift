import numpy as np
import torch
from torch import nn
import pandas as pd
from sklearn.model_selection import train_test_split
from torch.autograd import Variable
import matplotlib.pyplot as plt 
import seaborn as sns
import sys

def main(atomName:str):
    '''
    训练模型
    '''
    batch_size = 256
    #数据处理
    data = pd.read_csv('Working/Sparta+/data/train_set_'+atomName+'.csv',sep=',',header=0)
    data = data.loc[:,data.columns!='name']
    data = data.astype(np.float32)

    for key in data.columns:
        data.fillna({key:data[key].median()},inplace=True)
    shift_name = 'shift2'
    data = data_Normalization(data,shift_name)
    y = data.loc[:,data.columns==shift_name].values
    
    data = data.loc[:,data.columns!='shift1']
    data = (data.loc[:,data.columns!=shift_name]-data.loc[:,data.columns!=shift_name].min())/\
        (data.loc[:,data.columns!=shift_name].max()-data.loc[:,data.columns!=shift_name].min())
    #如果你要去除某一些特征，请在此处修改，并且修改神经网络的输入层神经元个数为新的个数
    #drop_list = ['HB_1O_Energy','HB_2HN_Energy','HB_2HA_Energy','HB_2O_Energy','HB_3HN_Energy']
    #for feature_name in drop_list:
    #    data = data.loc[:,data.columns!=feature_name]
    x = data.loc[:,data.columns!=shift_name].values
    #data = (data-data.mean())/(data.std())
    
    #数据集划分
    x_train,x_test,y_train,y_test = train_test_split(x,y,test_size=0.2,random_state=2021)

    #数据格式转化
    x_train = torch.from_numpy(x_train)
    y_train = torch.from_numpy(y_train)#.type(torch.LongTensor)
    x_test = torch.from_numpy(x_test)
    y_test = torch.from_numpy(y_test)#.type(torch.LongTensor)

    train = torch.utils.data.TensorDataset(x_train, y_train)
    test = torch.utils.data.TensorDataset(x_test, y_test)

    train_loader = torch.utils.data.DataLoader(train, batch_size=batch_size, shuffle=True)
    test_loader = torch.utils.data.DataLoader(test, batch_size=batch_size, shuffle=True)

    ##模型加载
    model = ANN_prediction(train,test,num_inputs=114)  #114
    ErrorLoss = nn.MSELoss()

    #模型训练
    pltx = []
    plty = [[],[]]
    loss_list = []
    if atomName=='CA' or atomName=='CB' or atomName=='C' or atomName=='N':
        lr = 0.001   ##0.001 200
        L = 1000
    else:
        lr = 0.0001
        L = 1500
    optimizer = torch.optim.Adam(model.parameters(),lr=lr,weight_decay=0.003)  ##0 0 0.003
    for iteration in range(L):
        for j,(features,labels) in enumerate(train_loader):
            train = Variable(features)
            labels = Variable(labels)
            optimizer.zero_grad()
            outputs = model(train)
            loss = ErrorLoss(outputs,labels)
            loss.backward()
            optimizer.step()
            if j % 50 == 0:
                correct = 0
                total = 0
                for features, labels in test_loader:
                    test = Variable(features)
                    outputs = model(test)
                    dist = ErrorLoss(outputs,labels)
                    
                    total += 1
                    correct += dist
                t_loss = correct/total
                loss_list.append(loss.data)
            
        if iteration % 50 == 0:
            print('Epoch:{} Loss:{} test_loss:{}'.format(iteration, loss.data, t_loss))
            pltx.append(iteration)
            plty[0].append(loss.data.tolist()),plty[1].append(t_loss.tolist())

    #绘制loss值变化的图像
    plt.figure(figsize=(8,8))
    plt.plot(pltx,plty[0],'-r',label='training-set')
    plt.plot(pltx,[np.median(plty[0])]*len(pltx),'-*r')
    plt.plot(pltx,plty[1],'-b',label='testing-set')
    plt.plot(pltx,[np.min(plty[1])]*len(pltx),'-*b')
    plt.title('loss value in training-set & testing-set('+atomName+')')
    plt.xlabel('interation(epoch)')
    plt.ylabel('loss value')
    plt.legend()
    plt.savefig('Working/Sparta+/results/loss/figure_'+atomName+'_loss_2.png')
    plt.close()
    torch.save(model,\
        'Working/Sparta+/results/ANN_'+atomName+'_models.pth.tar')
    print('model for ' +atomName+ ' is constructed!')
    print()
    print('#'*50)
    print()

def data_Normalization(data:pd.DataFrame,name:str):
    '''
    处理离群值，可能是正常的离群值，也可能是shift检测出现问题导致的离群值
    '''
    min_threshold = data[name].quantile(0.005)
    max_threshold = data[name].quantile(0.995)
    data = data.loc[(data[name] >= min_threshold) & (data[name] <= max_threshold)]
    return data

class ANN_prediction(nn.Module):
    '''
    ANN预测模型
    '''
    def __init__(self,train_set,test_set,batch_size:int=256,\
        num_inputs:int=114,num_outputs:int=1,num_hiddens:int=30,num_epochs:int=1500,lr:float=0.001):
        '''
        '''
        super(ANN_prediction,self).__init__()
        self.fc1 = nn.Linear(num_inputs,num_hiddens)
        self.relu1 = nn.Tanh()
        self.fc2 = nn.Linear(num_hiddens,num_hiddens)
        self.relu2 = nn.ReLU()
        self.fc3 = nn.Linear(num_hiddens,num_hiddens)
        self.relu3 = nn.ReLU()
        self.fc4 = nn.Linear(num_hiddens,num_outputs)
        """
        self.train_set = train_set
        self.test_set = test_set
        self.batch_size = batch_size
        self.num_inputs, self.num_outputs, self.num_hiddens = num_inputs,num_outputs,num_hiddens
        self.W1 = torch.tensor(np.random.normal(0, 0.01, (num_inputs, num_hiddens)), dtype=torch.float)
        self.b1 = torch.zeros(num_hiddens, dtype=torch.float)
        self.W2 = torch.tensor(np.random.normal(0, 0.01, (num_hiddens, num_outputs)), dtype=torch.float)
        self.b2 = torch.zeros(num_outputs, dtype=torch.float)
        self.num_epochs = num_epochs #5
        self.lr = lr#100.0
        """

    def forward(self,x):
        out = self.fc1(x)
        out = self.relu1(out)
        out = self.fc2(out)
        out = self.relu2(out)
        out = self.fc3(out)
        out = self.relu3(out)
        out = self.fc4(out)
        return out

    def Run_ANN(self):
        
        """
        #train_iter,test_iter = d2l.load_data_fashion_mnist(batch_size)
        params = [self.W1, self.b1, self.W2, self.b2]
        for param in params:
            param.requires_grad_(requires_grad=True)
        loss = torch.nn.CrossEntropyLoss()
        d2l.train_ch3(self.net,self.train_set,self.test_set,loss,self.num_epochs,self.batch_size,params,self.lr)
        #d2l.train_ch3(self.net,self.train_set,self.test_set,loss,self.num_epochs,self.batch_size,params,self.lr)
        """
    """
    def tanh_i(self,X):
        return torch.tanh(input=X,other=torch.tensor(0,0))

    def net(self,X):
        H = self.tanh_i(torch.matmul(X,self.W1)+self.b1)
        return torch.matmul(H,self.W2)+self.b2
    """
    
if __name__ == '__main__':
    '''
    '''
    atomNameList = ['HA','HN','CB','CA','C','N']
    for atomName in atomNameList:
        main(atomName)
