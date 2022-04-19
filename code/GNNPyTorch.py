##
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from torch.autograd import Variable
import sys

class TwoLayerNet:
    def __init__(self, input_size, hidden_size, output_size, weight_init_std=0.01):
        # 权重初期化，函数参数依次为：输入层的神经元数，隐藏层的神经元数，输出层的神经元数
        self.params = {}
        self.params['W1'] = weight_init_std * np.random.randn(input_size, hidden_size)
        self.params['b1'] = np.zeros(hidden_size)
        self.params['W2'] = weight_init_std * np.random.randn(hidden_size, output_size)
        self.params['b2'] = np.zeros(output_size)
    
    def predict(self, x):
        W1, W2 = self.params['W1'], self.params['W2']
        b1, b2 = self.params['b1'], self.params['b2']
        a1 = np.dot(x, W1) + b1
        z1 = tanh_i(a1)
        a2 = np.dot(z1, W2) + b2
        y = softmax(a2)
        return y

    def loss(self, x, t):
        y = self.predict(x)
        return cross_entropy_error(y, t)

    def accuracy(self, x, t):
        y = self.predict(x)
        y = np.argmax(y, axis=1)
        t = np.argmax(t, axis=1)
        accuracy = np.sum(y == t) / float(x.shape[0])
        return accuracy
    # x:输入数据， t:目标向量的真实值

    def numerical_gradient(self, x, t):
        loss_W = lambda W: self.loss(x, t)
        grads = {}
        grads['W1'] = n_gradient(loss_W, self.params['W1'])
        grads['b1'] = n_gradient(loss_W, self.params['b1'])
        grads['W2'] = n_gradient(loss_W, self.params['W2'])
        grads['b2'] = n_gradient(loss_W, self.params['b2'])
        return grads

def tanh_i(X):
    return np.tanh(X)

def sigmoid(x):
    return 1 / (1 + np.exp(-x))

def softmax(x):
    if x.ndim == 2:
        x = x.T
        x = x - np.max(x, axis=0)
        y = np.exp(x) / np.sum(np.exp(x), axis=0)
        return y.T
    x = x - np.max(x)
    return np.exp(x) / np.sum(np.exp(x))

# 批处理版本的交叉熵误差
def cross_entropy_error(y, t):
    if y.ndim == 1:
        t = t.reshape(1, t.size)
        y = y.reshape(1, y.size)
    #如果t是ont-hot向量, 将其转换为正确答案标签的索引
    if t.size == y.size:
        t = t.argmax(axis=1)
    batch_size = y.shape[0]
    return -np.sum(np.log(y[np.arange(batch_size), t] + 1e-7))/batch_size

def n_gradient(f, x):
    h = 1e-4                  # 0.0001
    grad = np.zeros_like(x)
    it = np.nditer(x, flags=['multi_index'], op_flags=['readwrite'])
    while not it.finished:
        idx = it.multi_index
        tmp_val = x[idx]
        x[idx] = float(tmp_val) + h
        fxh1 = f(x)                      # f(x+h)
        x[idx] = tmp_val - h
        fxh2 = f(x)                      # f(x-h)
        grad[idx] = (fxh1 - fxh2) / (2*h)
        x[idx] = tmp_val # 还原值
        it.iternext()
    return grad


if __name__ == '__main__':
    ann_model = TwoLayerNet(114,30,1)
    batch_size = 256
    learning_rate = 0.1
    iters_num = 10000 # 迭代次数
    #数据处理
    train = pd.read_csv('Working/Sparta+/data/train_set.csv',sep=',',header=0,dtype=np.float32)
    for key in train.columns:
        train.fillna({key:train[key].mean()},inplace=True)
    #print(train)
    y = train.loc[:,train.columns=='shift'].values
    x = train.loc[:,train.columns!='shift'].values
    x_train,x_test,y_train,y_test = train_test_split(x,y,test_size=0.2,random_state=2021)
    train_size = x_train.shape[0]

    train_loss_list = []
    train_acc_list = []
    test_acc_list = []
    iter_per_epoch = max(train_size / batch_size, 1)
    for i in range(iters_num):
        batch_mask = np.random.choice(train_size, batch_size)
        x_batch = x_train[batch_mask]
        y_batch = y_train[batch_mask]
        # 梯度计算
        grad = ann_model.numerical_gradient(x_batch, y_batch)
        #grad = network.gradient(x_batch, t_batch)
        # 参数更新
        for key in ('W1', 'b1', 'W2', 'b2'):
            ann_model.params[key] -= learning_rate * grad[key]
            # 记录学习过程
        loss = ann_model.loss(x_batch, y_batch)
        train_loss_list.append(loss)
        if i % iter_per_epoch == 0:
            train_acc = ann_model.accuracy(x_train, y_train)
            test_acc = ann_model.accuracy(x_test, y_test)
            train_acc_list.append(train_acc)
            test_acc_list.append(test_acc)
            print("train acc, test acc | " + str(train_acc) + ", " + str(test_acc))

