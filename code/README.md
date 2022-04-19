## 关于我的一部分工作的总结





## 纲目

这是关于文件路径以及各个文件作用的解释性文件

###### 文件夹

1. 我们的工作目录在white: /home/zhangjs/Working/Sparta+路径下，最终结果以及总结在/home/zhangjs/Documents/spartaPlusResults路径下

2. 对于工作目录

   data文件夹是训练数据集以及测试数据集的来源，同时，我们的程序将这部分中的原本来自于shiftx的pdb文件以及shift文件转化为ANN可以输入的矩阵形式，保存为csv格式

   parameters文件夹下存储了sparta+用于计算一些特征所需要的数据

   results文件夹下包括了我们输出的模型，以及绘制的散点图，直方分布图，小提琴图以及loss变化的折线图

   sparta+文件夹下包括了sparta+原本的ANN程序，采用C++

3. 对于最终结果目录

   01是我们的初始参数设定，没有去除离群值，没有采取L2范式化，包含了所有的114种特征值

   02是我们在01基础上，去除上下各0.05的离群值产生的结果

   03是我们在02基础上，采取L2范式化的结果

   04是我们在03基础上，将离群值去除减少到上下各0.005产生的结果

   05是我们在04的基础上，放弃氢键能量产生的结果

   06是我们在04的基础上，放弃计算ring current shift以及EF shift的结果

   所有结果均包括生成的模型，绘图，以及10次训练验证所产生的RMSD均值。PPT内时work report汇报文件可作为参考，excel文件内包括了各个预测器之间的比较。

   loss文件夹内包含了我之前采取的lr和epoch所产生的loss在训练集和测试集中的变化

   violin文件夹内没有什么很有意义的东西，自己画着玩的

###### python文件

​		/home/zhangjs/Working/Sparta+路径下有我们最终的的所有code文件，采用python书写

请保证您的工作环境包含以下这些安装包：

* pytorch

* pandas

* numpy

* matplotlib

* seaborn


具体的安装手段请自行查找，其中pytorch会比较麻烦，请仔细阅读官网的文件并且保证import不会出错（如果出错，按照错误安装对应的包即可）

setting.py是结构解释文件，你可以认为这是两个对象的定义文件，包括PDB结构和环结构

PDB.py是PDB的读取和计算文件，我们的特征计算算法基本上传承自sparta+，如果你希望采用更先进或者正确的算法，还请在这其中进行。你需要阅读的仅是\_init\_()函数以及每个其他函数的注释代码，请务必弄清楚所有的self.name存储的结构形式。请根据注释寻找需要替换的算法函数。主要是氢键计算部分，ring current部分，HNS2部分以及EF部分。

![image-20211128140948165](C:\Users\sunnyflyingsky\AppData\Roaming\Typora\typora-user-images\image-20211128140948165.png)

![image-20211128140958259](C:\Users\sunnyflyingsky\AppData\Roaming\Typora\typora-user-images\image-20211128140958259.png)

![image-20211128141014898](C:\Users\sunnyflyingsky\AppData\Roaming\Typora\typora-user-images\image-20211128141014898.png)

dataTransform.py这是特征值提取的脚本文件，他会将已经加氢的PDB文件转化问ANN可以解释的特征矩阵形式，请直接到脚本末尾找到 if  \_\_name\_\_ 并且从此处开始阅读。其中函数pred\_set\_()和train\_set_()均包含了输入原子，输入文件夹路径，输出文件夹路径这三个参数。

get_shift.py这是shift值提取文件，我们在dataTransform.py中调用它。

ANNPyTorch.py是我们模型训练部分。请直接到脚本末尾找到 if  \_\_name\_\_ 并且从此处开始阅读。所有的超参数我均已放入main()函数的参数表中，请在调用main()函数时修改这些参数。

predict.py是我们的模型验证和预测部分。请直接到脚本末尾找到 if  \_\_name\_\_ 并且从此处开始阅读。其中的数据处理流程同ANNPyTorch.py的部分，注释同样参考那边。



###### 数据文件以及图片文件

我们的特征矩阵文件保存在data路径下，它的结构如下

![image-20211128142213909](C:\Users\sunnyflyingsky\AppData\Roaming\Typora\typora-user-images\image-20211128142213909.png)

其中shift1是完整的shift值也就是 shift - random coil - ring current - EF

shift2是不完整的shift值，是shift - random coil



我们的loss折线图，红线代表training的loss变化，蓝线代表testing的loss变化。其中我们截取了training loss值的中位数作为其波动的参考线，testing loss值的最小值作为其趋近的参考线。



我们的散点图已做了一些命名优化，原本的四个图1,2,3,4分别表示验证集的散点图，直方图，训练集的散点图，直方图，现已优化为train\_1，train\_2，pred\_1，pred\_2。散点图表示预测值和实验值的比较，x轴是实验值，y轴是预测值，红线代表x=y



model.pth.tar文件是由pytorch保存的模型文件，可以用pytorch调用，可以阅读相关文档了解如何保存一个训练好的模型以及调用一个已有的模型。



## 细节性tips

* 脚本文件中包含了很多注释，希望对您有用
* 目前没有优化到可以在终端进行的输入输出，但是已经尽可能地剥离出相关的参数供代码读者调整，一般都在最表层的函数中可以看到。
* 每次进行训练以及验证后，请备份模型和图片，否则我们的脚本将会覆盖这些结果
* 如果有新的数据集，可以考虑直接加入到对应的data路径下，注意加氢问题



## 一些后续的想法（在我咨询了了解更多的前辈之后）

* 关于模型的改进方面，可以考虑加入卷积层，具体如何进行卷积，可以参考之前老师给的pytorch学习网页。
* 采用一些机器学习方法来考证这114个参数的重要程度，类似于shiftx中做的，或者用反向传播算法计算他们的权重
* 前辈有提出一个观点，可能在这训练集中的两万个三氨基酸对的特征中，存在很多相似性很高的例子，这些例子可能对于模型训练出现过拟合现象比较明显，尤其是对于那部分去除ring current和EF的模型而言。可以考虑进行降维或者聚类，保证每个例子具备参考性价值而不至于过多重复。
* 基于之前的一些数据结果去思考原因是一个比较好的出发点，或许可以考究里面存在这些变化的原因。
* 可以加入一些epoch的灵活结束手段，比如预测集的loss在10次迭代中不再减少之类的，这样可以保证过拟合现象的降低。
* 我之前并没有考察过训练集和验证集的交叉性，一般来说，为了学习器的泛化性能更强，验证集内尽可能不包含与训练集特征相同的例子，也就是说，后续可以考虑将训练集和验证集先整合起来，做完聚类后，拿出其中一些类别做验证，训练时要完全避免这部分特征的例子。



