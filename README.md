# CA_SEIR
元胞自动机模拟病毒传染（SEIR模型）可视化
#### 功能描述
设置一定大小的人数、四种人群（四种状态）、传染概率、潜伏时间、治愈时间、免疫时间，即可对其进行模拟也可以更改规则和人数、人群，实现自己想要的模拟

#### 环境配置&需要的包
* Python 3.x+Anaconda+Pycharm
* ramdom、numpy、sys
* pygame
* matplotlib

#### 如何更改参数或制定自己想要的规则
- 1.更改传染病模型人群种类（状态）代码2-6行注释
- 2.更改规则，对每一个类别的人群，确定相应规则 代码9-12行，98-128行
- 3.更改人群总数 代码226行（后两个参数的乘积）+246行（退出循环条件）
- 4.更改概率，代码20-26行

#### 相关示例



![fig.1](https://github.com/Windxy/CA_SEIR/blob/image/1.png)
![fig.2](https://github.com/Windxy/CA_SEIR/blob/image/2.png)
![fig.3](https://github.com/Windxy/CA_SEIR/blob/image/3.png)
![fig.4](https://github.com/Windxy/CA_SEIR/blob/image/4.png)
![fig.5](https://github.com/Windxy/CA_SEIR/blob/image/5.png)
