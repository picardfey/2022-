p = [2 5; 3 6; 12 2; 1 6; 9 2; 8 12; 4 7; 7 9]'; % 特征数据X1,X2 
t = [10 18 24 6 18 96 28 63]; % 样本值 
net = newff(p, t, 20); % 创建一个BP神经网络 ff=FeedForward 
net = train(net, p, t); % 用p,t数据来训练这个网络