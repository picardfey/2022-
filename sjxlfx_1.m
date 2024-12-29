clc,clear
%% 画随机游走图
y(1)=0;n=100;t=10; %按比例缩小随机游走过程
e=normrnd(0,1/sqrt(n),1,n*t);
plot(e)
for i=1:n*t
    y(i+1)=y(i)+e(i);
end
% plot(1:n*t+1,y)
xlabel('i');ylabel('y(i)');
%% y(t)=y(t-1)+u(t),u(t)~ARMA(1,2)
clc,clear
y(1)=0;u(1)=0;n=10000;
noise=normrnd(0,1,1,n);
for t=1:n-1
    u(t+1)=u(t)+noise(t+1)+noise(t);
    y(t+1)=y(t)+u(t+1);
end
plot(1:n,y)
adftest(y)
%% y(t)=0.6*y(t-1)+u(t),u(t)~ARMA(1,2)
clc,clear
y(1)=0;u(1)=0;n=10000;
noise=normrnd(0,1,1,n);
for t=1:n-1
    u(t+1)=u(t)+noise(t+1)+noise(t);
    y(t+1)=0.6*y(t)+u(t+1);
end
plot(1:n,y)
adftest(y)