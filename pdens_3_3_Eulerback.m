%%
%书上的算例
clear,clc;
%抛物方程向后欧拉格式求解
%uq=e^x+t
format long e
h=1/40;t=1/1600;r=t/h^2;m=1/h;n=1/t;
A=-r*ones(m-2,1);
u(1,:)=exp(h:h:1-h);
 A=-r*ones(m-2,1);
 B=(1+2*r)*ones(m-1,1);
 Y=-r*ones(m-2,1);
for k=2:n+1
    d(k-1,:)=u(k-1,:)+t*0;d(k-1,1)=d(k-1,1)+r*exp((k-2)*t);d(k-1,m-1)=d(k-1,m-1)+r*exp(1+t*(k-2));
    u(k,:)=Thomas(A,B,Y,d(k-1,:));
end
for i=1:n+1
    for j=1:m-1
        uq(i,j)=exp(h*j+t*(i-1));
    end
end
uz=abs(u-uq);
% [x,y]=meshgrid(h:h:1-h,0:t:1);
% meshz(x,y,uz)
%%
%作业题第8题的算例
%需要用到Thomas函数，请一并放在同一文件夹目录下运行
clear,clc
format long e
for l=1:3
    u=[];d=[];uq=[];erp=[];
    m=5*2^l;n=25*4^l;a=2;
    h=1/m;t=1/n; r=a*t/h^2;
    A=-r*ones(m-2,1);
    u(1,:)=exp(h:h:1-h)*sin(1/2);
    A=-r*ones(m-2,1);
    B=(1+2*r)*ones(m-1,1);
    Y=-r*ones(m-2,1);
    for k=2:n+1
        d(k-1,:)=u(k-1,:)-t*exp(h:h:1-h)*(cos(1/2-(k-1)*t)+2*sin(1/2-(k-1)*t));
        d(k-1,1)=d(k-1,1)+r*sin(1/2-(k-2)*t);d(k-1,m-1)=d(k-1,m-1)+r*exp(1)*sin(1/2-(k-2)*t);
        u(k,:)=Thomas(A,B,Y,d(k-1,:));
    end
    for i=1:n+1
        for j=1:m-1
            uq(i,j)=exp(j*h)*sin(1/2-t*(i-1));
        end
    end
    erp=abs(u(end,:)-uq(end,:));
    plot(h:h:1-h,erp);hold on
end
legend('h,\tau=(1/10,1//100)','h,\tau=(1/20,1//400)','h,\tau=(1/40,1//1600)');
title('误差解曲线，当t=1时');
%%
%作业第8提
clear,clc
%m=100;n=10000;
m=100;n=40000;
%
a=2;
h=1/m;t=1/n; r=a*t/h^2;
A=-r*ones(m-2,1);
u(1,:)=exp(h:h:1-h)*sin(1/2);
A=-r*ones(m-2,1);
B=(1+2*r)*ones(m-1,1);
Y=-r*ones(m-2,1);
for k=2:n+1
    d(k-1,:)=u(k-1,:)-t*exp(h:h:1-h)*(cos(1/2-(k-1)*t)+2*sin(1/2-(k-1)*t));
    d(k-1,1)=d(k-1,1)+r*sin(1/2-(k-2)*t);d(k-1,m-1)=d(k-1,m-1)+r*exp(1)*sin(1/2-(k-2)*t);
    u(k,:)=Thomas(A,B,Y,d(k-1,:));
end
for i=1:n+1
    for j=1:m-1
        uq(i,j)=exp(j*h)*sin(1/2-t*(i-1));
    end
end
spant=(n/10)*(1:10)+1;
tab(:,1)=u(spant,m/2);%表3.29的三列
tab(:,2)=uq(spant,m/2);
tab(:,3)=abs(tab(:,1)-tab(:,2));

