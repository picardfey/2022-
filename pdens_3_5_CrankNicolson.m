%%
%为了得到结果，请在需要的地方点击运行节
%书上例子
%抛物方程Crank-nicolson格式求解
%uq=e^x+t
clear,clc
format long e
h=1/400;t=1/400;r=t/h^2;m=1/h;n=1/t;
u(1,:)=exp(h:h:1-h);
 A=-r/2*ones(m-2,1);
 B=(1+r)*ones(m-1,1);
 Y=-r/2*ones(m-2,1);
 APY=diag(repmat([1-r],1,m-1))+diag(repmat([r/2],1,m-2),1)...
        +diag(repmat([r/2],1,m-2),-1);
for k=2:n+1
    d(k-1,:)=zeros(1,m-1);d(k-1,1)=d(k-1,1)+r/2*(exp((k-2)*t)+exp((k-1)*t));
    d(k-1,m-1)=d(k-1,m-1)+r/2*(exp(1+t*(k-2))+exp(1+t*(k-1)));
    d(k-1,:)=d(k-1,:)+(APY*u(k-1,:)')';
    u(k,:)=Thomas(A,B,Y,d(k-1,:));
end
for i=1:n+1
    for j=1:m-1
        uq(i,j)=exp(h*j+t*(i-1));
    end
end
uz=abs(u-uq);
[x,y]=meshgrid(h:h:1-h,0:t:1);
meshz(x,y,uz)
%%
%作业题3.11 uq=e^xsin(1/2-t) 
clear,clc
format long e
for l=1:5
    m=5*2^l;n=5*2^l;h=1/m;t=1/n;r=2*t/h^2;
    u=[];uq=[];d=[];
    u(1,:)=exp(h:h:1-h)*sin(1/2);
    A=-r/2*ones(m-2,1);
    B=(1+r)*ones(m-1,1);
    Y=-r/2*ones(m-2,1);
    APY=diag(repmat([1-r],1,m-1))+diag(repmat([r/2],1,m-2),1)...
        +diag(repmat([r/2],1,m-2),-1);
    for k=2:n+1
        d(k-1,:)=-t*exp(h:h:1-h)*(cos(1/2-t*(k-3/2))+2*sin(1/2-t*(k-3/2)));
        d(k-1,1)=d(k-1,1)+r/2*(sin(1/2-(k-2)*t)+sin(1/2-(k-1)*t));
        d(k-1,m-1)=d(k-1,m-1)+r/2*exp(1)*(sin(1/2-t*(k-2))+sin(1/2-t*(k-1)));
        d(k-1,:)=d(k-1,:)+(APY*u(k-1,:)')';
        u(k,:)=Thomas(A,B,Y,d(k-1,:));
    end
    for i=1:n+1
        for j=1:m-1
            uq(i,j)=exp(h*j)*sin(1/2-(i-1)*t);
        end
    end
    erp=abs(u(end,:)-uq(end,:));
    einf(l)=max(max(u-uq));
    plot(h:h:1-h,erp);hold on
end
legend('h,\tau=(1/10,1//10)','h,\tau=(1/20,1//20)','h,\tau=(1/40,1//40)');
title('误差解曲线，当t=1时');
for i=1:4
    bili(i)=einf(i)/einf(i+1);%结果显示接近于6 表3.36
end
%%
clear,clc
format long e
for l=1:2        %开始填表外推
    m=50*2^l;n=50*2^l;h=1/m;t=1/n;r=2*t/h^2;
    u=[];uq=[];d=[];
    u(1,:)=exp(h:h:1-h)*sin(1/2);
    A=-r/2*ones(m-2,1);
    B=(1+r)*ones(m-1,1);
    Y=-r/2*ones(m-2,1);
    APY=diag(repmat([1-r],1,m-1))+diag(repmat([r/2],1,m-2),1)...
        +diag(repmat([r/2],1,m-2),-1);
    for k=2:n+1
        d(k-1,:)=-t*exp(h:h:1-h)*(cos(1/2-t*(k-3/2))+2*sin(1/2-t*(k-3/2)));
        d(k-1,1)=d(k-1,1)+r/2*(sin(1/2-(k-2)*t)+sin(1/2-(k-1)*t));
        d(k-1,m-1)=d(k-1,m-1)+r/2*exp(1)*(sin(1/2-t*(k-2))+sin(1/2-t*(k-1)));
        d(k-1,:)=d(k-1,:)+(APY*u(k-1,:)')';
        u(k,:)=Thomas(A,B,Y,d(k-1,:));
    end
    for i=1:n+1
        for j=1:m-1
            uq(i,j)=exp(h*j)*sin(1/2-(i-1)*t);
        end
    end
    for i=1:10
        erphalf(i,l)=u(1+i*n/10,m/2);% 表3.33和3.34数值解h=1/100、/h=1/200
        erphalf(i,3)=uq(1+i*n/10,m/2);%精确解
        erphalf(i,l+3)=abs(erphalf(i,l)-erphalf(i,3));%精确解-数值解的绝对值
    end
end
erphalf(:,6)=4/3*erphalf(:,2)-1/3*erphalf(:,1);%外推结果数值解
erphalf(:,7)=abs(erphalf(:,6)-erphalf(:,3));%外推误差
plot(0.1:0.1:1,erphalf(:,4),0.1:0.1:1,erphalf(:,5),0.1:0.1:1,erphalf(:,7));
legend('h=1/100,/tau=1/100','h=1/200,/tau=2/100','外推误差')

%%
%表3.35
clear,clc
format long e

m=200;n=2000;h=1/m;t=1/n;r=2*t/h^2;
u=[];uq=[];d=[];
u(1,:)=exp(h:h:1-h)*sin(1/2);
A=-r/2*ones(m-2,1);
B=(1+r)*ones(m-1,1);
Y=-r/2*ones(m-2,1);
APY=diag(repmat([1-r],1,m-1))+diag(repmat([r/2],1,m-2),1)...
    +diag(repmat([r/2],1,m-2),-1);
for k=2:n+1
    d(k-1,:)=-t*exp(h:h:1-h)*(cos(1/2-t*(k-3/2))+2*sin(1/2-t*(k-3/2)));
    d(k-1,1)=d(k-1,1)+r/2*(sin(1/2-(k-2)*t)+sin(1/2-(k-1)*t));
    d(k-1,m-1)=d(k-1,m-1)+r/2*exp(1)*(sin(1/2-t*(k-2))+sin(1/2-t*(k-1)));
    d(k-1,:)=d(k-1,:)+(APY*u(k-1,:)')';
    u(k,:)=Thomas(A,B,Y,d(k-1,:));
end
for i=1:n+1
    for j=1:m-1
        uq(i,j)=exp(h*j)*sin(1/2-(i-1)*t);
    end
end
for i=1:10
    erp(i,1)=u(1+n*i/10,100);
    erp(i,2)=uq(1+n*i/10,100);
    erp(i,3)=abs(erp(i,1)-erp(i,2));
end
%结果不如RICHARDSON外推法的误差小