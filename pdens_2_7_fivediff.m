clc,clear
%u(0,y)=siny+cosy,u(1,y)=e(siny+cosy) 0<=y<=1
%u(x,0)=e^x,u(x,1)=e^x(sin1+cos1) 0<x<1
%精确解为e^x(siny+cos(y)
%取(h1,h2)=(1/4,1/4),(1/8,1/8),(1/16,1/160,(1/32,1/32),(1/64,1/64)
%% Gauss-Seidel迭代法
tic
clear,clc
N=20000;a=0.4;
for l=0:4
    figure(l+1)
    m1=4*2^l;m2=4*2^l;h1=1/m1;h2=1/m2;
    u_temp=zeros(m1+1,m2+1);
    %犯错，我只给了角点的初值条件
    u_temp(1:m1+1,1)=exp(0:h1:1);u_temp(1:m1+1,m2+1)=exp(0:h1:1)*(sin(1)+cos(1));
    u_temp(1,1:m2+1)=sin(0:h2:1)+cos(0:h2:1);u_temp(m1+1,1:m2+1)=exp(1)*(sin(0:h2:1)+cos(0:h2:1));
    u=u_temp;
    for k=1:N
        for i=2:m1
            for j=2:m2
                u_temp(i,j)=0.5*(u(i,j-1)/h2^2+...
                u(i-1,j)/h1^2+u(i+1,j)/h1^2+u(i,j+1)/h2^2)/(1/h1^2+1/h2^2);
                
            end
        end
        u_temp=a*u_temp+(1-a)*u;
        oeps=max(max(abs(u_temp-u)));
        if oeps<10^-10/2
            fprintf('迭代到第%d次时，误差为%d',k,oeps);
            break;
        else
            u=u_temp;
        end
    end
     [x,y]=meshgrid(0:h1:1);
     meshz(x,y,u_temp);%数值解曲面图
     for i=0:m1
         for j=0:m2
               uq(i+1,j+1)=exp(i*h1)*(sin(j*h2)+cos(j*h2));
         end
     end
end
z=abs(u_temp-uq);
figure(6)
meshz(x,y,z) %误差曲面图
mesh(x,y,uq) %精确解曲面图
toc
%% 直接构造矩阵求解
clear,clc
% pic_num=1;
for l=2:4
    m1=4*2^l;m2=4*2^l;h1=1/m1;h2=1/m2;
    C=diag(repmat([2/h1^2+2/h2^2],1,m1-1))+diag(repmat([-1/h1^2],1,m1-2),1)...
        +diag(repmat([-1/h1^2],1,m1-2),-1);a=-ones(m1-1,1)/h2^2;D=diag(a);
    A(1:m1-1,1:2*m1-2)=[C,D];
    for j=2:m2-2
        A((m1-1)*(j-1)+1:(m1-1)*j,1+(j-2)*(m1-1):(j+1)*(m1-1))=[D,C,D];
    end
    A((m2-2)*(m1-1)+1:(m2-1)*(m1-1),(m2-3)*(m1-1)+1:(m2-1)*(m1-1))=[D,C];
    for i=1:m1-1
        xi=i*h1;
        u0(i)=exp(xi);um2(i)=exp(xi)*(sin(1)+cos(1));
    end
    for j=1:m2-1
        yj=j*h2;
        ff(1:m1-1,j)=zeros(m1-1,1);
        ff(1,j)=(sin(yj)+cos(yj))/h1^2;
        ff(m1-1,j)=(exp(1)*(sin(yj)+cos(yj)))/h1^2;
    end
    f=ff(:);
    f(1:m1-1)=f(1:m1-1)-D*u0';
    f(end-m1+2:end)=f(end-m1+2:end)-D*um2';
    u=A^-1*f;
    [x,y]=meshgrid(h1:h1:1-h1);
    uq=exp(x).*(sin(y)+cos(y));
    u=reshape(u,m1-1,m2-1);
    u=u';
    
%     meshz(x,y,u);
    z=abs(u-uq);
    surf(x,y,z)
    hold on
%     a1='误差曲面图,当h=';a2='时';
%     title(sprintf('%s%d%s',a1,h1,a2));
%     drawnow;
%     F=getframe(gcf);
%     I=frame2im(F);
%     [I,map]=rgb2ind(I,256);
%     if pic_num == 1
%         imwrite(I,map,'test.gif','gif', 'Loopcount',inf,'DelayTime',0.2);
%     else
%         imwrite(I,map,'test.gif','gif','WriteMode','append','DelayTime',0.2);
%     end
%     pic_num = pic_num + 1;
uu(l+1,1)=u(2^l,2^l);uu(l+1,2)=u(2^l,2^l*2);uu(l+1,3)=u(2^l,2^l*3); %部分节点
uu(l+1,4)=u(2^l*3,2^l);uu(l+1,5)=u(2^l*3,2^l*2);uu(l+1,6)=u(2^l*3,2^l*3);
end
figure(6)
surf(x,y,u)
figure(7)
surf(x,y,uq)