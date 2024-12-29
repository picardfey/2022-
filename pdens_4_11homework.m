clear,clc
%紧差分格式
%c=1,f(x,t)=(t^2-x^2)sin(xt);0,x,0,sint L=1,T=1;uq=sin(xt)
format long e
for l=1:4
    m=5*2^l;n=m^2;h=1/m;t=1/n;s=1*t/h;u=[];
    u(1,:)=zeros(m-1,1);
    A0=diag(repmat([10-6*s^2],1,m-1))+diag(repmat([1+3*s^2],1,m-2),1)...
        +diag(repmat([1+3*s^2],1,m-2),-1);
    for i=1:m-1
        f0(i)=1/2*t^2*(0+0+0)+t*(10*i*h+(i-1)*h+(i+1)*h);
    end
    f0(1)=f0(1)+(3*s^2-1)*0+(1+3*s^2)*0;
    f0(m-1)=f0(m-1)+(3*s^2-1)*sin(t);
    u(2,:)=Thomas((1-3*s^2)*ones(m-2,1),(10+6*s^2)*ones(m-1,1),(1-3*s^2)*ones(m-2,1),f0);
    A1=diag(repmat([5/3],1,m-1))+diag(repmat([1/6],1,m-2),1)...
        +diag(repmat([1/6],1,m-2),-1);
    A2=diag(repmat([-5/6-s^2],1,m-1))+diag(repmat([-1/12+s^2/2],1,m-2),1)...
        +diag(repmat([-1/12+s^2/2],1,m-2),-1);
    for k=1:n-1
        tk=k*t;f=zeros(m-1,1);
        for i=1:m-1
            xi=i*h;xj=(i-1)*h;xk=(i+1)*h;
            f(i)=t^2/12*((tk^2-xj^2)*sin(tk*xj)+10*(tk^2-xi^2)*sin(tk*xi)+(tk^2-xk^2)*sin(tk*xk));
        end
        f(m-1)=f(m-1)+(s^2/2-1/12)*sin(tk+t)+1/6*sin(tk)+(-1/12+s^2/2)*sin(tk-t);
        d=f+A1*u(k+1,:)'+A2*u(k,:)';
        u(k+2,:)=Thomas((1/12-s^2/2)*ones(m-2,1),(5/6+s^2)*ones(m-1,1),(1/12-s^2/2)*ones(m-2,1),d);
    end
    for i=1:n+1
        for j=1:m-1
            uq(i,j)=sin(h*j*t*(i-1));
        end
    end
    z=abs(u-uq);
    umax(l)=max(max(z));
    [x,y]=meshgrid(0:t:1,h:h:1-h);
    
    if l<=2
        figure(2)
        surf(x,y,z')
        
        title('误差曲面图')
        hold on
        if l==1
            figure(1)
            surf(x,y,u')
            title('数值解的曲面图')
            for i=1:10
                table1(i,1)=u(1+10*i,m/2);  %第一个表
                table1(i,2)=uq(1+10*i,m/2);
                table1(i,3)=z(1+10*i,m/2);
            end
        else
            for i=1:10
                table2(i,1)=u(1+10*i,m/2);%第二个表
                table2(i,2)=uq(1+10*i,m/2);
                table2(i,3)=z(1+10*i,m/2);
            end
        end
    end
end
figure(4)
surf(x,y,uq')
title('精确解的曲面图')
for i=1:3
    uyz(i)=umax(i)/umax(i+1); %第三个表
end