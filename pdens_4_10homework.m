clear,clc
%隐式差分格式
%c=1,f(x,t)=(t^2-x^2)sin(xt);0,x,0,sint L=1,T=1;uq=sin(xt)
for l=1:6
    u=[];
    m=5*2^l;n=5*2^l;
    h=1/m;t=1/n;s=1*t/h;
    u(1,:)=zeros(m-1,1);
    A0=diag(repmat([2-s^2],1,m-1))+diag(repmat([1/2*s^2],1,m-2),1)...
        +diag(repmat([1/2*s^2],1,m-2),-1);
    x0=h:h:1-h;f0=2*t*x0;f0(1)=f0(1)+1/2*s^2*(0+0);
    f0(m-1)=f0(m-1)+1/2*s^2*(sin(t)+0);
    AL0=4*eye(m-1)-A0;
    u(2,:)=AL0^-1*(A0*u(1,:)'+f0');
    A=diag(repmat([1+s^2],1,m-1))+diag(repmat([-1/2*s^2],1,m-2),1)...
        +diag(repmat([-1/2*s^2],1,m-2),-1);
    for k=1:n-1
        tk=k*t;f=[];
        f=2*u(k+1,:)+t^2*(tk^2-x0.^2).*sin(x0*tk);
        f(1)=f(1)+1/2*s^2*(0+0);f(m-1)=f(m-1)+1/2*s^2*(sin(tk+t)+sin(tk-t));
        u(k+2,:)=A^-1*(-A*u(k,:)'+f');
    end
    for i=1:n+1
        for j=1:m-1
            uq(i,j)=sin(h*j*t*(i-1));
        end
    end
    z=abs(u-uq);
        %figure(l)
        [x,y]=meshgrid(0:t:1,h:h:1-h);
        surf(x,y,z');
        hold on
    umax(l)=max(z(:));  %解的最大误差
end
%解它们的比
for i=1:5
    u_yz(i)=umax(i)/umax(i+1);
end
%%
clear,clc
%m=100;n=200;
m=200;n=100;
h=1/m;t=1/n;s=1*t/h;
u(1,:)=zeros(m-1,1);
A0=diag(repmat([2-s^2],1,m-1))+diag(repmat([1/2*s^2],1,m-2),1)...
    +diag(repmat([1/2*s^2],1,m-2),-1);
x0=h:h:1-h;f0=2*t*x0;f0(1)=f0(1)+1/2*s^2*(0+0);
f0(m-1)=f0(m-1)+1/2*s^2*(sin(t)+0);
AL0=4*eye(m-1)-A0;
u(2,:)=AL0^-1*(A0*u(1,:)'+f0');
A=diag(repmat([1+s^2],1,m-1))+diag(repmat([-1/2*s^2],1,m-2),1)...
    +diag(repmat([-1/2*s^2],1,m-2),-1);
for k=1:n-1
    tk=k*t;f=[];
    f=2*u(k+1,:)+t^2*(tk^2-x0.^2).*sin(x0*tk);
    f(1)=f(1)+1/2*s^2*(0+0);f(m-1)=f(m-1)+1/2*s^2*(sin(tk+t)+sin(tk-t));
    u(k+2,:)=A^-1*(-A*u(k,:)'+f');
end
for i=1:n+1
    for j=1:m-1
        uq(i,j)=sin(h*j*t*(i-1));
    end
end
z=abs(u-uq);

for i=1:10
    ux(i,1)=u(1+n/10*i,m/2);
    ux(i,2)=uq(1+n/10*i,m/2);
    ux(i,3)=z(1+n/10*i,m/2);
end

