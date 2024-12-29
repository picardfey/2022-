clear,clc
%二维波动方程的初边值问题
f=@(x,y,t)(-3*exp(x+y)*sin(t));posai=@(x,y)(exp(x+y));
uqleft=@(y,t)(exp(y)*sin(t));uqright=@(y,t)(exp(1+y)*sin(t));
uqdown=@(x,t)(exp(x)*sin(t));uqup=@(x,t)(exp(1+x)*sin(t));
uz=@(x,y,t)(exp(x+y)*sin(t));
rouf=@(x,y)(-exp(x+y));%求偏导赋值
m1=20;m2=20;n=20;
%m1=40;m2=40;n=40;
h1=1/m1;h2=1/m2;t=1/n;rx=t^2/h1^2;ry=t^2/h2^2;
u(:,:,1)=zeros(m1+1,m2+1);
uq=zeros(m1+1,m2+1,n+1);

for i=0:m1
    for j=0:m2
        for k=0:n
            u(m1+1,j+1,k+1)=uqright(j*h2,k*t);
            u(1,j+1,k+1)=uqleft(j*h2,k*t);
            uq(i+1,j+1,k+1)=uz(i*h1,j*h2,k*t);
            u(i+1,1,k+1)=uqdown(i*h1,k*t);
            u(i+1,m2+1,k+1)=uqup(i*h1,k*t);
        end
    end
end

ustar=zeros(m1+1,m2-1);
%先求边界起步层书5.68
for j=1:m2-1
    ustar(1,j)=(1+ry)*uqleft(j*h2,t)-ry/2*(uqleft((j-1)*h2,t)+uqleft((j+1)*h2,t));
    ustar(m1+1,j)=(1+ry)*uqright(j*h2,t)-ry/2*(uqright((j-1)*h2,t)+uqright((j+1)*h2,t));
    Bi_1=(1+rx)*ones(m1-1,1);Ai_1=-rx/2*ones(m1-2,1);Yi_1=Ai_1;
    fvector1=zeros(1,m1-1);
    for i=1:m1-1      %右边系数项
        fvector1(i)=u(i+1,j+1,1)+t*posai(i*h1,j*h2)-t^3/3*rouf(i*h1,j*h2)+t^2/2*f(i*h1,j*h2,t);
    end
    fvector1(1)=fvector1(1)+rx/2*ustar(1,j);
    fvector1(m1-1)=fvector1(m1-1)+rx/2*ustar(m1+1,j);
    ustar(2:m1,j)=Thomas(Ai_1,Bi_1,Yi_1,fvector1);
end
    %第二部求第一层的值
    for i=1:m1-1
        u(i+1,1,2)=uqdown(i*h1,t);
        u(i+1,m2+1,2)=uqup(i*h1,t);
        Bi_2=(1+ry)*ones(m2-1,1);Ai_2=-ry/2*ones(m2-2,1);Yi_2=Ai_2;
        fvector2=ustar(i+1,:);
        fvector2(1)=fvector2(1)+ry/2*u(i+1,1,2);
        fvector2(m2-1)=fvector2(m2-1)+ry/2*u(i+1,m2+1,2);
        u(i+1,2:m2,2)=Thomas(Ai_2,Bi_2,Yi_2,fvector2);
    end
%求剩下的层值

for k=1:n-1
    tk=t*k;ustar2=zeros(m1+1,m2-1);
    for j=1:m2-1
        yj=j*h2;
        ustar2(1,j)=1/2*((1+ry)*uqleft(yj,tk-t)-ry/2*(uqleft(yj-h2,tk-t)+uqleft(yj+h2,tk-t))...
            +(1+ry)*uqleft(yj,tk+t)-ry/2*(uqleft(yj-h2,tk+t)+uqleft(yj+h2,tk+t)));
        ustar2(m1+1,j)=1/2*((1+ry)*uqright(yj,tk-t)-ry/2*(uqright(yj-h2,tk-t)+uqright(yj+h2,tk-t))...
            +(1+ry)*uqright(yj,tk+t)-ry/2*(uqright(yj-h2,tk+t)+uqright(yj+h2,tk+t)));
        B2_1=(1+rx)*ones(m1-1,1);A2_1=-rx/2*ones(m1-2,1);Y2_1=A2_1;
         fvector3=zeros(1,m1-1);
        for i=1:m1-1
            fvector3(i)=u(i+1,j+1,k+1)+t^2/2*f(i*h1,yj,tk);
        end
        fvector3(1)=fvector3(1)+rx/2*ustar2(1,j);
        fvector3(m1-1)=fvector3(1)+rx/2*ustar2(m1+1,j);
        ustar2(2:m1,j)=Thomas(A2_1,B2_1,Y2_1,fvector3);
    end
    umid=zeros(m1-1,m2+1);
    for i=1:m1-1
        xi=i*h1;
        umid(i,1)=(uqdown(xi,tk+t)+uqdown(xi,tk-t))/2;
        umid(i,m2+1)=(uqup(xi,tk+t)+uqup(xi,tk-t))/2;
        B2_2=(1+ry)*ones(m2-1,1);A2_2=-ry/2*ones(m2-2,1);Y2_2=A2_2;
        fvector4=ustar2(i+1,:);
        fvector4(1)=fvector4(1)+ry/2*umid(i,1);
        fvector4(m2-1)=fvector4(m2-1)+ry/2*umid(i,m2+1);
        umid(i,2:m2)=Thomas(A2_2,B2_2,Y2_2,fvector4);
    end
   u(2:m1,2:m2,k+2)=2*umid(:,2:m2)-u(2:m1,2:m2,k);
end
uerror=abs(u-uq);
for i=1:10  %填表，找部分数值结果
    table(i,1)=u(m1/2+1,m2/2+1,i*n/10+1);
    table(i,2)=uq(m1/2+1,m2/2+1,i*n/10+1);
    table(i,3)=uerror(m1/2+1,m2/2+1,i*n/10+1);
end
