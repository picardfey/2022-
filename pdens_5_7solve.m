clear,clc
%a=1,f(x,y,t)=-3/2e^(1/2(x+y)-t);初值条件；边界条件
f=@(x,y,t)((x^2+y^2)*t^2*sin(x*y*t)+x*y*cos(x*y*t));
uqz=@(x,y,t)(sin(x*y*t));uqright=@(y,t)(sin(y*t));uqup=@(x,t)(sin(x*t));
%m1=20;m2=20;n=20;
m1=40;m2=40;n=40;
h1=1/m1;h2=1/m2;t=1/n;rx=t/h1^2;ry=t/h2^2;
u=zeros(m1+1,m2+1,n+1);
for j=1:m2+1
    for k=1:n+1
        u(m1+1,j,k)=uqright((j-1)*h2,(k-1)*t);
    end
end
for k=2:n+1
    tk=(k-1)*t; ustar=[];
    for j=1:m2-1
        yj=j*h2;
        ustar(1,j)=0;
        ustar(m1+1,j)=(1+ry)*uqright(yj,tk)-ry/2*uqright(yj-h2,tk)-ry/2*uqright(yj+h2,tk);%左右边界
        %构建三对角矩阵向量用追赶法求解
        B1=(1+rx)*ones(m1-1,1);A1=-rx/2*ones(m1-2,1);Y1=A1;
        %右边系数确定
        Pmatrix=diag((1-rx-ry)*ones(m1-1,1))+diag(rx/2*ones(m1-2,1),-1)+diag(rx/2*ones(m1-2,1),1);
        for i=1:m1-1
            xi=i*h1;
            fvector(i)=t/2*(f(xi,yj,tk)+f(xi,yj,tk-t))+ry/2*(u(i+1,j,k-1)+u(i+1,j+2,k-1));%tmd
        end
        fvector(1)=fvector(1)+rx/2*ustar(1,j)+rx/2*u(1,j+1,k-1);
        fvector(m1-1)=fvector(m1-1)+rx/2*ustar(m1+1,j)+rx/2*u(m1+1,j+1,k-1);
        fcomplete=fvector'+Pmatrix*u(2:m1,j+1,k-1);
        %求出中间过渡层
        ustar(2:m1,j)=Thomas(A1,B1,Y1,fcomplete);
    end 
    %下一步求第k+1层
    for i=1:m1-1
        xi=i*h1;
        u(i+1,1,k)=0;
        u(i+1,m2+1,k)=uqup(xi,tk);
        %构建三对角矩阵向量用追赶法求解
        B2=(1+ry)*ones(m2-1,1);A2=-ry/2*ones(m2-2,1);Y2=A2;
        %三对角矩阵的头尾移向
        fcomplete2=ustar(i+1,:);
        fcomplete2(1)=fcomplete2(1)+ry/2*u(i+1,1,k);
        fcomplete2(m2-1)=fcomplete2(m2-1)+ry/2*u(i+1,m2+1,k);
        u(i+1,2:m2,k)=Thomas(A2,B2,Y2,fcomplete2);
    end 
end
for i=0:m1
    for j=0:m2
        for k=0:n
            uq(i+1,j+1,k+1)=uqz(i*h1,j*h2,k*t);
        end
    end
end
uerror=abs(u-uq);
for i=1:10  %填表，找部分数值结果
    table(i,1)=u(m1/2+1,m2/2+1,i*n/10+1);
    table(i,2)=uq(m1/2+1,m2/2+1,i*n/10+1);
    table(i,3)=uerror(m1/2+1,m2/2+1,i*n/10+1);
end 