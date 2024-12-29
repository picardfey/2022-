clear,clc
%紧致交替方向隐格式
f=@(x,y,t)((x^2+y^2)*t^2*sin(x*y*t)+x*y*cos(x*y*t));
uqz=@(x,y,t)(sin(x*y*t));uqright=@(y,t)(sin(y*t));uqup=@(x,t)(sin(x*t));
m1=10;m2=10;n=100;
m1=20;m2=20;n=400;
h1=1/m1;h2=1/m2;t=1/n;rx=t/h1^2;ry=t/h2^2;
u=zeros(m1+1,m2+1,n+1);
for j=0:m2
    for k=0:n
        u(m1+1,j+1,k+1)=uqright(j*h2,k*t);
        for i=0:m1
            uq(i+1,j+1,k+1)=uqz(i*h1,j*h2,k*t);
        end
    end
end
for k=1:n
    tk=k*t;
    for j=1:m2-1
        yj=j*h2;
        ustar(1,j)=0; %取边界条件
        ustar(m1+1,j)=(1/12-ry/2)*(uqright(yj-h2,tk)+uqright(yj+h2,tk))+(10/12+ry)*uqright(yj,tk);
        B1=(10/12+rx)*ones(m1-1,1);A1=(1/12-rx/2)*ones(m1-2,1);Y1=A1; %稀疏矩阵
        %最难的部分
        Pmatrixrow=[1/12+rx/2 10/12-rx 1/12+rx/2];
        Pmatrixcol=[1/12+ry/2;10/12-ry;1/12+ry/2];
        for i=1:m1-1
            for i1=1:3    %为了简洁被迫构造矩阵
                for i2=1:3
                    fmatrix(i1,i2)=f((i+i1-2)*h1,(j+i2-2)*h2,tk-t)+f((i+i1-2)*h1,(j+i2-2)*h2,tk);
                end
            end
            midmatrix=u(i:i+2,j:j+2,k);
            fvector(i)=t/2*[1/12 10/12 1/12]*fmatrix*[1/12;10/12;1/12]+Pmatrixrow*midmatrix*Pmatrixcol;
        end
        %头尾的移项
        fvector(1)=fvector(1)+(rx/2-1/12)*ustar(1,j);
        fvector(m1-1)=fvector(m1-1)+(rx/2-1/12)*ustar(m1+1,j);
        %求得中间层
        ustar(2:m1,j)=Thomas(A1,B1,Y1,fvector);
    end  
    %求第k+1层
    for i=1:m1-1
        u(i+1,1,k+1)=0;
        u(i+1,m2+1,k+1)=uqup(i*h1,tk);
        %三对角线
        B2=(10/12+ry)*ones(m2-1,1);A2=(1/12-ry/2)*ones(m2-2,1);Y2=A2;
        fvector2=ustar(i+1,:);
        fvector2(1)=fvector2(1)+(ry/2-1/12)*u(i+1,1,k+1);
        fvector2(m2-1)=fvector2(m2-1)+(ry/2-1/12)*u(i+1,m2+1,k+1);
        u(i+1,2:m2,k+1)=Thomas(A2,B2,Y2,fvector2);
    end
end
uerror=abs(u-uq);
for i=1:10  %填表，找部分数值结果
    table(i,1)=u(m1/2+1,m2/2+1,i*n/10+1);
    table(i,2)=uq(m1/2+1,m2/2+1,i*n/10+1);
    table(i,3)=uerror(m1/2+1,m2/2+1,i*n/10+1);
end
  