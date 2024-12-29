%抛物方程的Euler显示格式
clear,clc
a=1;
for i=0:3
    for j=0:3
%         figure(i+j+1)
        m=10*2^i;n=200*4^j;h=1/m;t=1/n;
        u=zeros(m-1,n+1);
        r=a*t/h^2;A=diag(repmat([1-2*r],1,m-1))+diag(repmat([r],1,m-2),1)...
        +diag(repmat([r],1,m-2),-1);
        u(:,1)=exp(h:h:1-h);
        for k=1:n
            b=zeros(m-1,1);b(1)=r*exp((k-1)*t);b(m-1)=r*exp(1+(k-1)*t);
            u(:,k+1)=A*u(:,k)+b;
        end
       [x,y]=meshgrid(0:t:1,h:h:1-h);
       uq=exp(x+y);
       z=abs(u-uq);
       mesh(x,y,z)
    end
end

mesh(x,y,u)
        