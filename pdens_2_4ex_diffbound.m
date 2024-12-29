clc,clear
%创建稀疏矩阵
%精确解u(x,y)=e^xsin(piy)
xx=2;yy=1;m1=16;m2=8;
h1=xx/m1;h2=yy/m2;
f=zeros((m1+1)*(m2+1),1);
f(1:m1+1)=(pi^2-1)*exp(0:h1:xx)*sin(pi*0)+2/h2*(-pi*exp(0:h1:xx));
f(1)=f(1)+2/h1*0;f(m1+1)=f(m1+1)+2/h1*exp(2)*(1+2*0)*sin(pi*0);
for j=1:m2-1
    for i=0:m1
        f((m1+1)*j+i+1)=(pi^2-1)*exp(h1*i)*sin(pi*j*h2);
    end
    f((m1+1)*j+1)=f((m1+1)*j+1)+2/h1*0;f((m1+1)*(j+1))=f((m1+1)*(j+1))+2/h1*exp(2)*(1+2*j*h2)*sin(pi*j*h2);
end
for i=0:m1
    f(end-m1+i)=(pi^2-1)*exp(i*h1)*sin(pi*1)+2/h2*(-pi*exp(i*h1));
end
f(end-m1)=f(end-m1)+2/h1*0;f(end)=f(end)+2/h1*exp(2)*(1+2*1)*sin(pi*1);
row_index=[];colum_index=[];num=[];
for i=0:m2
    row_index(1+i*(3*m1+1):2+i*(3*m1+1))=i*(m1+1)+1;
    row_index((i+1)*(3*m1+1)-1:(i+1)*(3*m1+1))=(i+1)*(m1+1);
    for j=2:m1
         row_index(3*j-3+i*(3*m1+1):3*j-1+i*(3*m1+1))=j+i*(m1+1);
    end
end
% A=sparse(row_index,colum_index,num);
