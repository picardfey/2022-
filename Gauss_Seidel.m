function [u,k,obj]=Gauss_Seidel(A,b,eps)
[~,n]=size(A);
x0=zeros(n,1);u=zeros(n,1);k=0;
obj=1/2*u'*A*u-b'*u;
for i=1:n
        u(i)=b(i);
        for j=1:i-1
            u(i)=u(i)-A(i,j)*u(j);
        end
         for j=i+1:n
            u(i)=u(i)-A(i,j)*x0(j);
         end
        u(i)=u(i)/A(i,i);
end
while norm(u-x0)>eps
    k=k+1;
    x0=u;
    for i=1:n
        u(i)=b(i);
        for j=1:i-1
            u(i)=u(i)-A(i,j)*u(j);
        end
         for j=i+1:n
            u(i)=u(i)-A(i,j)*x0(j);
         end
        u(i)=u(i)/A(i,i);
    end
    obj=[obj,1/2*u'*A*u-b'*u];
end