function [u,i,obj]=pcg_solver(A,b,epsilon)
u0=zeros(length(b),1);
u=u0;i=0;
[L,~]=LU(A,length(b));
D=diag(diag(A));
B=L*D*L';
b=B^-1/2*b;A=B^-1/2*A*B^-1/2;
coa=cond(A);
r=b-A*u;
p=r;%前进方向
resid_new=r'*r; %残差由内积确定
resid_new0=resid_new;
obj=1/2*u'*A*u-b'*u;
while resid_new/resid_new0>epsilon
    alpha=resid_new/(p'*A*p);
    u=u+alpha*p;
    obj=[obj,1/2*u'*A*u-b'*u];
    r=r-alpha*A*p;%新的残差
    resid_old=resid_new;
    resid_new=r'*r;
    beta=resid_new/resid_old;
    p=r+beta*p;  %下一步的方向
    i=i+1;
end
u=B^-1/2*u;
end