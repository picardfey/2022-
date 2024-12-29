clear,clc;
epsilon=1e-8;
m1=16;m2=16;h1=2/m1;h2=2/m2;
C=diag(repmat([2/h1^2+2/h2^2],1,m1-1))+diag(repmat([-1/h1^2],1,m1-2),1)...
    +diag(repmat([-1/h1^2],1,m1-2),-1);a=-ones(m1-1,1)/h2^2;D=diag(a);
A(1:m1-1,1:2*m1-2)=[C,D];
for j=2:m2-2
    A((m1-1)*(j-1)+1:(m1-1)*j,1+(j-2)*(m1-1):(j+1)*(m1-1))=[D,C,D];
end
A((m2-2)*(m1-1)+1:(m2-1)*(m1-1),(m2-3)*(m1-1)+1:(m2-1)*(m1-1))=[D,C];
for j=1:m2-1
    yj=j*h2;
    for i=1:m1-1
        xi=i*h1;
        u0(i)=xi^2;um2(i)=-xi^2;
        ff(i,j)=(xi^2*pi^2/4-2)*cos(pi*yj/2);
        uq(i,j)=xi^2*cos(pi*yj/2);
    end
    ff(1,j)=ff(1,j)+0/h1^2;
    ff(m1-1,j)=ff(m1-1,j)+4*cos(pi*yj/2)/h1^2;   
end
f=ff(:);
f(1:m1-1)=f(1:m1-1)-D*u0';
f(end-m1+2:end)=f(end-m1+2:end)-D*um2';
[x,y]=meshgrid(h1:h1:2-h1);
surf(x,y,uq)
figure(1)
tic
[u_pcg,k1,obj1]=pcg_solver(A,f,epsilon);
toc
u_pcg=reshape(u_pcg,m1-1,m2-1);
z1=abs(u_pcg-uq);
meshz(x,y,z1)
tic
[u_gauss,k2,obj2]=Gauss_Seidel(A,f,epsilon);
toc
u_gauss=reshape(u_gauss,m1-1,m2-1);
z2=abs(u_gauss-uq);
figure(2)
% meshz(x,y,z2);
u=A^-1*f;

u=reshape(u,m1-1,m2-1);
z3=abs(u-uq);
meshz(x,y,z3);
% plot(1:30,obj2(1:30))
% hold on
% plot(1:30,obj1)
% legend('gauss迭代的目标函数','pcg迭代的目标函数')