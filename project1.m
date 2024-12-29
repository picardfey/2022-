clear,clc
%4阶导数边值的紧差分格式
%精确解u(x)=xsinx+cosx
lamda1=2;lamda2=3;%f(x)=2xsinx;
format long e
for l=1:4
    m=5*2^(l-1);h=pi/m;  %网格剖分
    c0=-2*lamda1*h-h^3/3*lamda1+(24+10*h^2)/(h^2-12);
    d1=1-(12-h^2)/(h^2-12);
    dm=-1+(12-h^2)/(h^2-12);
    cm=-(24+10*h^2)/(h^2-12)+2*lamda2*h+h^3/3*lamda2;
    f1=-2*2*h-h^3/3*2+2*h^2*((-h)*sin(-h)+h*sin(h))/(h^2-12);
    fm=-2*(3+pi)*h+h^3/3*(-(3+pi)+2*pi)-2*h^2*((pi-h)*sin(pi-h)+(pi+h)*sin(pi+h))/(h^2-12);
    xq=0:h:1-h;
    A=(1/12-1/h^2)*ones(m,1);A(m)=dm;
    Y=(1/12-1/h^2)*ones(m,1);Y(1)=d1;
    B=(5/6+2/h^2)*ones(m+1,1);B(1)=c0;B(m+1)=cm;
    d=zeros(m+1,1);
    for i=0:m
        xi=(i-1)*h;xj=i*h;xk=(i+1)*h;
        d(i+1)=1/12*2*(xi*sin(xi)+10*xj*sin(xj)+xk*sin(xk));
        uq(i+1)=xj*sin(xj)+cos(xj);
    end
    d(1)=f1;d(m+1)=fm;
    u=Thomas(A,B,Y,d);
    z=abs(u'-uq);
    plot(0:h:pi,z);hold on
    for i=1:6     %取出需要的节点值
        un(l,i)=u(1+(i-1)*2^(l-1));
        uqn(i)=uq(1+(i-1)*2^(l-1));
    end
end
legend('h=pi/5','h=pi/10','h=pi/20','h=pi/40')
for i=1:4
    erro(i,:)=abs(un(i,:)-uqn);
end
figure(2)
plot(0:h:pi,uq)