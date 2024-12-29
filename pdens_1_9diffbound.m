clc,clear
%导数边值问题格式求解
format long e
%%function 精确解u(x)=e^x;q(x)=1+sin(x),f(x)=e^xsin(x)
a=0;b=1;alp=0;beta=3*exp(1);
for j=0:4
    m=4*2^j;h=(b-a)/m;
    d(1)=h*alp;dia(1)=1+h*1+h^2*(1+sin(a))/2;
    g(1)=d(1)/dia(1);w(1)=-1/dia(1);
    d(m+1)=h*beta+(h^2*sin(b)*exp(b))/2;dia(m+1)=1+h*2+h^2*(1+sin(b))/2;
    for i=2:m
        x(i-1)=a+(i-1)*h;
        d(i)=h^2*exp(x(i-1))*sin(x(i-1));
        dia(i)=2+(1+sin(x(i-1)))*h^2;
        g(i)=(d(i)+g(i-1))/(dia(i)+w(i-1));
        w(i)=-1/(dia(i)+w(i-1));
    end
    g(m+1)=(d(m+1)+g(m))/(dia(m+1)+w(m));
    u(m+1)=g(m+1);uq(m+1)=exp(b);
    for i=m:-1:1
        u(i)=g(i)-w(i)*u(i+1);
        uq(i)=exp(a+(i-1)*h);
    end
    for k=1:3
        un(j+1,k)=u(k*2^j+1);
        uqn(j+1,k)=uq(k*2^j+1);
    end
end
A=abs(un-uqn);
for i=1:5
    einf(i)=max(A(i,:));
    plot(1/4:1/4:3/4,A(i,:))
    hold on
end
plot(0:h:1,uq)
hold on
plot(0:h:1,u,'or')
legend('精确解','数值解')