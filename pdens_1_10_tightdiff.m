clc,clear
%紧差分格式求解
format long e
%%function e^x*sin(x) f(x)=e^x(sinx-2cosx)
a=0;b=pi;
for j=0:4
    m=10*2^j;h=(b-a)/m;
    for i=0:m
        x(i+1)=a+i*h;
        f(i+1)=exp(x(i+1))*(sin(x(i+1))-2*cos(x(i+1)));
    end
    for i=1:m-1
        d(i)=(f(i)+10*f(i+1)+f(i+2))/12;
    end
    dia=2/h^2+10/12;ga=-1/h^2+1/12;
    g(1)=d(1)/ga;w(1)=ga/dia;
    for i=2:m-2
        g(i)=(d(i)-ga*g(i-1))/(dia-ga*w(i-1));
        w(i)=ga/(dia-ga*w(i-1));
    end
    g(m-1)=(d(m-1)-ga*g(m-2))/(dia-ga*w(m-2));
    u(m-1)=g(m-1);uq(m-1)=exp(x(m))*sin(x(m));
    for i=m-2:-1:1
        u(i)=g(i)-w(i)*u(i+1);
        uq(i)=exp(x(i+1))*sin(x(i+1));
    end
end
u-uq;