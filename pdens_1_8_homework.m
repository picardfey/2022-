clear,clc
format long e
a=0;b=pi/2;
%for j=0:3
m=8;h=(b-a)/m;
A(1)=a;A(m+1)=b
for i=1:m-1
    A(i+1)=a+i*h;
    diag0(i)=2+h^2*(a+i*h-1/2)^2;
    d(i)=h^2*((a+i*h)^2-(a+i*h)+5/4)*sin(a+i*h);
    uq(i)=sin(a+i*h);
end
d(m-1)=d(m-1)+1;
g(1)=d(1)/diag0(1);w(1)=-1/diag0(1);
for i=2:m-2
    g(i)=(d(i)+g(i-1))/(diag0(i)+w(i-1));
    w(i)=-1/(diag0(i)+w(i-1));
end
g(m-1)=(d(m-1)+g(m-2))/(diag0(m-1)+w(m-2));
u(m-1)=g(m-1);
for i=m-2:-1:1
    u(i)=g(i)-w(i)*u(i+1);
end
for i=1:7
    un(i)=u(i*m/8);
    unq(i)=uq(i*m/8);
end
plot(a:h:b,[0,un,1],'-o')
legend('')
hold on
x=0:0.01:pi/2;
plot(x,sin(x))
legend('uh(x)','u(x)')
xlabel('x');
ylabel('u');
% plot(a:2^j*h:b,[0,abs(un-unq),0],'o-')
% hold on
%end
% xlabel('x')
% ylabel('|u(x)-uh(x)|')
% legend('h=pi/16','h=pi/32','h=pi/64','h=pi/128');
