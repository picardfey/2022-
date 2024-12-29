clear,clc
for i=1:4
    m=80*2^i;
    a=0;b=pi;h=(b-a)/m;
    alpha=-1;beta=-exp(pi);
    lamda1=0;lamda2=0;
    A=-1*ones(m,1);Y=A;
    xm=a:h:b;
    B=(2+h^2)*ones(m+1,1);
    B(1)=1+lamda1*h;B(end)=1+lamda2*h;
    uq=exp(xm).*sin(xm);
    d=h^2*(exp(xm).*(sin(xm)-2*cos(xm)));
    d(1)=h*alpha;d(end)=h*beta;
    u=Thomas(A,B,Y,d)';
    z=abs(u-uq);
    plot(xm,z)
    hold on
end
 legend('h=pi/160','h=pi/320','h=pi/640','h=pi/1280')