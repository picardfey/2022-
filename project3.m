clc,clear
f=@(x,t)(x.^2-x-(4*x+1).*t);
uqz=@(x,t)(x.*(x-1).*t);ax=@(x)(x+1);bx=@(x)(1./(x+1).^2);
cx=@(x)(1./(x+1));af=@(x,t)((x+1).*(x.^2-x-(4*x+1).*t));
for l=1:3
    m=5*2^l;n=2*4^l;h=1/m;t=1/n;xx=h:h:1-h;
    u=[];uq=[];
    u(1,:)=zeros(1,m-1);uq(1,:)=zeros(1,m-1);
    for k=1:n
        A=(1+h/2)/12/t+h/24/t*cx(xx-h)+bx(xx-h/2)/24-ax(xx-h/2)/2/h^2;
        B=5/6/t+(ax(xx+h/2)+ax(xx-h/2))/2/h^2-(bx(xx+h/2)+bx(xx-h/2))/24;
        C=(1-h/2)/12/t-h/24/t*cx(xx+h)+bx(xx+h/2)/24-ax(xx+h/2)/2/h^2;
        D=(1+h/2)/12/t+h/24/t*cx(xx-h)-bx(xx-h/2)/24+ax(xx-h/2)/2/h^2;
        E=5/6/t-(ax(xx+h/2)+ax(xx-h/2))/2/h^2+(bx(xx+h/2)+bx(xx-h/2))/24;
        F=(1-h/2)/12/t-h/24/t*cx(xx+h)-bx(xx+h/2)/24+ax(xx+h/2)/2/h^2;
        fx=5/6*f(xx,(k-1/2)*t)+(f(xx-h,(k-1/2)*t)+f(xx+h,(k-1/2)*t))/12+...
            h^2/24*(af(xx+h,(k-1/2)*t)-af(xx-h,(k-1/2)*t));
        uq(k+1,:)=uqz(xx,k*t);
        AL=diag(A(2:end),-1)+diag(B)+diag(C(1:end-1),1);
        AR=diag(D(2:end),-1)+diag(E)+diag(F(1:end-1),1);
        d=AR*u(k,:)'+fx';
        u(k+1,:)=Thomas(A(2:end),B,C(1:end-1),d);
    end
    z=abs(u-uq);
    %figure(l) %画误差曲面图
    [x,y]=meshgrid(h:h:1-h,0:t:1);
    surf(x,y,z)
    xlabel('x');ylabel('t');
    title('误差曲面图')
    hold on
    umax(l)=max(max(z));
end
umax(1:end-1)./umax(2:end)