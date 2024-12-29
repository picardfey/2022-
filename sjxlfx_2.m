clear,clc
N=1000;
for i=1:N
    T=50;
    white_noise=normrnd(0,2,1,T);
    fz=0;pfh=0;y=[];y(1)=0;
    for j=1:T
        y(j+1)=y(j)+white_noise(j);
        fz=fz+y(j)*white_noise(j);
        pfh=pfh+y(j)^2;
    end
    [p,pint,~,~,~]=regress(y(2:T+1)',y(1:T)');
    sigm2=0;
    for j=1:T
        sigm2=sigm2+((y(j+1)-p*y(j))^2)/(T-1);
    end
    tjlt(i)=fz/(sqrt(sigm2*pfh));
end
cdfplot(tjlt)
hold on
yzt=normrnd(0,1,1,N);
cdfplot(yzt)
legend('模拟结果','标准正态分布')