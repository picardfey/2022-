clc,clear
f=@(x)((1/20*x^4-6)*x);
h=1/4;
A=[32+10/12 -16+1/12 0
    -16+1/12 32+10/12 -16+1/12
    0 -16+1/12 32+10/12];
A=[ 2/h^2+10/12 1/12-1/h^2 0
    1/12-1/h^2 2/h^2+10/12 1/12-1/h^2
    0 1/12-1/h^2 2/h^2+10/12 ];
b=[(f(0)+10*f(1/4)+f(2/4))/12
    (f(1/4)+10*f(2/4)+f(3/4))/12
    (f(2/4)+10*f(3/4)+f(1))/12+21/20*(1/h^2-1/12)];
for i=1/4:1/4:3/4
    x(i*4)=1/20*i^5+i^3;
end
xx=A^-1*b;
error=abs(xx'-x)