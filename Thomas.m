function u =Thomas(A, B, Y, d)% 解三对角矩阵
u 三对角线性方程组的解
% A 对角线下vector
% B 对角线vector
% Y 对角线上vector
% d right side vector
N = length(B);
u = zeros(N,1);
g(1) = d(1)/B(1);
w(1) = Y(1)/B(1);
for i = 2:N-1
    g(i) = (d(i)-A(i-1)*g(i-1))/(B(i)-A(i-1)*w(i-1));
    w(i) = Y(i)/(B(i)-A(i-1)*w(i-1));
end
g(N) = (d(N)-A(N-1)*g(N-1))/(B(N)-A(N-1)*w(N-1));
u(N) = g(N);
for i = N-1:-1:1
    u(i) = g(i)-w(i)*u(i+1);
end
