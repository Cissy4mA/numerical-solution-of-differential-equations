% testHyperbolic.m 双曲型方程有限差分法测试函数

clear; clc; close all
% 各个函数
sol = @(x,t)[sin(pi*(x-t))+sin(pi*(x+t))];% 真实解
fun = @(x,t)[x.*0];
phi = @(x)[2.*sin(pi.*x)];  psi = @(x)[(0*x)];
alpha = @(t)[t.*0];  beta = @(t)[t.*0];

a = 1; T=2;
m = 10; tau = 1/20;

% 用显式方法求解
[u1, x, t] = Hyperbolic2( fun, phi, psi, alpha, beta, a, T, m, tau, 'Explicit');
[u2, x, t] = Hyperbolic2( fun, phi, psi, alpha, beta, a, T, m, tau, 'Implicit');

[xx,tt] = meshgrid(x,t);
u_exact = sol(xx',tt'); % 精确解
mesh(xx',tt',u1)%显式解画图
xlabel('x'); ylabel('t'); zlabel('u');

figure();
n = length(t);
M = moviein(n);% 建立一个 n 列的大矩阵
for k = 1:n
    plot(x,u_exact(:,k),'-k','Linewidth',3);
    hold on;
    axis([0,1,-1.05,1.05]);
    grid on
    hold on
    plot(x,u1(:,k),'-r','Linewidth',2);% 显式解比较
    plot(x,u2(:,k),'--b');% 显式解比较
    legend('真实解','数值解-显式','数值解-隐式')
%     sprintf('t=t_%d',k)
    hold off
    M(k) = getframe; % 将图形保存到M矩阵
end
movie(M,1);  % 播放画面1次
