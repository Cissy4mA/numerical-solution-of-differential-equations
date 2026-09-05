% testParabolic

clear; clc; close all
% 各个函数
sol = @(x,t)[exp(-2*(pi^2)*t).*cos(pi*x)+1-cos(t)];% 真实解
fun = @(x,t)[x.*0];
g = @(x)[cos(pi*x)];
f1 = @(t)[t.*0];
f2 = @(t)[t.*0];
a = 1; T=0.6;
m = 40; tau = 1/3200;

% 采用向前差分，注意网格比 r = a*tau/h^2<=1/2
% [u, x, t] = ForwardParabolic( fun, g, f1, f2, a, T, m, tau);
[u, x, t] = GrankNicolsonParabolic( fun, g, f1, f2, a, T, m, tau);
% [u, x, t] = GrankNicolsonParabolic( fun, g, f1, f2, a, T, m, tau);

% [u1, x, t] = Parabolic( fun, g, f1, f2, a, T, m, tau, 'Forward');
% [u2, x, t] = Parabolic( fun, g, f1, f2, a, T, m, tau, 'Backward');
% [u3, x, t] = Parabolic( fun, g, f1, f2, a, T, m, tau, 'GrankNicolson');

[xx,tt] = meshgrid(x,t);
u_exact = sol(xx',tt');
% u=u1;

mesh(xx',tt',u)
xlabel('x'); ylabel('t'); zlabel('u');
view(30, 45);

figure();
n = length(t);
M = moviein(n);% 建立一个 n 列的大矩阵
for k = 1:n
    plot(x,u_exact(:,k),'-k');
    hold on;
    axis([0,1,0,1.05]);
    grid on
    hold on
    plot(x,u(:,k),'ks','MarkerSize',6,...
        'MarkerEdgeColor','g',...
        'MarkerFaceColor','r');
    hold off
    legend('真实解','数值解')
    M(k) = getframe; % 将图形保存到M矩阵
end
movie(M,1);  % 播放画面1次
