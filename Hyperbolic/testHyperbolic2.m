% testHyperbolic2.m 
% 用隐式迎风格式求解对流占优扩散方程
% 本段代码直接将所有求解写在同一文件中

clear; clc; close all
% 方程阐述设置
a = 1; b = 1; u0 =2;
% 精确解
sol = @(x,t)[sin(pi(x-t))+sin(pi(x+t))];

% 数值解参数设置
h = 0.5; tau = 0.02;
L = 50; % 空间截断
T = 50; % 时间截断

% 生成网格
x = 0:h:L;   m = length(x)-1;
t = 0:tau:T;  n =length(t)-1;

% 解函数初始化，x=L 的地方看成 \infty 处
u = zeros(m+1,n+1);
u(1,:) = u0+t*0; 

r1 = b*tau/h; r2 = a*tau/(h^2);
w1 = -r1-r2; w2 =1+r1+2*r2; w3 = -r2;
for k = 1:n
    d = u(2:m,k);
    d(1) = d(1)-w1*u(1,k+1);
    d(m-1) = d(m-1)-w3*u(m+1,k+1);
    u(2:m,k+1) = Chasing3(w1,w2,w3,d);
end

[xx,tt] = meshgrid(x,t);
u_exact = sol(xx',tt'); % 精确解
mesh(xx',tt',u)%显式解画图
xlabel('x'); ylabel('t'); zlabel('u');

figure();
k = 25/tau+1;
plot(x,u_exact(:,k),'-b',x,u(:,k),'--r')
xlabel('x'); ylabel('u');
% title('t=25 时解的比较')
legend('精确解','数值解');
grid on

figure()
k=25/h+1;
plot(t,u_exact(k,:),'-b',t,u(k,:),'--r')
xlabel('t'); ylabel('u');
% title('x=25 时解的比较')
legend('精确解','数值解');
grid on

figure();
n = length(t);
M = moviein(n);% 建立一个 n 列的大矩阵
for k = 1:n
    plot(x,u_exact(:,k),'-b');
    hold on;
    axis([0,L,0,2]);
    grid on
    hold on
    plot(x,u(:,k),'--r');% 显式解比较
    hold off
    legend('精确解','数值解');
    M(k) = getframe; % 将图形保存到M矩阵
end
movie(M,1);  % 播放画面1次
