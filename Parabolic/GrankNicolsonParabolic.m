function [u, x, t] = GrankNicolsonParabolic( fun, g, f1, f2, a, T, m, tau)
% Solve u\_t = a*u\_xx + fun(x,t)
%   over the region D = [0, 1, 0, T] = {(x,t) |0 <= x <= 1, 0 <= t <= T}
%   with the initial and boundary Conditions:
%          u(x,0) = g(x)
%          u(0,t) = f1(t), u(1,t) = f2(t)

% 生成格点
t = 0:tau:T;
n = length(t)-1;
hx = 1/m;
r = a*tau/hx/hx;
x = linspace(0,1,m+1);
u = zeros(m+1,n+1);
f = zeros(m+1,n+1);
% 设定好边界和初值条件
u(:,1) = g(x'); % 列向量
u(1,:) = f1(t); u(m+1,:) = f2(t); % 行向量
f = fun( meshgrid(x,t) ); %右端项初始化
f = f';

for k = 1:n
    d = (u(1:m-1,k)+u(3:m+1,k) )*r/2 +u(2:m,k)*(1-r)...
        +( f(2:m,k)+f(2:m,k+1) )*tau/2;
    d(1) = d(1) + u(1,k+1)*r/2;
    d(m-1) = d(m-1)+u(m+1,k+1)*r/2;
    % 调用追赶法代码求线性方程组
    u(2:m,k+1) = Chasing3(-r/2,1+r,-r/2,d);
end

end
