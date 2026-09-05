function [u, x, t] = Hyperbolic2( fun, phi, psi, alpha, beta, a, T, m, tau, method)
% Solve u\_tt = a*u\_xx + fun(x,t)
%   over the region D = [0, 1, 0, T] = {(x,t) |0 <= x <= 1, 0 <= t <= T}
%   with the initial and boundary Conditions:
%          u(x,0) = phi(x),  u\_t(x,0)=psi(x),
%          u(0,t) = alpha(t), u(1,t) = beta(t)
%  method 可选包括
%     Explicit: 显格式
%     Implicit: 隐格式 

% 生成格点
t = 0:tau:T;
n = length(t)-1;
hx = 1/m;
r = a*tau/hx;

x = linspace(0,1,m+1);
u = zeros(m+1,n+1);
f = zeros(m+1,n+1);
% 设定好边界和初值条件
u(:,1) = phi(x'); % 列向量
ut = psi(x');
u(1,:) = alpha(t); u(m+1,:) = beta(t); % 行向量
f = fun( meshgrid(x,t) ); %右端项初始化
f = f';

switch method
    case 'Explicit'
        if(r>1)
            disp('网格比过大，算法不稳定');
        end
        u(2:m,2) = tau* ut(2:m) + tau^2/2*f(2:m,1)...
            + r^2/2*(u(1:m-1,1) + u(3:m+1,1) )...
            + (1-r^2)*u(2:m,1); 
        c = [r^2, 2*(1-r^2), r^2];
        for k = 2:n
            u(2:m,k+1) = tau^2*f(2:m,k) - u(2:m,k-1)...
                +r^2* ( u(1:m-1,k) + u(3:m+1,k) )...
                +2*(1-r^2) *u(2:m,k);
        end
    case 'Implicit'
        d = 2*tau*ut(2:m)+tau^2*f(2:m,1)...
            +r^2/2*( u(1:m-1,1) +u(3:m+1,1) )...
            +(2-r^2)* u(2:m,1);
        d(1) = d(1)+u(1,2)*r^2/2;
        d(m-1) = d(m-1) + u(m,2)*r^2/2;
        u(2:m,2) = Chasing(2+r^2,-r^2/2,d);
        for k = 2:n
            d = 2*u(2:m,k)+tau^2*f(2:m,k)...
                +r^2/2*( u(1:m-1,k-1)+ u(3:m+1,k-1) )...
                -(1+r^2)*u(2:m,k-1);
            d(1) = d(1) + r^2/2*u(1,k+1);
            d(m-1) = d(m-1)+r^2/2*u(m+1,k+1);
            % 调用追赶法代码求线性方程组
            u(2:m,k+1) = Chasing(1+r^2, -r^2/2, d);
        end
    otherwise
        disp('选择方法错误');
end

end