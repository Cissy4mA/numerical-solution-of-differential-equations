%%向前差分格式
%网格划分数量
m = 80;%空间节点数（x方向）
n = 3200;% 时间节点数（t方向）
h = 1/m;%空间步长
q = 1/n;%时间步长
a = 1;%扩散系数
%r=a*q/(h^2);

%matrix_R
R_=(1-2*r)*diag(ones(1,m-1));%主对角线元素（中心差分系数）
r1=(r)*diag(ones(1,m-2),1);%上对角线元素
r2=(r)*diag(ones(1,m-2),-1);%下对角线元素
R=r1+R_+r2;

% 构造块三对角系统矩阵
main_diag = kron(speye(n), speye(m-1));
lower_diag = kron(spdiags(ones(n,1), -1, n, n), -R);
RI = main_diag + lower_diag;

% 右端项（含初始条件）
F = zeros((m-1)*n, 1);
for j = 1:n
    t = q*j;
    f = q * sin(t);%右端项乘以Δt
    F((j-1)*(m-1)+1: j*(m-1)) = f*ones(m-1,1);%均匀分布源项
end

%初始条件（t=0时u=cos(πx)）
u0 = cos(pi*(h:h:1-h)');%空间离散初始值
F(1:m-1) = F(1:m-1) + R * u0;%第一层受初始条件影响

%求逆计算u
u = RI \ F;
%转换为U矩阵
U = reshape(u, m-1, n);
U = [zeros(1,n); U; zeros(1,n)];%添加边界
U=U(2:m,:);

%数值解与精确解的对比
figure;
%数值解
x=[1:(m-1)*n];
u_=zeros(1,(m-1)*n);
U1=U';
for i=1:m-1
    u_(n*(i-1)+1:n*i)=U1(:,i);
end
plot(x,u_,'-b','MarkerSize',5,'LineWidth',1);
hold on;
%精确解
x = linspace(0, 1, m-1);
y = linspace(0, 1, n);
[X, Y] = meshgrid(x, y);
Z=exp(-2*(pi^2)*Y).*cos(pi*X)+1-cos(Y);
z=zeros(1,(m-1)*n);
for i=1:m-1
    z(n*(i-1)+1:n*i)=Z(:,i);
end
x=[1:(m-1)*n];
plot(x,z,'-r','MarkerSize',5,'LineWidth',1);
xlabel('节点编号');
ylabel('函数值');
title('数值解与精确解的对比');
legend('数值解','精确解');

%误差
figure;
error=abs(z-u_);
plot(x,error,'-k','MarkerSize',5,'LineWidth',1);
xlabel('节点编号');
ylabel('误差绝对值');

%向前差分格式三维图像
x = linspace(0, 1, m-1);
t = linspace(0, 1, n);
[X, T] = meshgrid(x, t);%网格坐标
figure;
surf(X, T, U', 'EdgeColor','none','FaceColor','interp');
alpha(0.6);
xlabel('x'); ylabel('t'); zlabel('u(x,t)');
title('向前差分格式数值解');
view(30, 45);
colorbar;