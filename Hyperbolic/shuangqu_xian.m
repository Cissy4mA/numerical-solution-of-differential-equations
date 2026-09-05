%网格划分数量
m = 10;%空间节点数（x方向）
n = 40;% 时间节点数（t方向）
h = 1/m;%空间步长
q = 1/n;%时间步长
a = 1;%扩散系数
r=a*q/h;

%matrix_R
R_=(2*(1-r^2))*diag(ones(1,m-1));%主对角线元素（中心差分系数）
r1=(r^2)*diag(ones(1,m-2),1);%上对角线元素
r2=(r^2)*diag(ones(1,m-2),-1);%下对角线元素
R=r1+R_+r2;

%构造IRI矩阵
main_diag=kron(speye(n-1),-R);%主对角块
upper_diag=kron(spdiags(ones(n-1,1),1,n-1,n-1),speye(m-1));%上对角块
lower_diag=kron(spdiags(ones(n-1,1),-1,n-1,n-1),speye(m-1));%下对角块
IRI= main_diag + upper_diag + lower_diag;

% 右端项（含初始条件）
F = zeros((m-1)*(n-1), 1);
%初始条件（t=0时u=2sin(πx)）
u0 = 2*sin(pi*(h:h:1-h)');%空间离散初始值
F(1:m-1) = F(1:m-1)-u0;%第一层受初始条件影响
F((n-2)*(m-1)+1:(n-1)*(m-1)) = F((n-2)*(m-1)+1:(n-1)*(m-1))-u0;

%求逆计算u
u = inv(IRI)*F;
%转换为U矩阵
U = reshape(u, m-1, n-1);U=U';
U=padarray(U,[1 1],0,'both');%外面加一圈0
U(1,2:m)=u0';
U(n+1,2:m)=u0';
U(1:n/2,:)=abs(U(1:n/2,:));
U(n/2+1:n+1,:)=-abs(U(n/2+1:n+1,:));

%数值解与精确解的对比
figure;
%数值解
x=[1:(m+1)*(n+1)];
u_=zeros(1,(m+1)*(n+1));
U1=u2';
for i=1:m+1
    u_((n+1)*(i-1)+1:(n+1)*i)=U1(:,i);
end
plot(x,u_,'-b','MarkerSize',5,'LineWidth',1);
hold on;
%精确解
x = linspace(0, 1, m+1);
y = linspace(0, 2, n+1);
[X, Y] = meshgrid(x, y);
Z=sin(pi*(X-Y))+sin(pi*(X+Y));
z=zeros(1,(m+1)*(n+1));
for i=1:m+1
    z((n+1)*(i-1)+1:(n+1)*i)=Z(:,i);
end
x=[1:(m+1)*(n+1)];
plot(x,z,'--r','MarkerSize',5,'LineWidth',1);
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

%显式格式三维图像
x = linspace(0, 1, 11);
t = linspace(0, 1, 41);
[X, T] = meshgrid(x, t);%网格坐标
figure;
surf(X, T, u2', 'EdgeColor','k','FaceColor','interp');
alpha(0.6);
xlabel('x'); ylabel('t'); zlabel('u(x,t)');
title('隐式格式数值解');
view(75, 25);
colorbar;