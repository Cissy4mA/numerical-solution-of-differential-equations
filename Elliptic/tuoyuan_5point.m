%网格划分数量
n=64;
h1=1/n; h2=1/n;

%matrix_C
C_=(2/h1.^2+2/h2.^2)*diag(ones(1,n-1));%主对角线元素（中心差分系数）
c1=-(1/h1.^2)*diag(ones(1,n-2),1);%上对角线元素
c2=-(1/h1.^2)*diag(ones(1,n-2),-1);%下对角线元素
C=c1+C_+c2;
%matrix_D
D=-(1/h2.^2)*diag(ones(1,n-1));

%构造DCD矩阵
main_diag=kron(speye(n-1),C);%主对角块：C矩阵的Kronecker积
upper_diag=kron(spdiags(ones(n-1,1),1,n-1,n-1),D);%上对角块：D矩阵的Kronecker积
lower_diag=kron(spdiags(ones(n-1,1),-1,n-1,n-1),D);%下对角块：D矩阵的Kronecker积
DCD= main_diag + upper_diag + lower_diag;

%f向量
F=[];
for j=1:n-1
    for i=1:n-1
        x=h1*i;
        y=h2*j;
        %poisson方程右端项
        f=2*pi.^2*exp(pi*(x+y)).*(sin(pi*x).*cos(pi*y)+cos(pi*x).*sin(pi*y));
        F=[F,f];
    end
end

%求逆计算u
u=-inv(DCD)*F';u=u';
%转换为U矩阵
U=zeros(n-1,n-1);
for i=1:n-1
    U(:,i)=u((n-1)*(i-1)+1:(n-1)*i)';
end
U=padarray(U,[1 1],0,'both');%外面加一圈0

%数值解与精确解的对比
figure;
%数值解
x=[1:n*n];
u_=zeros(1,n*n);
for i=1:n
    u_(n*(i-1)+1:n*i)=U(1:n,i);
end
plot(x,u_,'-b','MarkerSize',5,'LineWidth',1);
hold on;
%精确解
x = linspace(0, 1, n);
y = linspace(0, 1, n);
[X, Y] = meshgrid(x, y);
Z=exp(pi*(X+Y)).*sin(pi*X).*sin(pi*Y);
z=zeros(1,n*n);
for i=1:n
    z(n*(i-1)+1:n*i)=Z(:,i);
end
x=[1:n*n];
plot(x,z,'--r','MarkerSize',5,'LineWidth',1);
xlabel('节点编号');
ylabel('函数值');
title('h=1/128数值解与精确解的对比');
legend('数值解','精确解');

%误差
figure;
error=abs(z-u_);
xlabel('节点编号');
ylabel('误差绝对值');
plot(x,error,'-k','MarkerSize',5,'LineWidth',1);

%五点差分格式三维图像
[X, Y]=meshgrid(0:h1:1,0:h2:1);%生成网格坐标
figure;
mesh(X, Y, U, 'EdgeColor','k','FaceColor','interp');
alpha(0.7);
hold on;
xlabel('x');
ylabel('y');
zlabel('u(x,y)');
colorbar;
title('五点差分格式 h=1/64');
view(30, 45);%调整视角（方位角45°，俯仰角30°）
grid on;