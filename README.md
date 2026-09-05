# Numerical Solution of Differential Equations — MATLAB

> Classical finite-difference and linear-multistep methods for ODE initial-value problems and second-order PDEs (parabolic, elliptic, hyperbolic), with step-size and scheme comparisons, error analysis, and convergence observations. The figures and conclusions from the original report are embedded below.

## Overview

| Category | Prototype | Numerical methods | Directory |
|----------|-----------|-------------------|-----------|
| ODE IVP | `y' = f(t, y), y(t0) = y0` | 2nd-order explicit Adams, trapezoidal (implicit Adams), Euler; step-size and scheme error comparison | `ODE/` |
| Parabolic PDE | `u_t = a u_xx + f(x,t)` | Forward (explicit), Crank–Nicolson, upwind; Thomas/chasing for tridiagonal systems | `Parabolic/` |
| Elliptic PDE | Poisson `u_xx + u_yy = f(x,y)` | Five-point stencil with Kronecker-product block matrix; grids `h = 1/64` and `1/128` | `Elliptic/` |
| Hyperbolic PDE | Wave equation `u_tt = a^2 u_xx` | Explicit, implicit, and upwind schemes; block-tridiagonal solver | `Hyperbolic/` |

## Directory Structure

```
numerical-solution-of-differential-equations/
├── ODE/            # ODE initial-value problems
├── Parabolic/      # Parabolic PDEs
├── Elliptic/       # Elliptic PDEs
├── Hyperbolic/     # Hyperbolic PDEs
├── images/         # figures extracted from the report
├── README.md
└── LICENSE
```

## Report — Problems, Methods, Results & Conclusions

实验一  常微分方程数值解法
实验试题
求解下列初值问题的数值解：
要求：选择不同的步长及不同的求解方法（收敛阶不同），对比不同步长下采用同一求解方法时产生的误差，以及同一步长下采用不同求解方法产生的误差。
试题分析及求解过程
不同步长下采用同一求解方法时产生的误差：采用二阶显式Adams法(1.2)，步长取h=0.1和h=0.05，由实验试题可知，。使用Euler法算出，再通过迭代依次算出到
接着通过符号计算求解了常微分方程的精确解，将其转换为可调用的函数，并在时间区间内生成1000个等间距点计算对应的数值解，最终绘制精确解曲线（黑色实线），同时添加图例、坐标轴标签和标题。
```matlab
h=0.1;
y=1;Y2=[1,1.1];%y1=y0+hy0=1.1
for i=2:20
y=Y2(i)+0.5*h*(3*(4*(h*i)*(Y2(i)).^(1/2))-4*(h*(i-1))*(Y2(i-1)).^(1/2));
Y2=[Y2,y];
end
x=[0:0.1:2];
plot(x,Y2,'-bx','MarkerSize',5,'LineWidth',1);
hold on;
h=0.05;
y=1;Y2=[1,1.1];%y1=y0+hy0=1.1
for i=2:40
y=Y2(i)+0.5*h*(3*(4*(h*i)*(Y2(i)).^(1/2))-4*(h*(i-1))*(Y2(i-1)).^(1/2));
Y2=[Y2,y];
end
x=[0:0.05:2];
plot(x,Y2,'-gd','MarkerSize',5,'LineWidth',1);
hold on;
%精确解
syms y(t) t
dy=diff(y,t);
ode=dy==4*t*y.^(1/2);
conds=y(0)==1;
ysol=dsolve(ode,conds);
y_sol=matlabFunction(ysol);
t=linspace(0,2,1000);
y_vals=y_sol(t);
plot(t,y_vals(2,:),'k','LineWidth',1);
legend('h=0.1','h=0.05','精确解');
xlabel('t');
ylabel('y');
title('不同h的二阶显式Adams对比');
xticks(0:0.1:2);
```
同一步长下采用不同求解方法产生的误差：采用二阶显式Adams法和梯形法(1.3)，步长取h=0.1。方法同上。
```matlab
h=0.1;
%y_n+1=y_n+h/2(f(t_n,y_n)+f(t_n+1,y_n+1))
%f(t)=4t(y).^(1/2)
y=1;t=0;%y(0)=1,t_0=0
Y1=[1];%储存y_n值
for i=2:21
syms x;
eq= x==Y1(i-1)+(h/2)*(4*t*(Y1(i-1)).^(1/2)+4*(t+h)*x.^(1/2));
y=solve(eq,x);
Y1=[Y1,y];
t=t+h;
end
x=[0:0.1:2];
plot(x,Y1,'-ro','MarkerFaceColor','r','LineWidth',1);
hold on;
%显式Adams
y=1;Y2=[1,1.1];%y1=y0+hy0=1.1
for i=2:20
y=Y2(i)+0.5*h*(3*(4*(h*i)*(Y2(i)).^(1/2))-4*(h*(i-1))*(Y2(i-1)).^(1/2));
Y2=[Y2,y];
end
x=[0:0.1:2];
plot(x,Y2,'-bx','MarkerSize',5,'LineWidth',1);
hold on;
%精确解
syms y(t) t
dy=diff(y,t);
ode=dy==4*t*y.^(1/2);
conds=y(0)==1;
ysol=dsolve(ode,conds);
y_sol=matlabFunction(ysol);
t=linspace(0,2,1000);
y_vals=y_sol(t);
plot(t,y_vals(2,:),'k','LineWidth',1);
legend('梯形','显式Adams','精确解');
xlabel('t');
ylabel('y');
title('h=0.1 梯形与二阶显式Adams对比');
xticks(0:0.1:2);
```
实验结果及分析
不同步长下采用同一求解方法时产生的误差：通过图1不难看出，步长为0.05比步长为0.1的显式Adams法求得的数值解更接近精确解曲线。所以可以得出结论步长越小，求得的数值解更精确。通过tictoc记录不同步长下的用时，得到步长为0.1时“历时 0.080657 秒”， 步长为0.05时“历时 0.078674 秒”，运行时长相差不大。

![图 11 不同步长同一求解方法](images/fig_001.jpg)
*图 11 不同步长同一求解方法*

同一步长下采用不同求解方法产生的误差：通过图2可知，同一步长下，梯形法几乎和精确解重合，而显式Adams法与精确解还有一定距离。通过tictoc记录两种方法下的用时，梯形法“历时 6.401789 秒”， 显式Adams法“历时 0.003981 秒”。由此可见，梯形法比显式Adams法精确度更高，但用时更久。

![图 12 同一步长不同求解方法](images/fig_002.jpg)
*图 12 同一步长不同求解方法*

总结
实验结果表明，梯形法虽然在精度上显著优于显式Adams法，但其计算成本也更高，这表明在实际应用中，需要根据问题的具体需求选择合适的数值方法。对于对精度要求极高的问题，梯形法是更好的选择；而对于对效率要求较高的问题，显式Adams法则更为合适。此外，步长的选择对数值解的精度有直接影响，步长越小，数值解越接近精确解，但步长的减小并不总是导致计算时间的增加，这可能与实验环境或具体实现有关。通过tictoc记录运行时间，验证了不同方法和步长下的计算效率，这为数值方法的选择提供了实际依据。总体而言，本次实验不仅加深了我对数值方法的理解，还提供了在实际应用中选择和优化数值方法的案例。
实验二  椭圆型方程的数值解法
实验试题
求下列边值问题的数值解：
其精确解为.
要求：取步长 , 作五点差分格式。
试题分析及求解过程
对于二维Poisson方程的Dirichlet边值问题
在内节点处用二阶中心差商，并引入差分算子
可用差分方程(2.4)求解
将（2.3）和（2.4）改写为未知量的线性方程组
记
则差分方程可以写成
进一步可以写成
由公式2.1可知，边界条件都为0，所以和都为0。编写程序时，先构造由C和D构成的分块系数矩阵，再通过网格节点数依次计算Poisson方程右端项的值，构造f向量。用矩阵求逆算出内一列节点u的值，重新排列u，并在周围加一圈外节点，获得网格函数的矩阵。使用meshgrid函数绘制三维图像。同时把网格矩阵降维成一维向量，和精确解绘制在同一张二维图像中，进行对比。
程序只需修改步长（n的倒数），便可获得不同步长下的五点差分格式。以下以步长为1/64作为示例。
```matlab
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
title('h=1/64数值解与精确解的对比');
legend('数值解','精确解');
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
```
最后，绘制精确解三维图像，为后续对比做准备。
```matlab
n=64;
x = linspace(0, 1, n);  % 示例范围：0 ≤ x ≤ 1
y = linspace(0, 1, n);  % 示例范围：0 ≤ y ≤ 1
[X, Y] = meshgrid(x, y);  % 生成网格点
% 计算函数值 u(x, y)
Z = exp(pi * (X + Y)) .* sin(pi * X) .* sin(pi * Y);
% 绘制三维曲面图
figure;
surf(X, Y, Z);
alpha(0.5);
title('精确解：u(x,y)=exp(pi*(X+Y)).*sin(pi*X).*sin(pi*Y)');
xlabel('x');ylabel('y');zlabel('u(x,y)');
colorbar;% 显示颜色条
%shading interp;  % 平滑颜色过渡
%shading flat;
colormap jet;
% 调整视角（可选）
view(30, 45);     % 方位角30°，俯仰角45°
```
实验结果及分析
3.1 步长=1/64的五点差分格式

![图 23精确解三维图像](images/fig_003.jpeg)
*图 23精确解三维图像*


![图 24 步长为1/64的五点差分格式三维图像](images/fig_004.jpeg)
*图 24 步长为1/64的五点差分格式三维图像*

图2-3和图2-4分别为该Poisson方程的精确解和步长为1/64的数值解，三维图像可以看清整体的趋势，但不便于比较精确解和数值解。下图将精确解和数值解降维到二维平面进行对比。

![图 25 步长为1/64的数值解与精确解对比](images/fig_005.jpeg)
*图 25 步长为1/64的数值解与精确解对比*

数值解总体近似与精确解，整体趋势一致。下面将局部放大观察，数值解在前半部分（3000节点编号）之前的波动幅度更小，后半部分波动更大，峰值更高。

![图 26 步长为1/64的数值解与精确解对比（局部）](images/fig_006.jpeg)
*图 26 步长为1/64的数值解与精确解对比（局部）*


![图 26 步长为1/64的数值解与精确解对比（局部）](images/fig_007.jpeg)
*图 26 步长为1/64的数值解与精确解对比（局部）*

3.2 步长=1/128的五点差分格式
首先绘制该补偿下的五点差分格式三维图像，观察到与精确解的图像近似，步长更小，网格更密集，说明该数值解计算正确。

![图 27步长为1/128的五点差分格式三维图像](images/fig_008.jpg)
*图 27步长为1/128的五点差分格式三维图像*


![图 28步长为1/128的数值解与精确解对比](images/fig_009.jpeg)
*图 28步长为1/128的数值解与精确解对比*


![图 29 步长为1/128的数值解与精确解对比（局部）](images/fig_010.jpg)
*图 29 步长为1/128的数值解与精确解对比（局部）*


![图 29 步长为1/128的数值解与精确解对比（局部）](images/fig_011.jpg)
*图 29 步长为1/128的数值解与精确解对比（局部）*

总结
实验三  抛物型方程的数值解法
实验试题
分别用向前差分格式和Grank-Nicolson格式求下列一维抛物方程的初边值问题的数值解：
其精确解为.
要求：网格分别取两组以及 比较不同数值方法之间的差异。
试题分析及求解过程
2.1 向前差分格式
对于一维抛物方程的初边值问题
将求解区域
作剖分。将区间作等分，将区间作n等分，记
分别称和为空间步长和时间步长。
将（3.3）作向前差分格式：
其中
令，将（3.5）改写为线性方程组
记
则（3.6）可以写成
其中表示单位矩阵，进一步写成矩阵形式：
先构造矩阵，接着构造由和单位阵组成的分块系数矩阵。根据（3.1）可知，，因此由构成的方程右端项矩阵最后还需添加初始条件。最后，矩阵求逆算出节点、精确解与数值解对比、网格函数可视化这部分思路与实验二类似，可以沿用。
以下以作为示例。
```matlab
m = 40;%空间节点数（x方向）
n = 3200;% 时间节点数（t方向）
h = 1/m;%空间步长
q = 1/n;%时间步长
a = 1;%扩散系数
r=a*q/(h^2);
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
```
绘制精确解的三维图像，以验证数值解的准确性：
```matlab
n=100;
x = linspace(0, 1, n);  % 示例范围：0 ≤ x ≤ 1
y = linspace(0, 1, n);  % 示例范围：0 ≤ y ≤ 1
[X, Y] = meshgrid(x, y);  % 生成网格点
% 计算函数值 u(x, y)
Z = exp(-2*(pi^2)*Y).*cos(pi*X)+1-cos(Y);
% 绘制三维曲面图
figure;
surf(X, Y, Z);
alpha(0.5);
title('精确解：u(x,y)=exp(-2*(pi^2)*t).*cos(pi*x)+1-cos(t)');
xlabel('x');
ylabel('t');
zlabel('u(x,t)');
colorbar;% 显示颜色条
%shading interp;  % 平滑颜色过渡
%shading flat;
colormap jet;
view(30, 45);     % 方位角30°，俯仰角45°
```
Grank-Nicolson格式
将（3.2）作Grank-Nicolson格式：
其中
令，将（3.11）改写为线性方程组
记
则（3.13）可以写成
其中表示单位矩阵，进一步写成矩阵形式：
先构造和矩阵，接着构造由和组成的分块系数矩阵。根据（3.1）可知，，因此由构成的方程右端项矩阵最后添加初始条件。最后，矩阵求逆算出节点、精确解与数值解对比、网格函数可视化。其中添加了可视化精确解与数值解的误差的绝对值。
以下以作为示例。
```matlab
%网格划分数量
m = 40;%空间节点数（x方向）
n = 3200;% 时间节点数（t方向）
h = 1/m;%空间步长
q = 1/n;%时间步长
a = 1;%扩散系数
r=a*q/(h^2);
%matrix_R1
R_=(1+r)*diag(ones(1,m-1));%主对角线元素（中心差分系数）
r1=(-r/2)*diag(ones(1,m-2),1);%上对角线元素
r2=(-r/2)*diag(ones(1,m-2),-1);%下对角线元素
R1=r1+R_+r2;
%matrix_R2
R_=(1-r)*diag(ones(1,m-1));%主对角线元素（中心差分系数）
r1=(r/2)*diag(ones(1,m-2),1);%上对角线元素
r2=(r/2)*diag(ones(1,m-2),-1);%下对角线元素
R2=r1+R_+r2;
% 构造块三对角系统矩阵
main_diag = kron(speye(n), R1);
lower_diag = kron(spdiags(ones(n,1), -1, n, n), -R2);
RI = main_diag + lower_diag;
% 右端项（含初始条件）
F = zeros((m-1)*n, 1);
for j = 1:n
t1 = q*j; t2=q*(j+1);
f = q/2 * (sin(t1)+sin(t2));%右端项乘以Δt
F((j-1)*(m-1)+1: j*(m-1)) = f*ones(m-1,1);
end
%初始条件（t=0时u=cos(πx)）
u0 = cos(pi*(h:h:1-h)');%空间离散初始值
F(1:m-1) = F(1:m-1) + R2*u0;%第一层受初始条件影响
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
title('Grank-Nicolson格式数值解');
view(30, 45);
colorbar;
```
实验结果及分析
3.1 向前差分格式
3.1.1 的向前差分格式

![图 3-1 精确解三维图像](images/fig_012.jpeg)
*图 3-1 精确解三维图像*


![图 3-2 时的向前差分格式](images/fig_013.jpeg)
*图 3-2 时的向前差分格式*

由图3-1和图3-2的可视化，该向前差分数值解与精确解趋势大致一致，说明数值解计算准确。

![图 3-3 时的数值解与精确解对比](images/fig_014.jpeg)
*图 3-3 时的数值解与精确解对比*

3.1.2 的向前差分格式

![图 3-4 时的向前差分格式](images/fig_015.jpeg)
*图 3-4 时的向前差分格式*


![图 3-5 时的向前差分格式](images/fig_016.jpeg)
*图 3-5 时的向前差分格式*

对比图3-4和图3-5，发现减小x的步长，得到的数值解图像变得更加圆滑，可能说明拟合效果更好。接下来用函数值对比观察数值解的误差。

![图 3-6 时的数值解与精确解对比](images/fig_017.jpeg)
*图 3-6 时的数值解与精确解对比*

3.2 Grank-Nicolson格式
3.2.1 的Grank-Nicolson格式

![图 3-7 精确解三维图像](images/fig_018.jpeg)
*图 3-7 精确解三维图像*


![图 3-8 的Grank-Nicolson格式](images/fig_019.jpeg)
*图 3-8 的Grank-Nicolson格式*

由图3-7和图3-8的可视化，该向前差分数值解与精确解趋势大致一致，说明数值解计算准确。

![图 3-9 时的数值解与精确解对比](images/fig_020.jpeg)
*图 3-9 时的数值解与精确解对比*


![图 3-10 时的误差绝对值](images/fig_021.jpeg)
*图 3-10 时的误差绝对值*

3.2.2 的Grank-Nicolson格式

![图 3-11 的Grank-Nicolson格式](images/fig_022.jpeg)
*图 3-11 的Grank-Nicolson格式*


![图 3-12 的Grank-Nicolson格式](images/fig_023.jpeg)
*图 3-12 的Grank-Nicolson格式*

步长变小对Grank-Nicolson格式的影响看上去不大，主要是因为该格式具有无条件稳定性和较高的精度。即使步长减小，数值解的精度提升有限，且整体趋势保持一致。

![图 3-9 时的数值解与精确解对比](images/fig_024.jpeg)
*图 3-9 时的数值解与精确解对比*


![图 3-10 时的误差绝对值](images/fig_025.jpeg)
*图 3-10 时的误差绝对值*

总结
通过本次实验，深刻体会到不同数值方法的特点和适用场景。向前差分格式简单易实现，但在稳定性和精度上受限较大；Grank-Nicolson格式虽实现稍复杂，但其高精度和无条件稳定性使其在求解抛物型方程时更具优势。
实验四  双曲型方程的数值解法
实验试题
分别用显式和隐式差分格式求下列波动方程混合边值问题的数值解：
其精确解为.
要求：
比较时在处两种方法得到的数值解；
比较 时在处两种方法得到的数值解。
试题分析及求解过程
2.1 显式格式
对于波动方程的初边值问题
将求解区域
作剖分。将区间作等分，将区间作n等分，记
分别称和为空间步长和时间步长。
将（4.3）作显式格式：
其中
令，将（4.5）改写为线性方程组
记
则（4.6）可以写成
其中表示单位矩阵，进一步写成矩阵形式：
2.2 隐式格式
将（4.3）作隐式格式：
其中当时为显格式.
考虑的情况，得到隐式差分格式
令，将（4.12）改写为线性方程组
记
则（3.13）可以写成
其中表示单位矩阵，进一步写成矩阵形式：
实验结果及分析
3.1的数值解
该双曲型方程的精确解如下图所示：

![图 4-1 精确解图像](images/fig_026.jpg)
*图 4-1 精确解图像*

显式和隐式数值解分别为图4-2和图4-3,数值解与精确解的对比分别为图4-4和图4-5：

![图 4-2时的显式数值解](images/fig_027.jpg)
*图 4-2时的显式数值解*


![图 4-3 时的隐式数值解](images/fig_028.jpg)
*图 4-3 时的隐式数值解*


![图 4-4显式数值解与精确解对比](images/fig_029.jpg)
*图 4-4显式数值解与精确解对比*


![图 4-5隐式数值解与精确解对比](images/fig_030.jpg)
*图 4-5隐式数值解与精确解对比*

图片难以观察出显式和隐式数值解的精确度和差别，接下来选取指定时间点的数值作对比，当时，分别对应k=5，10，15，20时的精确解和数值解：
表 4-1 ,时的精确解
分别对比显式数值解和精确解，隐式数值解和精确解，发现显式数值解更加接近精确解的值。再计算每个数值解与精确解的误差的绝对值的总和，以此验证两种差分格式的精确度。
显式数值解的误差总和：1.4766e-14，
隐式数值解的误差总和：2.4293，
显式差分格式的误差接近于0，总误差优于隐式差分格式。
3.2的数值解
此时的显式和隐式数值解分别为图4-6和图4-7，数值解与精确解的对比分别为图4-8和图4-9：

![图 4-6时的显式数值解](images/fig_031.jpeg)
*图 4-6时的显式数值解*


![图4-7 时的隐式数值解](images/fig_032.jpeg)
*图4-7 时的隐式数值解*


![图 4-8显式数值解与精确解对比](images/fig_033.jpeg)
*图 4-8显式数值解与精确解对比*


![图4-9隐式数值解与精确解对比](images/fig_034.jpeg)
*图4-9隐式数值解与精确解对比*

当时，分别对应k=10，20，30，40时的精确解和数值解：
表 4-4,时的精确解
表 4-5 ,时的显式数值解
表 4-6 ,时的隐式数值解
显式数值解的误差总和：0.2881，
隐式数值解的误差总和：0.8534，
因此，在这种步长情况下，显式总误差小于隐式总误差。
总结

## Running the Code

- Requires **MATLAB** with basic numerical and plotting support; no external toolboxes are needed.
- Each `.m` file is self-contained. Scripts that start with `test*.m` are entry points that run a complete method and produce comparison plots.
- For the Parabolic and Hyperbolic implicit schemes, make sure the corresponding `Chasing*.m` files are on the MATLAB path.

## License

MIT — see [LICENSE](LICENSE).
