%梯形
h=0.1;
%y_n+1=y_n+h/2(f(t_n,y_n)+f(t_n+1,y_n+1))
%f(t)=4t(y).^(1/2)
tic;
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
toc;

tic;
%显式Adams
y=1;Y2=[1,1.1];%y1=y0+hy0=1.1
for i=2:20
    y=Y2(i)+0.5*h*(3*(4*(h*i)*(Y2(i)).^(1/2))-4*(h*(i-1))*(Y2(i-1)).^(1/2));
    Y2=[Y2,y];
end
x=[0:0.1:2];
plot(x,Y2,'-bx','MarkerSize',5,'LineWidth',1);
hold on;
toc;

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