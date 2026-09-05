%显式Adams
tic;
h=0.1;
y=1;Y2=[1,1.1];%y1=y0+hy0=1.1
for i=2:20
    y=Y2(i)+0.5*h*(3*(4*(h*i)*(Y2(i)).^(1/2))-4*(h*(i-1))*(Y2(i-1)).^(1/2));
    Y2=[Y2,y];
end
x=[0:0.1:2];
plot(x,Y2,'-bx','MarkerSize',5,'LineWidth',1);
hold on;
toc;

tic;
h=0.05;
y=1;Y2=[1,1.1];%y1=y0+hy0=1.1
for i=2:40
    y=Y2(i)+0.5*h*(3*(4*(h*i)*(Y2(i)).^(1/2))-4*(h*(i-1))*(Y2(i-1)).^(1/2));
    Y2=[Y2,y];
end
x=[0:0.05:2];
plot(x,Y2,'-gd','MarkerSize',5,'LineWidth',1);
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
legend('h=0.1','h=0.05','精确解');
xlabel('t');
ylabel('y');
title('不同h的二阶显式Adams对比');
xticks(0:0.1:2);