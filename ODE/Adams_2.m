h=0.1;
%显式Adams
y=2;Y1=[2,2.2];%y1=y0+hy0=2.2
for i=2:20
    y=Y1(i)+0.5*h*(3*(2*Y1(i)-3*(0.1*(i-1)).^2)-2*Y1(i-1)+3*(0.1*(i-2)).^2);
    Y1=[Y1,y];
end
x=[0:0.1:2];
plot(x,Y1,'ro','MarkerFaceColor','r');
hold on;
%隐式Adams
y=2;Y2=[2];
for i=1:20
    y=((1+0.1)*y-1.5*h*((0.1*(i-1)).^2+(0.1*i).^2))/(1-h);
    Y2=[Y2,y];
end
plot(x,Y2,'bx','MarkerSize',5,'LineWidth',3);
hold on;
%精确解
syms y(t) t
dy=diff(y,t);
ode=dy==2*y-3*t.^2;
conds=y(0)==2;
ysol=dsolve(ode,conds);
y_sol=matlabFunction(ysol);
t=linspace(0,2,1000);
y_vals=y_sol(t);
plot(t,y_vals,'k','LineWidth',1);
legend('显式Adams','隐式Adams','精确解');
xlabel('x');
ylabel('y');
title('Adams数值解与精确解对比');
xticks(0:0.1:2);

Y3=[];
for i=0:0.1:2
    y=y_sol(i);
    Y3=[Y3,y];
end
Y=[Y1',Y2',Y3'];