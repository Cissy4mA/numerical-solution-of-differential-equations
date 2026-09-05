% 定义 x 和 y 的范围（根据需求调整范围）
x = linspace(0, 1, 11);  % 示例范围：0 ≤ x ≤ 1
y = linspace(0, 2, 41);  % 示例范围：0 ≤ y ≤ 1
[X, Y] = meshgrid(x, y);  % 生成网格点

% 计算函数值 u(x, y)
Z = sin(pi*(X-Y))+sin(pi*(X+Y));

% 绘制三维曲面图
figure;
surf(X, Y, Z);
alpha(0.5);
title('精确解：u(x,y)=sin(pi(x-t))+sin(pi(x+t))');
xlabel('x');
ylabel('t');
zlabel('u(x,t)');
colorbar;% 显示颜色条
%shading interp;  % 平滑颜色过渡
%shading flat;
colormap jet;

% 调整视角（可选）
view(75, 25);     % 方位角30°，俯仰角45°