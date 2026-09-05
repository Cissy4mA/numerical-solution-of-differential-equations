function u = Chasing3( a, b, c, f )
%追赶法求线性方程组，系数矩阵为三对角矩阵
%   [  b,  c,  0, ..., 0,  0, 0 ]
%   [  a,  b,  c, ..., 0,  0, 0 ]
%   [  ....................]
%   [  0,  0,  0, ..., a,  b,  c]
%   [  0,  0,  0, ..., 0,  a,  b]

n = length(f); % 矩阵维数
d = zeros(n-1,1)+b;   % n 维行向量，存储对角线元素
u = zeros(n,1);   % 解向量，不覆盖原右端项
y = f;  % 用于更新右端项

for i=2:n
    l = a/d(i-1);
    d(i) = b-l*c;
    y(i)=y(i)-l*y(i-1);
end

% 等价方程组的求解
u(n) = y(n)/d(n);
for i = n-1:-1:1
    u(i) = (y(i)-c*u(i+1) ) / d(i);
end

end