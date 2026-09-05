# Numerical Solution of Differential Equations — MATLAB Coursework

> 本科课程《微分方程数值解》的 MATLAB 实现。涵盖**常微分方程初值问题、抛物型、椭圆型、双曲型**四类方程的数值解法，包含不同步长、不同格式下的误差对比与收敛性分析。

## 内容概览

| 类别 | 典型问题 | 数值方法 | 目录 |
|------|----------|----------|------|
| 常微分方程初值问题 | $y' = f(t,y),\ y(t_0)=y_0$ | 二阶显式 / 隐式 Adams、梯形法、Euler；步长与格式误差对比 | `ODE/` |
| 抛物型 PDE | $u_t = a\,u_{xx} + f(x,t)$ | 向前差分（显式）、Crank–Nicolson（调用追赶法解三对角）、迎风格式 | `Parabolic/` |
| 椭圆型 PDE | $u_{xx} + u_{yy} = f(x,y)$（Poisson） | 五点差分格式（DCD 矩阵 + Kronecker 积），步长 $h=1/64$ 与 $1/128$ 对比 | `Elliptic/` |
| 双曲型 PDE | $u_{tt} = a^2 u_{xx}$（波动方程） | 显式 / 隐式格式、迎风格式，三对角 / 块三对角方程组求解 | `Hyperbolic/` |

## 目录结构

```
numerical-solution-of-differential-equations/
├── ODE/            # 常微分方程初值问题
├── Parabolic/      # 抛物型 PDE
├── Elliptic/       # 椭圆型 PDE
├── Hyperbolic/     # 双曲型 PDE
├── report/         # 课程设计报告（docx，已匿名）
├── README.md
└── LICENSE
```

## 文件说明

### `ODE/` — 常微分方程初值问题
- `Adams_2.m`：二阶显式 / 隐式 Adams 与符号精确解对比（问题 $y'=2y-3t^2,\ y(0)=2$）。
- `h1h2Adams.m`：两种步长 $h_1/h_2$ 下二阶显式 Adams 的误差对比。
- `ladder_adams.m`：梯形法（隐式单步 Adams 类）实现与对比。

### `Parabolic/` — 抛物型 PDE
- `GrankNicolsonParabolic.m`：Crank–Nicolson 格式函数，内部调用 `Chasing3` 求解三对角方程组。
- `Chasing3.m`：追赶法（Thomas 算法）求解三对角线性方程组。
- `paowu.m`：抛物型方程向前差分（显式）格式。
- `paowu_GN.m`：抛物型方程 Crank–Nicolson 格式。
- `paowu_xiangqian.m`：抛物型方程迎风（向前）格式。
- `testParabolic.m`：抛物型集成测试 / 绘图脚本。

### `Elliptic/` — 椭圆型 PDE
- `tuoyuan.m`：椭圆型 Poisson 方程五点差分主程序。
- `tuoyuan_5point.m`：五点差分格式，步长 $h=1/64$ 与 $1/128$ 的数值解、误差及三维曲面对比。

### `Hyperbolic/` — 双曲型 PDE
- `Chasing.m` / `Chasing3.m`：追赶法解三对角 / 块三对角方程组（双曲隐式格式用）。
- `Hyperbolic2.m`：双曲型另一种格式实现。
- `shuangqu.m`：双曲型波动方程数值求解。
- `shuangqu_xian.m`：双曲型显式 / 隐式格式（初始条件 $u(x,0)=2\sin(\pi x)$，精确解 $\sin\pi(x-t)+\sin\pi(x+t)$）。
- `testHyperbolic.m` / `testHyperbolic2.m`：双曲型集成测试 / 绘图脚本。

## 运行环境

- 需要 **MATLAB**（代码使用基础数值计算与绘图，无额外工具箱依赖）。
- 脚本类 `.m` 直接在命令行运行即可；函数类（如 `GrankNicolsonParabolic`）按帮助注释传入参数调用；`test*.m` 为集成测试与绘图入口。

## 课程设计报告

完整报告见 [`report/differential-equations-numerical-report.docx`](report/differential-equations-numerical-report.docx)。
GitHub 不直接在线渲染 `.docx`，请下载后查看。

## License

MIT — 见 [LICENSE](LICENSE)。
