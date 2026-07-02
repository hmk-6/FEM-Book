# 2-4 平衡方程组求解程序设计作业（MATLAB R2024a）

学号：20252210045  
姓名：韩梦坤

## 1. 运行方法

1. 将整个文件夹保持原有目录结构解压。
2. 在 MATLAB R2024a 中把“当前文件夹”切换到本项目根目录。
3. 运行：

```matlab
main_2_4
```

4. 程序会自动运行全部算例，并生成：
   - `results/2-4_all_results.txt`
   - 各算例的 `.mat` / `.csv` 结果
   - `figures/poisson_solution_50.png`
   - `figures/poisson_error_50.png`
   - `figures/poisson_convergence.png`

Poisson 网格 `nx=ny=200` 可能需要数十秒至数分钟，取决于电脑性能。

## 2. 程序结构

- `main_2_4.m`：总入口。
- `src/ldlt_factor.m`：自行实现稠密 LDL^T 分解，检测零主元与非正主元。
- `src/ldlt_solve.m`：前代、对角求解、回代，支持多个右端项。
- `src/solve_equilibrium.m`：统一求解接口。
- `src/dense_to_skyline.m`、`skyline_ldlt_factor.m`、`skyline_ldlt_solve.m`：活动列/轮廓存储与求解。
- `src/solve_sparse_system.m`：稀疏求解器自动检测与后备接口。
- `src/assemble_poisson_t3.m`：二维 Poisson 方程 T3 有限元装配。
- `cases/`：全部验证算例。
- `input/`：JSON 示例输入文件。
- `report/`：作业报告 DOCX 与 PDF。

## 3. 验证算例

程序包括：

1. 2-3 一维两单元杆结构接口验证。
2. 2-3 二维两杆桁架接口验证。
3. `n=10,100,500,1000` 三对角 SPD 矩阵。
4. 非正定矩阵检测。
5. 病态矩阵双精度、4 位有效数字、单精度和右端项扰动分析。
6. 多载荷工况效率比较。
7. Skyline/活动列存储验证。
8. T3 有限元 Poisson 方程：`nx=ny=50,100,200`。

## 4. 稀疏求解器说明

`solve_sparse_system.m` 会按以下顺序工作：

1. 检测常见 PARDISO Project MATLAB 接口：
   `pardisoinit`、`pardisoreorder`、`pardisofactor`、`pardisosolve`、`pardisofree`。
2. 若接口可用，则调用 PARDISO。
3. 若未安装或调用失败，则使用：

```matlab
p = symamd(K);
xp = K(p,p) \ b(p);
x(p) = xp;
```

即 `MATLAB sparse backslash + SYMAMD` 后备方案。运行结果中会明确输出实际求解器名称，不会把后备方案冒充为 PARDISO。

若课程硬性要求提交 PARDISO/MUMPS 的真实运行截图，需要先在本机安装相应 MATLAB 接口，再重新运行主程序。

## 5. 2-3 与 2-4 的衔接

2-3 作业负责：

- 节点与单元输入；
- 单元刚度矩阵；
- 对号矩阵 LM；
- 总体刚度矩阵 K；
- 边界条件分块与缩减方程。

2-4 作业替换并强化“方程求解模块”：

```matlab
[dF, info] = solve_equilibrium(K_FF, rhs, 'ldlt');
```

求得自由位移后，继续复用 2-3 的完整位移重构、约束反力和单元应力/轴力后处理。

## 6. 提交前需要填写

- 报告封面姓名。
- `environment.txt` 中的 CPU、内存、操作系统和线程数。
- 在本人电脑运行后，将 MATLAB 实测时间更新到报告中。
- 若安装 PARDISO/MUMPS，加入安装方式、版本号与运行截图。

## 7. 数组下标

MATLAB 数组下标从 1 开始；程序中的节点号、单元号和自由度号均采用 1 起始编号。
