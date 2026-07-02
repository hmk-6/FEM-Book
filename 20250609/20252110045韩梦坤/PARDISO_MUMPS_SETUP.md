# PARDISO / MUMPS 外部求解器实验说明

课程要求大规模算例至少调用一种 MKL/PARDISO/MUMPS 外部稀疏直接求解器。
本项目的 `src/solve_sparse_system.m` 已预留并自动检测 PARDISO Project 常见 MATLAB 接口：

- `pardisoinit`
- `pardisoreorder`
- `pardisofactor`
- `pardisosolve`
- `pardisofree`

安装接口并加入 MATLAB path 后，重新运行 `main_2_4.m`，输出中的求解器名称应变为
`PARDISO Project MATLAB interface`。若接口不存在或调用失败，程序会明确显示
`MATLAB sparse backslash + SYMAMD (fallback)`。

## 建议记录内容

1. 求解器名称与版本。
2. 安装或编译方式。
3. CPU、内存、操作系统、MATLAB版本和线程数。
4. 矩阵阶数、非零元数、稀疏格式。
5. 符号分解、数值分解与回代总时间。
6. 相对残差。
7. 与MATLAB后备求解的时间、内存和误差对比。

## Intel oneMKL PARDISO关键概念

- 实对称正定矩阵类型通常对应 `mtype=2`。
- 输入采用CSR结构：数值数组、列索引数组和行指针数组。
- 分析/重排序、数值分解和求解可分阶段执行。
- 多个右端项应复用分析和数值分解结果。
- 重排序用于减少填充元和因子内存。

注意：不同MATLAB接口的函数参数可能存在差异，应以所安装接口自带示例为准。
