% ============================================================
% 2-4 平衡方程组求解程序设计作业 - 主程序
% 功能：
% 1. 稠密矩阵 LDL^T 分解与求解
% 2. 2-3 一维杆与二维桁架接口验证
% 3. 三对角 SPD 矩阵规模测试
% 4. 非正定矩阵检测
% 5. 病态矩阵误差、残差与条件数分析
% 6. 多载荷工况“先分解、后多次回代”效率比较
% 7. 活动列/轮廓（Skyline）存储与求解
% 8. T3 有限元求解二维 Poisson 方程
% 9. 稀疏求解器自动检测：PARDISO Project / MATLAB sparse fallback
% ============================================================
clear; clc; close all;

projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(projectRoot, 'src')));
addpath(genpath(fullfile(projectRoot, 'cases')));

resultsDir = fullfile(projectRoot, 'results');
figuresDir = fullfile(projectRoot, 'figures');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end
if ~exist(figuresDir, 'dir'), mkdir(figuresDir); end

resultFile = fullfile(resultsDir, '2-4_all_results.txt');
if exist(resultFile, 'file'), delete(resultFile); end

diary(resultFile);
diary on;

fprintf('============================================================\n');

fprintf('============================================================\n\n');

try
    case_json_input(projectRoot,resultsDir);
    case0_bar(resultsDir);
    case0_truss(resultsDir);
    case1_tridiagonal(resultsDir);
    case2_indefinite(resultsDir);
    case3_ill_conditioned(resultsDir);
    case_multiload(resultsDir);
    case_skyline(resultsDir);
    case4_poisson(resultsDir, figuresDir, [50, 100, 200]);

    fprintf('\n============================================================\n');
    fprintf('全部算例运行完成。\n');
    fprintf('结束时间：%s\n', char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
    fprintf('总输出文件：%s\n', resultFile);
    fprintf('============================================================\n');
catch ME
    fprintf('\n============================================================\n');
    fprintf('程序运行发生错误：%s\n', ME.message);
    for k = 1:length(ME.stack)
        fprintf('  文件 %s，第 %d 行\n', ME.stack(k).name, ME.stack(k).line);
    end
    fprintf('============================================================\n');
    diary off;
    rethrow(ME);
end

diary off;
fprintf('\n运行结束，结果已保存至：\n%s\n', resultFile);
