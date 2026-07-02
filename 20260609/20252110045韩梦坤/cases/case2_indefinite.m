function result = case2_indefinite(resultsDir)
%CASE2_INDEFINITE 检测非正定矩阵的非正主元。

    fprintf('\n\n============================================================\n');
    fprintf('算例2：非正定矩阵检测\n');
    fprintf('============================================================\n');
    K=[1,2;2,1]; R=[1;1];
    detected=false; message='';
    try
        [L,D]=ldlt_factor(K); %#ok<NASGU>
        fprintf('警告：程序未检测到非正定性。\n');
    catch ME
        detected=true; message=ME.message;
        fprintf('成功检测到非正主元，LDL^T 求解停止。\n');
        fprintf('错误提示：%s\n',ME.message);
    end
    fprintf(['说明：有限元模型缺少足够位移边界条件时会保留刚体模态，' ...
             '缩减刚度矩阵可能出现零主元；材料或单元设置不合理也可能导致非正主元。\n']);
    result.detected=detected; result.message=message; result.K=K; result.R=R;
    save(fullfile(resultsDir,'case2_indefinite.mat'),'result');
end
