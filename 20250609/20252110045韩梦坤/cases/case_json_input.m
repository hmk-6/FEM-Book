function result = case_json_input(projectRoot,resultsDir)
%CASE_JSON_INPUT 验证建议JSON输入格式。
    fprintf('\n\n============================================================\n');
    fprintf('附加验证：JSON方程组输入\n');
    fprintf('============================================================\n');
    filename=fullfile(projectRoot,'input','ldlt_test.json');
    data=read_equation_json(filename);
    [x,info]=solve_equilibrium(data.K,data.R,'ldlt');
    [r,absR,relR]=residual_norm(data.K,x,data.R);
    fprintf('标题：%s\n',data.Title);
    fprintf('矩阵阶数：%d\n',size(data.K,1));
    fprintf('解向量：\n'); disp(x);
    fprintf('残差向量：\n'); disp(r);
    fprintf('残差范数：%.6e，相对残差：%.6e\n',absR,relR);
    fprintf('最小主元：%.6e\n',info.minPivot);
    result.x=x; result.relativeResidual=relR; result.info=info;
    save(fullfile(resultsDir,'case_json_input.mat'),'result');
end
