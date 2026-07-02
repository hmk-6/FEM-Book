function result = case_skyline(resultsDir)
%CASE_SKYLINE 活动列/轮廓存储方法验证。

    fprintf('\n\n============================================================\n');
    fprintf('任务：活动列/轮廓（Skyline）存储与求解\n');
    fprintf('============================================================\n');

    n=1000; e=ones(n,1);
    Ksp=spdiags([-e,2*e,-e],[-1,0,1],n,n);
    K=full(Ksp); exact=ones(n,1); R=Ksp*exact;
    sky=dense_to_skyline(K);
    [skyL,D,fi]=skyline_ldlt_factor(sky);
    [x,si]=skyline_ldlt_solve(skyL,D,R);
    [~,~,rr]=residual_norm(Ksp,x,R);
    relErr=norm(x-exact)/norm(exact);

    denseEntries=n*n; skyEntries=length(sky.values); sparseEntries=nnz(Ksp);
    denseMB=8*denseEntries/1024^2;
    skyMB=8*skyEntries/1024^2;
    sparseMB=estimate_sparse_memory(Ksp)/1024^2;

    fprintf('n=%d\n',n);
    fprintf('稠密存储元素数：%d，约 %.3f MB\n',denseEntries,denseMB);
    fprintf('Skyline存储元素数：%d，约 %.6f MB\n',skyEntries,skyMB);
    fprintf('稀疏矩阵非零元：%d，CSC估算 %.6f MB\n',sparseEntries,sparseMB);
    fprintf('Skyline分解时间：%.6f s，求解时间：%.6f s\n',fi.factorTime,si.solveTime);
    fprintf('相对残差：%.6e，相对误差：%.6e\n',rr,relErr);
    fprintf(['说明：Skyline仅在每一行的首个非零项至对角线之间循环，' ...
             '避免对轮廓外零元素进行存储和计算。\n']);

    result.n=n; result.denseEntries=denseEntries; result.skyEntries=skyEntries;
    result.sparseEntries=sparseEntries; result.relativeResidual=rr;
    result.relativeError=relErr; result.factorTime=fi.factorTime; result.solveTime=si.solveTime;
    save(fullfile(resultsDir,'case_skyline.mat'),'result');
end
