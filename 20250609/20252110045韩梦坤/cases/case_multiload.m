function result = case_multiload(resultsDir)
%CASE_MULTILOAD 比较多载荷工况重复分解与一次分解多次回代。

    fprintf('\n\n============================================================\n');
    fprintf('任务：多载荷工况求解效率比较\n');
    fprintf('============================================================\n');

    rng(2026);
    n=250; nrhs=12;
    e=ones(n,1); K=full(spdiags([-e,2.5*e,-e],[-1,0,1],n,n));
    R=randn(n,nrhs);

    t=tic;
    Xrepeat=zeros(n,nrhs);
    for j=1:nrhs
        [L,D]=ldlt_factor(K);
        Xrepeat(:,j)=ldlt_solve(L,D,R(:,j));
    end
    repeatedTime=toc(t);

    t=tic;
    [L,D]=ldlt_factor(K);
    Xreuse=ldlt_solve(L,D,R);
    reuseTime=toc(t);

    relDifference=norm(Xrepeat-Xreuse,'fro')/norm(Xreuse,'fro');
    speedup=repeatedTime/reuseTime;
    [~,~,relativeResidual]=residual_norm(K,Xreuse,R);

    fprintf('矩阵阶数 n=%d，载荷工况数=%d\n',n,nrhs);
    fprintf('每个右端项重复分解总时间：%.6f s\n',repeatedTime);
    fprintf('一次分解、多右端项求解时间：%.6f s\n',reuseTime);
    fprintf('加速比：%.3f\n',speedup);
    fprintf('两种结果相对差：%.6e\n',relDifference);
    fprintf('多右端项相对残差：%.6e\n',relativeResidual);

    result.n=n; result.nrhs=nrhs; result.repeatedTime=repeatedTime;
    result.reuseTime=reuseTime; result.speedup=speedup;
    result.relativeDifference=relDifference; result.relativeResidual=relativeResidual;
    save(fullfile(resultsDir,'case_multiload.mat'),'result');
end
