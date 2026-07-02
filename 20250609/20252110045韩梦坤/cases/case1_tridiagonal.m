function T = case1_tridiagonal(resultsDir)
%CASE1_TRIDIAGONAL n=10,100,500,1000 的三对角 SPD 方程测试。

    fprintf('\n\n============================================================\n');
    fprintf('算例1：三对角对称正定矩阵规模、时间和内存分析\n');
    fprintf('============================================================\n');

    nList=[10,100,500,1000];
    nCase=length(nList);
    denseTime=zeros(nCase,1); skyTime=zeros(nCase,1); sparseTime=zeros(nCase,1);
    denseRelRes=zeros(nCase,1); skyRelRes=zeros(nCase,1); sparseRelRes=zeros(nCase,1);
    denseRelErr=zeros(nCase,1); skyRelErr=zeros(nCase,1); sparseRelErr=zeros(nCase,1);
    denseMemoryMB=zeros(nCase,1); sparseMemoryMB=zeros(nCase,1); skylineMemoryMB=zeros(nCase,1);
    nnzK=zeros(nCase,1); sparseSolver=strings(nCase,1);

    for ic=1:nCase
        n=nList(ic);
        e=ones(n,1);
        Ksp=spdiags([-e,2*e,-e],[-1,0,1],n,n);
        aExact=ones(n,1); R=Ksp*aExact;
        nnzK(ic)=nnz(Ksp);

        K=full(Ksp);
        t=tic;
        [L,D]=ldlt_factor(K);
        aDense=ldlt_solve(L,D,R);
        denseTime(ic)=toc(t);
        [~,~,denseRelRes(ic)]=residual_norm(K,aDense,R);
        denseRelErr(ic)=norm(aDense-aExact)/norm(aExact);

        t=tic;
        sky=dense_to_skyline(K);
        [skyL,Dsky]=skyline_ldlt_factor(sky);
        aSky=skyline_ldlt_solve(skyL,Dsky,R);
        skyTime(ic)=toc(t);
        [~,~,skyRelRes(ic)]=residual_norm(Ksp,aSky,R);
        skyRelErr(ic)=norm(aSky-aExact)/norm(aExact);

        [aSparse,sInfo]=solve_sparse_system(Ksp,R,struct('Reorder',true,'PreferExternal',true));
        sparseTime(ic)=sInfo.solveTime;
        sparseSolver(ic)=string(sInfo.solverName);
        [~,~,sparseRelRes(ic)]=residual_norm(Ksp,aSparse,R);
        sparseRelErr(ic)=norm(aSparse-aExact)/norm(aExact);

        denseMemoryMB(ic)=8*n*n/1024^2;
        sparseMemoryMB(ic)=estimate_sparse_memory(Ksp)/1024^2;
        skylineMemoryMB(ic)=8*(2*n-1)/1024^2;

        fprintf(['n=%4d: dense %.4fs, skyline %.4fs, sparse %.4fs; ' ...
                 'relerr %.3e / %.3e / %.3e\n'], ...
                 n,denseTime(ic),skyTime(ic),sparseTime(ic), ...
                 denseRelErr(ic),skyRelErr(ic),sparseRelErr(ic));
    end

    T=table(nList.',nnzK,denseTime,skyTime,sparseTime, ...
        denseRelRes,skyRelRes,sparseRelRes, ...
        denseRelErr,skyRelErr,sparseRelErr, ...
        denseMemoryMB,skylineMemoryMB,sparseMemoryMB,sparseSolver, ...
        'VariableNames',{'n','nnzK','DenseLDLT_s','Skyline_s','Sparse_s', ...
        'DenseRelResidual','SkylineRelResidual','SparseRelResidual', ...
        'DenseRelError','SkylineRelError','SparseRelError', ...
        'DenseMemory_MB','SkylineMemory_MB','SparseMemory_MB','SparseSolver'});
    disp(T);
    writetable(T,fullfile(resultsDir,'case1_tridiagonal.csv'));
end
