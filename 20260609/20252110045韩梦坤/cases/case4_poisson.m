function summary = case4_poisson(resultsDir,figuresDir,nList)
%CASE4_POISSON T3有限元求解二维Poisson方程。

    fprintf('\n\n============================================================\n');
    fprintf('算例4：二维 Poisson 方程 T3 有限元稀疏求解\n');
    fprintf('============================================================\n');


    % 小规模网格与自行实现稠密 LDL^T 对比（满足作业中的小规模校核要求）。
    smallModel=assemble_poisson_t3(10,10);
    [Ls,Ds]=ldlt_factor(full(smallModel.KFF));
    usDense=ldlt_solve(Ls,Ds,smallModel.rhs);
    [usSparse,smallInfo]=solve_sparse_system(smallModel.KFF,smallModel.rhs, ...
        struct('Reorder',true,'PreferExternal',true));
    smallDifference=norm(usDense-usSparse)/max(norm(usDense),eps);
    fprintf('小规模 nx=ny=10：稠密LDL^T与稀疏求解相对差=%.6e，稀疏求解器=%s\n', ...
        smallDifference,smallInfo.solverName);

    nc=length(nList);
    nxCol=zeros(nc,1); nodes=zeros(nc,1); elements=zeros(nc,1); unknowns=zeros(nc,1);
    nnzK=zeros(nc,1); assemblyTime=zeros(nc,1); bcTime=zeros(nc,1);
    solveTime=zeros(nc,1); totalTime=zeros(nc,1); relResidual=zeros(nc,1);
    maxError=zeros(nc,1); l2RelativeError=zeros(nc,1); memoryMB=zeros(nc,1);
    solverName=strings(nc,1);

    for ic=1:nc
        nx=nList(ic); ny=nx;
        tAll=tic;
        model=assemble_poisson_t3(nx,ny);
        tBC=tic;
        KFF=model.KFF; rhs=model.rhs;
        bcTime(ic)=toc(tBC);

        [uF,sInfo]=solve_sparse_system(KFF,rhs,struct('Reorder',true,'PreferExternal',true));
        u=zeros(size(model.coordinates,1),1); u(model.free)=uF;
        uExact=sin(pi*model.coordinates(:,1)).*sin(pi*model.coordinates(:,2));
        err=abs(u-uExact);

        [~,~,relResidual(ic)]=residual_norm(KFF,uF,rhs);
        maxError(ic)=max(err);
        l2RelativeError(ic)=norm(u-uExact)/norm(uExact);

        nxCol(ic)=nx;
        nodes(ic)=size(model.coordinates,1);
        elements(ic)=size(model.elements,1);
        unknowns(ic)=length(model.free);
        nnzK(ic)=nnz(KFF);
        assemblyTime(ic)=model.assemblyTime;
        solveTime(ic)=sInfo.solveTime;
        totalTime(ic)=toc(tAll);
        memoryMB(ic)=estimate_sparse_memory(KFF)/1024^2;
        solverName(ic)=string(sInfo.solverName);

        fprintf(['nx=ny=%d: nodes=%d, elements=%d, unknowns=%d, nnz=%d, ' ...
                 'assembly=%.4fs, solve=%.4fs, residual=%.3e, maxErr=%.3e, L2rel=%.3e\n'], ...
                 nx,nodes(ic),elements(ic),unknowns(ic),nnzK(ic), ...
                 assemblyTime(ic),solveTime(ic),relResidual(ic), ...
                 maxError(ic),l2RelativeError(ic));
        fprintf('求解器：%s\n',sInfo.solverName);

        % 最小网格生成数值解与误差图；所有网格生成中心线误差数据。
        if ic==1
            tri=model.elements;
            fig=figure('Visible','off');
            trisurf(tri,model.coordinates(:,1),model.coordinates(:,2),u,'EdgeColor','none');
            xlabel('x'); ylabel('y'); zlabel('u_h'); title(sprintf('T3 FEM数值解 nx=ny=%d',nx));
            view(45,30); colorbar; grid on;
            exportgraphics(fig,fullfile(figuresDir,'poisson_solution_50.png'),'Resolution',200);
            close(fig);

            fig=figure('Visible','off');
            trisurf(tri,model.coordinates(:,1),model.coordinates(:,2),err,'EdgeColor','none');
            xlabel('x'); ylabel('y'); zlabel('|u_h-u_{exact}|'); title(sprintf('节点绝对误差 nx=ny=%d',nx));
            view(45,30); colorbar; grid on;
            exportgraphics(fig,fullfile(figuresDir,'poisson_error_50.png'),'Resolution',200);
            close(fig);
        end
    end

    summary=table(nxCol,nodes,elements,unknowns,nnzK,assemblyTime,bcTime,solveTime,totalTime, ...
        memoryMB,relResidual,maxError,l2RelativeError,solverName, ...
        'VariableNames',{'nx','Nodes','Elements','UnknownDOF','nnzK','Assembly_s','BC_s', ...
        'Solve_s','Total_s','SparseMemory_MB','RelativeResidual','MaxNodalError', ...
        'DiscreteL2RelativeError','Solver'});
    disp(summary);
    writetable(summary,fullfile(resultsDir,'case4_poisson_summary.csv'));

    fig=figure('Visible','off');
    loglog(nList,maxError,'o-','LineWidth',1.3); hold on;
    loglog(nList,l2RelativeError,'s-','LineWidth',1.3);
    xlabel('每方向单元数 n_x=n_y'); ylabel('误差');
    legend('最大节点误差','离散L2相对误差','Location','southwest');
    title('Poisson有限元误差收敛'); grid on;
    exportgraphics(fig,fullfile(figuresDir,'poisson_convergence.png'),'Resolution',200);
    close(fig);

    fprintf(['稀疏格式说明：MATLAB sparse采用CSC格式；每列记录非零值、行索引和列指针。' ...
             'SYMAMD重排序用于减少Cholesky/LDL^T分解中的填充元。\n']);
end
