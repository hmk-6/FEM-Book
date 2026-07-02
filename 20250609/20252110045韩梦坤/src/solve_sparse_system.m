function [x,info] = solve_sparse_system(K,b,options)
%SOLVE_SPARSE_SYSTEM 稀疏对称正定方程统一接口。
%优先尝试 PARDISO Project MATLAB 接口；若未安装，则使用 MATLAB
%稀疏反斜杠作为可运行的后备方案，并在结果中明确标注。

    if nargin < 3, options = struct(); end
    if ~isfield(options,'Reorder'), options.Reorder=true; end
    if ~isfield(options,'PreferExternal'), options.PreferExternal=true; end

    K = sparse(K);
    n = size(K,1);
    info.n = n;
    info.nnz = nnz(K);
    info.format = 'MATLAB CSC';
    info.externalAvailable = false;
    info.reordered = false;
    info.permutation = (1:n).';

    % 常见 PARDISO Project MATLAB 接口自动检测。
    hasPardisoProject = all(cellfun(@(f) exist(f,'file')~=0, ...
        {'pardisoinit','pardisoreorder','pardisofactor','pardisosolve','pardisofree'}));

    if options.PreferExternal && hasPardisoProject
        try
            verbose = false;
            % 2 表示实对称正定矩阵；部分版本使用 11/2 等参数约定。
            infoP = pardisoinit(2,0);
            Atri = tril(K);
            t = tic;
            infoP = pardisoreorder(Atri,infoP,verbose);
            infoP = pardisofactor(Atri,infoP,verbose);
            [x,infoP] = pardisosolve(Atri,b,infoP,verbose);
            solveTime = toc(t);
            pardisofree(infoP);
            info.solverName = 'PARDISO Project MATLAB interface';
            info.solveTime = solveTime;
            info.externalAvailable = true;
            [~,~,info.relativeResidual] = residual_norm(K,x,b);
            return;
        catch ME
            warning('PARDISO接口调用失败，改用MATLAB稀疏求解：%s',ME.message);
            info.externalFailure = ME.message;
        end
    end

    % MATLAB 稀疏矩阵采用 CSC 存储。对 SPD 矩阵，反斜杠会调用
    % MATLAB 可用的高性能稀疏直接求解内核。
    if options.Reorder
        p = symamd(K);
        info.reordered = true;
        info.permutation = p;
        Kp = K(p,p);
        bp = b(p,:);
        t = tic;
        xp = Kp \ bp;
        info.solveTime = toc(t);
        x = zeros(size(b));
        x(p,:) = xp;
        info.solverName = 'MATLAB sparse backslash + SYMAMD (fallback)';
    else
        t = tic;
        x = K \ b;
        info.solveTime = toc(t);
        info.solverName = 'MATLAB sparse backslash (fallback)';
    end

    [~,~,info.relativeResidual] = residual_norm(K,x,b);
end
