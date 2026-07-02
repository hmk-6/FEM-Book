function [L, D, info] = ldlt_factor(K, options)
%LDLT_FACTOR 自行实现的稠密对称正定矩阵 LDL^T 分解
%   K = L*diag(D)*L'
%   L 为单位下三角矩阵，D 为对角元素向量。
%
%   不调用 ldl、chol、mldivide 或线性方程组求解库完成核心分解。

    arguments
        K double
        options.SymmetryTolerance (1,1) double = 1.0e-12
        options.PivotTolerance (1,1) double = 1.0e-14
    end

    [n, m] = size(K);
    if n ~= m
        error('LDLT:NotSquare', '输入矩阵必须为方阵。');
    end
    if any(~isfinite(K), 'all')
        error('LDLT:InvalidValue', '输入矩阵包含 NaN 或 Inf。');
    end

    scale = max(1.0, norm(K, inf));
    symmetryError = norm(K - K', inf) / scale;
    if symmetryError > options.SymmetryTolerance
        error('LDLT:NotSymmetric', ...
            '矩阵不满足对称性要求，相对不对称量为 %.3e。', symmetryError);
    end

    % 消除仅由浮点误差造成的微小不对称。
    A = 0.5 * (K + K');
    L = eye(n);
    D = zeros(n, 1);
    pivotTol = options.PivotTolerance * scale;

    tFactor = tic;
    for j = 1:n
        % A(j,j) 已经是消去前 j-1 列后的 Schur 补主元。
        pivot = A(j,j);
        if pivot <= pivotTol
            error('LDLT:NonPositivePivot', ...
                ['矩阵非正定或存在零主元：第 %d 个主元 = %.6e，' ...
                 '允许下限 = %.6e。'], j, pivot, pivotTol);
        end

        D(j) = pivot;

        if j < n
            L(j+1:n, j) = A(j+1:n, j) / D(j);

            % 对剩余 Schur 补作秩一更新。
            v = L(j+1:n, j);
            A(j+1:n, j+1:n) = A(j+1:n, j+1:n) - D(j) * (v * v');

            % 保持数值对称，抑制累积舍入误差。
            A(j+1:n, j+1:n) = 0.5 * ...
                (A(j+1:n, j+1:n) + A(j+1:n, j+1:n)');
        end
    end

    info.success = true;
    info.n = n;
    info.factorTime = toc(tFactor);
    info.minPivot = min(D);
    info.maxPivot = max(D);
    info.symmetryError = symmetryError;
    info.reconstructionError = norm(K - L*diag(D)*L', 'fro') / max(norm(K,'fro'), eps);
end
