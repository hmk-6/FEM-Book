function [a, info] = ldlt_solve(L, D, R)
%LDLT_SOLVE 求解 L*D*L^T*a = R，支持多个右端项。

    n = size(L,1);
    if size(L,2) ~= n || length(D) ~= n || size(R,1) ~= n
        error('LDLT:DimensionMismatch', 'L、D 与 R 的维数不一致。');
    end

    D = D(:);
    nrhs = size(R,2);
    tSolve = tic;

    % 前代：L*y = R，L 对角线为 1。
    y = zeros(n, nrhs);
    for i = 1:n
        if i == 1
            y(i,:) = R(i,:);
        else
            y(i,:) = R(i,:) - L(i,1:i-1) * y(1:i-1,:);
        end
    end

    % 对角求解：D*z = y。
    z = y ./ D;

    % 回代：L^T*a = z。
    a = zeros(n, nrhs);
    for i = n:-1:1
        if i == n
            a(i,:) = z(i,:);
        else
            a(i,:) = z(i,:) - L(i+1:n,i)' * a(i+1:n,:);
        end
    end

    info.solveTime = toc(tSolve);
    info.nrhs = nrhs;
end
