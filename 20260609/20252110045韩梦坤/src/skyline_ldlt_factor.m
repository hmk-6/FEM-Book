function [skyL,D,info] = skyline_ldlt_factor(sky)
%SKYLINE_LDLT_FACTOR 采用活动列/轮廓范围完成 LDL^T 分解。
%输出 skyL 的非对角项存 L，对角项设置为 1；D 单独存储。

    n = sky.n;
    skyL = sky;
    D = zeros(n,1);
    tol = 1.0e-14 * max(1,max(abs(sky.values)));
    t = tic;

    for i = 1:n
        fi = sky.firstCol(i);

        for j = fi:i-1
            startK = max(fi, sky.firstCol(j));
            sumValue = 0;
            for k = startK:j-1
                sumValue = sumValue + ...
                    skyline_get(skyL,i,k) * D(k) * skyline_get(skyL,j,k);
            end
            aij = skyline_get(sky,i,j);
            lij = (aij - sumValue) / D(j);
            skyL = skyline_set(skyL,i,j,lij);
        end

        sumDiag = 0;
        for k = fi:i-1
            lik = skyline_get(skyL,i,k);
            sumDiag = sumDiag + lik*lik*D(k);
        end
        D(i) = skyline_get(sky,i,i) - sumDiag;
        if D(i) <= tol
            error('Skyline:NonPositivePivot', ...
                '矩阵非正定或存在零主元，第%d个主元=%.6e。',i,D(i));
        end
        skyL = skyline_set(skyL,i,i,1.0);
    end

    info.success = true;
    info.factorTime = toc(t);
    info.n = n;
    info.storedEntries = length(sky.values);
    info.minPivot = min(D);
end
