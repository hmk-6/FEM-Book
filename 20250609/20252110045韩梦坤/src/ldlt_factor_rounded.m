function [L,D,info] = ldlt_factor_rounded(K, digits)
%LDLT_FACTOR_ROUNDED 在每个基本运算后保留 digits 位有效数字。
%用于模拟有限精度/人工四舍五入对病态方程求解的影响。

    n = size(K,1);
    K = round_sig(K, digits);
    L = eye(n);
    D = zeros(n,1);

    for j = 1:n
        s = 0;
        for k = 1:j-1
            term = round_sig(round_sig(L(j,k)^2,digits) * D(k), digits);
            s = round_sig(s + term, digits);
        end
        D(j) = round_sig(K(j,j) - s, digits);
        if D(j) <= 0
            error('LDLT:RoundedNonPositivePivot', ...
                '四舍五入计算后第 %d 个主元非正：%.6e。', j, D(j));
        end

        for i = j+1:n
            s = 0;
            for k = 1:j-1
                term = round_sig(L(i,k)*L(j,k),digits);
                term = round_sig(term*D(k),digits);
                s = round_sig(s + term,digits);
            end
            numerator = round_sig(K(i,j)-s,digits);
            L(i,j) = round_sig(numerator/D(j),digits);
        end
    end

    info.success = true;
    info.digits = digits;
end
