function a = ldlt_solve_rounded(L,D,R,digits)
%LDLT_SOLVE_ROUNDED 采用指定有效数字完成前代、对角求解和回代。
    n = size(L,1);
    R = round_sig(R,digits);
    y = zeros(n,1);
    for i = 1:n
        s = 0;
        for j = 1:i-1
            s = round_sig(s + round_sig(L(i,j)*y(j),digits),digits);
        end
        y(i) = round_sig(R(i)-s,digits);
    end

    z = zeros(n,1);
    for i = 1:n
        z(i) = round_sig(y(i)/D(i),digits);
    end

    a = zeros(n,1);
    for i = n:-1:1
        s = 0;
        for j = i+1:n
            s = round_sig(s + round_sig(L(j,i)*a(j),digits),digits);
        end
        a(i) = round_sig(z(i)-s,digits);
    end
end
