function [x,info] = skyline_ldlt_solve(skyL,D,R)
%SKYLINE_LDLT_SOLVE 使用 Skyline 因子求解，支持多个右端项。
    n = skyL.n;
    nrhs = size(R,2);
    t = tic;

    y = R;
    for i = 1:n
        for j = skyL.firstCol(i):i-1
            y(i,:) = y(i,:) - skyline_get(skyL,i,j)*y(j,:);
        end
    end

    z = y ./ D;

    % L^T*x=z。采用逆序行传播，循环只访问轮廓内元素。
    x = z;
    for i = n:-1:1
        for j = skyL.firstCol(i):i-1
            x(j,:) = x(j,:) - skyline_get(skyL,i,j)*x(i,:);
        end
    end

    info.solveTime = toc(t);
    info.nrhs = nrhs;
end
