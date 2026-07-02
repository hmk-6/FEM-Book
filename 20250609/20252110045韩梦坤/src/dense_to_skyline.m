function sky = dense_to_skyline(K, tolerance)
%DENSE_TO_SKYLINE 将对称矩阵的下三角转换为行轮廓（Skyline）存储。
%每一行 i 仅存储 firstCol(i):i 之间的数值，主对角元位于 diagPtr(i)。

    if nargin < 2
        tolerance = 1.0e-14 * max(1,norm(K,inf));
    end
    n = size(K,1);
    firstCol = zeros(n,1);
    rowLength = zeros(n,1);

    for i = 1:n
        j = find(abs(K(i,1:i)) > tolerance, 1, 'first');
        if isempty(j), j = i; end
        firstCol(i) = j;
        rowLength(i) = i-j+1;
    end

    diagPtr = cumsum(rowLength);
    values = zeros(diagPtr(end),1);
    for i = 1:n
        j0 = firstCol(i);
        idx0 = diagPtr(i) - (i-j0);
        values(idx0:diagPtr(i)) = K(i,j0:i).';
    end

    sky.n = n;
    sky.values = values;
    sky.firstCol = firstCol;
    sky.diagPtr = diagPtr;
end
