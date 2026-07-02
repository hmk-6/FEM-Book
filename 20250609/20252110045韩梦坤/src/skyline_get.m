function value = skyline_get(sky, i, j)
%SKYLINE_GET 获取对称 Skyline 矩阵元素。
    if i < j
        tmp=i; i=j; j=tmp;
    end
    if j < sky.firstCol(i)
        value = 0;
    else
        idx = sky.diagPtr(i) - (i-j);
        value = sky.values(idx);
    end
end
