function sky = skyline_set(sky, i, j, value)
%SKYLINE_SET 设置 Skyline 下三角存储中的元素。
    if i < j
        tmp=i; i=j; j=tmp;
    end
    if j < sky.firstCol(i)
        if abs(value) > 0
            error('Skyline:OutsideProfile','试图在轮廓外写入非零元素。');
        end
        return;
    end
    idx = sky.diagPtr(i) - (i-j);
    sky.values(idx) = value;
end
