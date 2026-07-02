function K = skyline_to_sparse(sky)
%SKYLINE_TO_SPARSE 将 Skyline 数据恢复为 MATLAB 稀疏对称矩阵。
    n = sky.n;
    nzLower = length(sky.values);
    I = zeros(2*nzLower,1);
    J = zeros(2*nzLower,1);
    V = zeros(2*nzLower,1);
    p = 0;
    for i = 1:n
        for j = sky.firstCol(i):i
            value = skyline_get(sky,i,j);
            p=p+1; I(p)=i; J(p)=j; V(p)=value;
            if i~=j
                p=p+1; I(p)=j; J(p)=i; V(p)=value;
            end
        end
    end
    K=sparse(I(1:p),J(1:p),V(1:p),n,n);
end
