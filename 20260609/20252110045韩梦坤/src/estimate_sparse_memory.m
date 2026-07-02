function bytes = estimate_sparse_memory(A)
%ESTIMATE_SPARSE_MEMORY 估算 MATLAB CSC 稀疏矩阵存储字节数。
%double值8字节，行索引按8字节估计，列指针8字节。
    [~,n]=size(A);
    bytes = 16*nnz(A)+8*(n+1);
end
