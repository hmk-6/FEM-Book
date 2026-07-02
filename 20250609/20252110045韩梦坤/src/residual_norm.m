function [r, absResidual, relativeResidual] = residual_norm(K, a, R)
%RESIDUAL_NORM 计算残差向量、绝对残差范数与相对残差。
    r = R - K*a;
    absResidual = norm(r, 2);
    relativeResidual = absResidual / max(norm(R,2), eps);
end
