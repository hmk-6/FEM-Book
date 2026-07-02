function y = round_sig(x, digits)
%ROUND_SIG 将数值四舍五入到指定有效数字位数。
    y = zeros(size(x), 'like', x);
    idx = (x ~= 0) & isfinite(x);
    p = digits - 1 - floor(log10(abs(double(x(idx)))));
    scale = 10.^p;
    y(idx) = cast(round(double(x(idx)).*scale)./scale, 'like', x);
    y(~isfinite(x)) = x(~isfinite(x));
end
