function T = case3_ill_conditioned(resultsDir)
%CASE3_ILL_CONDITIONED 病态矩阵的残差、误差和条件数分析。

    fprintf('\n\n============================================================\n');
    fprintf('算例3：病态方程的残差、相对误差与条件数\n');
    fprintf('============================================================\n');

    K=[1.0000,1.0000;1.0000,1.0001];
    aExact=[1;1]; R=K*aExact;
    kappa=cond(K,2);

    [L,D]=ldlt_factor(K);
    aDouble=ldlt_solve(L,D,R);
    [~,~,rrDouble]=residual_norm(K,aDouble,R);
    reDouble=norm(aDouble-aExact)/norm(aExact);

    roundedSuccess=true; roundedMessage=''; aRound=[NaN;NaN]; rrRound=NaN; reRound=NaN;
    try
        Kr=round_sig(K,4); Rr=round_sig(R,4);
        [Lr,Dr]=ldlt_factor_rounded(Kr,4);
        aRound=ldlt_solve_rounded(Lr,Dr,Rr,4);
        rrRound=norm(Rr-Kr*aRound)/max(norm(Rr),eps);
        reRound=norm(aRound-aExact)/norm(aExact);
    catch ME
        roundedSuccess=false; roundedMessage=ME.message;
    end

    Ks=single(K); Rs=single(R);
    aSingle=double(Ks\Rs);
    rrSingle=norm(R-K*aSingle)/norm(R);
    reSingle=norm(aSingle-aExact)/norm(aExact);

    % 额外演示：右端项仅扰动 1e-6，线性方程残差仍可接近机器零，
    % 但解的相对误差被条件数放大到约 1%。
    Rp=R+[0;1e-6];
    ap=ldlt_solve(L,D,Rp);
    rrPert=norm(Rp-K*ap)/norm(Rp);
    rePert=norm(ap-aExact)/norm(aExact);

    fprintf('cond_2(K) = %.6e\n',kappa);
    fprintf('双精度解：[% .10f, % .10f]^T\n',aDouble);
    fprintf('双精度相对残差：%.6e，相对误差：%.6e\n',rrDouble,reDouble);
    if roundedSuccess
        fprintf('4位有效数字解：[% .10f, % .10f]^T\n',aRound);
        fprintf('4位有效数字相对残差：%.6e，相对误差：%.6e\n',rrRound,reRound);
    else
        fprintf('4位有效数字计算失败：%s\n',roundedMessage);
        fprintf('原因：1.0001 四舍五入为 1.000，矩阵变为奇异矩阵。\n');
    end
    fprintf('单精度解：[% .10f, % .10f]^T\n',aSingle);
    fprintf('单精度相对残差：%.6e，相对误差：%.6e\n',rrSingle,reSingle);
    fprintf('右端项扰动1e-6后解：[% .10f, % .10f]^T\n',ap);
    fprintf('扰动问题相对残差：%.6e，相对误差：%.6e\n',rrPert,rePert);
    fprintf(['结论：残差衡量的是计算解对“当前方程”的满足程度；病态矩阵会把输入和舍入误差放大，' ...
             '因此残差很小不能保证相对于真实解的误差也很小。\n']);

    Method=["double";"round4";"single";"rhs_perturbed"];
    RelativeResidual=[rrDouble;rrRound;rrSingle;rrPert];
    RelativeError=[reDouble;reRound;reSingle;rePert];
    Success=[true;roundedSuccess;true;true];
    T=table(Method,Success,RelativeResidual,RelativeError, ...
        'VariableNames',{'Method','Success','RelativeResidual','RelativeError'});
    writetable(T,fullfile(resultsDir,'case3_ill_conditioned.csv'));
end
