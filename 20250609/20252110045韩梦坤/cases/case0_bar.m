function result = case0_bar(resultsDir)
%CASE0_BAR 复用2-3算例1：一维两单元杆结构。

    fprintf('\n\n============================================================\n');
    fprintf('算例0A：2-3 一维两单元杆结构接口验证\n');
    fprintf('============================================================\n');

    K=[100,-100,0; -100,300,-200; 0,-200,200];
    F=[0;0;10];
    fixed=1; free=[2,3]; dE=0;
    rhs=F(free)-K(free,fixed)*dE;

    [dF,info]=solve_equilibrium(K(free,free),rhs,'ldlt');
    d=zeros(3,1); d(fixed)=dE; d(free)=dF;
    reaction=K*d-F;

    x=[0;1;2]; E=[100;200]; A=[1;1]; IEN=[1,2;2,3];
    stress=zeros(2,1); axial=zeros(2,1);
    for e=1:2
        n1=IEN(e,1); n2=IEN(e,2); L=x(n2)-x(n1);
        strain=(d(n2)-d(n1))/L;
        stress(e)=E(e)*strain;
        axial(e)=stress(e)*A(e);
    end

    [r,absR,relR]=residual_norm(K(free,free),dF,rhs);
    fprintf('缩减矩阵 K_FF：\n'); disp(K(free,free));
    fprintf('右端项 rhs：\n'); disp(rhs);
    fprintf('LDL^T 最小主元：%.6e\n',info.minPivot);
    fprintf('节点位移 d = [%.8f, %.8f, %.8f]^T\n',d);
    fprintf('理论要求：d2=0.1, d3=0.15\n');
    fprintf('节点1约束反力：%.8f\n',reaction(1));
    fprintf('单元1轴力：%.8f，单元2轴力：%.8f\n',axial(1),axial(2));
    fprintf('残差向量：\n'); disp(r);
    fprintf('相对残差：%.6e\n',relR);

    result.d=d; result.reaction=reaction; result.stress=stress;
    result.axial=axial; result.relativeResidual=relR;
    save(fullfile(resultsDir,'case0_bar.mat'),'result');
end
