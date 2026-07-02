function result = case0_truss(resultsDir)
%CASE0_TRUSS 复用2-3算例2：二维两杆桁架结构。

    fprintf('\n\n============================================================\n');
    fprintf('算例0B：2-3 二维两杆桁架接口验证\n');
    fprintf('============================================================\n');

    % 与作业给定校核值一致的模型：
    % 节点1(1,0)、节点2(0,0)、节点3(1,1)
    % 单元1:1-3（竖杆），单元2:2-3（斜杆）
    x=[1;0;1]; y=[0;0;1];
    IEN=[1,3;2,3]; E=[1;1]; A=[1;1];
    [K,LM,elementData]=assemble_truss_2d(x,y,IEN,E,A);

    F=zeros(6,1); F(5)=10;
    fixed=[1,2,3,4]; free=[5,6]; dE=zeros(4,1);
    rhs=F(free)-K(free,fixed)*dE;

    [dF,info]=solve_equilibrium(K(free,free),rhs,'ldlt');
    d=zeros(6,1); d(fixed)=dE; d(free)=dF;
    reaction=K*d-F;

    stress=zeros(2,1); axial=zeros(2,1);
    for e=1:2
        n1=IEN(e,1); n2=IEN(e,2);
        dofs=[2*n1-1,2*n1,2*n2-1,2*n2];
        de=d(dofs); L=elementData(e).L; c=elementData(e).c; s=elementData(e).s;
        strain=(-c*de(1)-s*de(2)+c*de(3)+s*de(4))/L;
        stress(e)=E(e)*strain;
        axial(e)=stress(e)*A(e);
    end

    [r,~,relR]=residual_norm(K(free,free),dF,rhs);
    fprintf('对号矩阵 LM：\n'); disp(LM);
    fprintf('总体刚度矩阵 K：\n'); disp(K);
    fprintf('缩减矩阵 K_FF：\n'); disp(K(free,free));
    fprintf('节点3位移：u3=%.6f, v3=%.6f\n',d(5),d(6));
    fprintf('理论校核：u3=38.284271, v3=-10.000000\n');
    fprintf('单元1应力：%.6f，理论约 -10.000000\n',stress(1));
    fprintf('单元2应力：%.6f，理论约 14.142136\n',stress(2));
    fprintf('单元轴力：[%.6f, %.6f]^T\n',axial);
    fprintf('约束反力：\n'); disp(reaction(fixed));
    fprintf('残差向量：\n'); disp(r);
    fprintf('相对残差：%.6e\n',relR);
    fprintf('求解器：%s\n',info.method);

    result.K=K; result.LM=LM; result.d=d; result.reaction=reaction;
    result.stress=stress; result.axial=axial; result.relativeResidual=relR;
    save(fullfile(resultsDir,'case0_truss.mat'),'result');
end
