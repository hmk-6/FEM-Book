function [K,LM,elementData] = assemble_truss_2d(x,y,IEN,E,A)
%ASSEMBLE_TRUSS_2D 复用 2-3 作业思想组装二维桁架总体刚度矩阵。
    nnp = length(x);
    nel = size(IEN,1);
    ndof = 2;
    K = zeros(nnp*ndof);
    LM = zeros(4,nel);
    elementData = struct([]);

    for e = 1:nel
        n1=IEN(e,1); n2=IEN(e,2);
        dx=x(n2)-x(n1); dy=y(n2)-y(n1);
        L=sqrt(dx^2+dy^2); c=dx/L; s=dy/L;
        ke=E(e)*A(e)/L*[c^2,c*s,-c^2,-c*s; ...
                         c*s,s^2,-c*s,-s^2; ...
                        -c^2,-c*s,c^2,c*s; ...
                        -c*s,-s^2,c*s,s^2];
        dofs=[2*n1-1,2*n1,2*n2-1,2*n2];
        LM(:,e)=dofs.';
        K(dofs,dofs)=K(dofs,dofs)+ke;
        elementData(e).L=L;
        elementData(e).c=c;
        elementData(e).s=s;
        elementData(e).ke=ke;
    end
end
