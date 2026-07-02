function model = assemble_poisson_t3(nx,ny)
%ASSEMBLE_POISSON_T3 在单位正方形上装配 T3 有限元 Poisson 方程。
% -Delta u = f, u=0 on boundary
% 制造解 u=sin(pi*x)sin(pi*y), f=2*pi^2*sin(pi*x)sin(pi*y)

    tTotal = tic;
    hx=1/nx; hy=1/ny;
    [X,Y]=meshgrid(0:hx:1,0:hy:1);
    coordinates=[X(:),Y(:)];
    nnode=size(coordinates,1);

    % 每个矩形划分成两个逆时针 T3 三角形。
    nelem=2*nx*ny;
    elements=zeros(nelem,3);
    e=0;
    for j=1:ny
        for i=1:nx
            n1=(j-1)*(nx+1)+i;
            n2=n1+1;
            n4=j*(nx+1)+i;
            n3=n4+1;
            e=e+1; elements(e,:)=[n1,n2,n3];
            e=e+1; elements(e,:)=[n1,n3,n4];
        end
    end

    % 每个 T3 单元贡献 9 个刚度条目。
    I=zeros(9*nelem,1); J=zeros(9*nelem,1); V=zeros(9*nelem,1);
    F=zeros(nnode,1);
    cursor=0;
    tAssembly=tic;

    % 三点面积坐标积分：每个积分点权重 Area/3。
    bary=[1/6,1/6,2/3; 1/6,2/3,1/6; 2/3,1/6,1/6];

    for ie=1:nelem
        nodes=elements(ie,:);
        xe=coordinates(nodes,1); ye=coordinates(nodes,2);
        detJ=(xe(2)-xe(1))*(ye(3)-ye(1))-(xe(3)-xe(1))*(ye(2)-ye(1));
        area=abs(detJ)/2;
        b=[ye(2)-ye(3); ye(3)-ye(1); ye(1)-ye(2)];
        c=[xe(3)-xe(2); xe(1)-xe(3); xe(2)-xe(1)];
        ke=(b*b'+c*c')/(4*area);

        fe=zeros(3,1);
        for q=1:3
            N=bary(q,:).';
            xq=N.'*xe; yq=N.'*ye;
            fq=2*pi^2*sin(pi*xq)*sin(pi*yq);
            fe=fe+(area/3)*N*fq;
        end

        [ii,jj]=ndgrid(nodes,nodes);
        range=cursor+(1:9);
        I(range)=ii(:); J(range)=jj(:); V(range)=ke(:);
        cursor=cursor+9;
        F(nodes)=F(nodes)+fe;
    end
    K=sparse(I,J,V,nnode,nnode);
    assemblyTime=toc(tAssembly);

    boundary = coordinates(:,1)==0 | coordinates(:,1)==1 | ...
               coordinates(:,2)==0 | coordinates(:,2)==1;
    free=find(~boundary);
    fixed=find(boundary);

    model.nx=nx; model.ny=ny;
    model.coordinates=coordinates;
    model.elements=elements;
    model.K=K; model.F=F;
    model.free=free; model.fixed=fixed;
    model.KFF=K(free,free); model.rhs=F(free);
    model.assemblyTime=assemblyTime;
    model.totalPreparationTime=toc(tTotal);
end
