function [x, info] = solve_equilibrium(K_FF, rhs, method, options)
%SOLVE_EQUILIBRIUM 统一平衡方程求解接口
%   method: 'ldlt'、'skyline'、'sparse'、'auto'

    arguments
        K_FF
        rhs
        method char = 'ldlt'
        options.Reorder (1,1) logical = true
        options.PreferExternal (1,1) logical = true
    end

    method = lower(strtrim(method));

    switch method
        case 'ldlt'
            [L,D,fInfo] = ldlt_factor(full(K_FF));
            [x,sInfo] = ldlt_solve(L,D,rhs);
            info = fInfo;
            info.solveTime = sInfo.solveTime;
            info.method = '自行实现稠密 LDL^T';

        case {'skyline','colsol'}
            sky = dense_to_skyline(full(K_FF));
            [skyL,D,fInfo] = skyline_ldlt_factor(sky);
            [x,sInfo] = skyline_ldlt_solve(skyL,D,rhs);
            info = fInfo;
            info.solveTime = sInfo.solveTime;
            info.method = 'Skyline/活动列 LDL^T';

        case {'sparse','mkl_pardiso','pardiso','mumps'}
            [x,info] = solve_sparse_system(sparse(K_FF),rhs, ...
                struct('Reorder',options.Reorder, ...
                       'PreferExternal',options.PreferExternal));

        case 'auto'
            if issparse(K_FF) || size(K_FF,1) > 1000
                [x,info] = solve_sparse_system(sparse(K_FF),rhs, ...
                    struct('Reorder',options.Reorder, ...
                           'PreferExternal',options.PreferExternal));
            else
                [L,D,fInfo] = ldlt_factor(full(K_FF));
                [x,sInfo] = ldlt_solve(L,D,rhs);
                info = fInfo;
                info.solveTime = sInfo.solveTime;
                info.method = '自行实现稠密 LDL^T';
            end

        otherwise
            error('未知求解方法：%s', method);
    end

    [~,~,info.relativeResidual] = residual_norm(K_FF,x,rhs);
end
