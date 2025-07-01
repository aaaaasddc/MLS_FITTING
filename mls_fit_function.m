% 步骤1：定义MLS拟合函数  高斯权函数、多项式基函数
function [RMSE, X_grid, Y_grid, Z_smooth] = mls_fit_function(collapse3D1, h, poly_order, num_grid_points)
    x = collapse3D1(:,1);
    y = collapse3D1(:,2);
    z = collapse3D1(:,3);
    
    % 生成网格
    x_min = min(x); x_max = max(x);
    y_min = min(y); y_max = max(y);
    [X_grid, Y_grid] = meshgrid(linspace(x_min, x_max, num_grid_points),...
                        linspace(y_min, y_max, num_grid_points));
    grid_points = [X_grid(:), Y_grid(:)];
    
    % 确定基函数
    if poly_order == 2
        min_pts = 6; % 原为4
        basis_func = @(x,y) [ones(size(x)), x, y, x.^2,y.^2,x.*y];

    else
        min_pts = 3;
        basis_func = @(x,y) [ones(size(x)), x, y];
    end
    
    % 构建k-d树
    data_points = [x, y];
    [idx, dist] = rangesearch(data_points, grid_points, h);
    
    % MLS拟合
    Z_grid = nan(size(X_grid));
    for k = 1:size(grid_points, 1)
        if numel(idx{k}) < min_pts, continue; end
        x_neigh = x(idx{k});
        y_neigh = y(idx{k});
        z_neigh = z(idx{k});
        d = dist{k}(:);
        
        % 计算权重
        w = exp(-(d.^2)/(h^2));
        W = diag(w);
        
        % 构建矩阵
        dx = x_neigh - grid_points(k,1);
        dy = y_neigh - grid_points(k,2);
        A = basis_func(dx, dy);
        
        % 加权最小二乘
%         A_weighted = W * A;
%         b_weighted = W * z_neigh;
%         beta = A_weighted \ b_weighted;

    % 加入正则化（在加权最小二乘部分）
    lambda = 1e-6; % 正则化系数
    A_weighted = W * A;
    b_weighted = W * z_neigh;
    beta = (A_weighted' * A_weighted + lambda * eye(size(A,2))) \ (A_weighted' * b_weighted);
        % 预测中心点
        Z_grid(k) = basis_func(0,0) * beta;
    end
    
    % 后处理
    Z_filled = fillmissing(Z_grid, 'movmedian', 5);  % 填充缺失值
    Z_smooth = imgaussfilt(Z_filled, 1);             % 高斯平滑
    
    % 插值评估
    valid_grid = ~isnan(Z_smooth);
    F = scatteredInterpolant(X_grid(valid_grid), Y_grid(valid_grid),...
                            Z_smooth(valid_grid), 'natural', 'nearest');
    Z_interp = F(x, y);
    
    % 计算RMSE
    valid = ~isnan(Z_interp);
    RMSE = sqrt(mean((z(valid) - Z_interp(valid)).^2));
end
