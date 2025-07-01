clc;
clear;
close all;

%% 时序聚类
% 事件序列号筛选
filename = 'XX.xlsx'; % Excel文件的路径
data = readtable(filename, 'PreserveVariableNames', true);  % 保留原始列名
% 确保时间列是 datetime 类型
data.Time = datetime(data.('时间基准点'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS'); % 根据实际时间格式调整
% 12月数据读取和去除离群点
startTime1 = datetime('2022-12-1 00:00:00.000', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');
endTime1 = datetime('2022-12-31 23:59:59.999', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');
[points1,time1, filteredData1] = extractdata(data, startTime1, endTime1);
filteredData1 = sortrows(filteredData1, 'Time'); % 按时间排序
sequenceNumbers1 = (1:height(filteredData1))'; % 生成事件序列号
% 调用聚类分析函数
[isolated1, clusters1, filteredClean] = clusterAnalysis(...
    filteredData1, ...
    points1, ...
    'Threshold', 12, ...
    'ShowAllClusters', true, ...
    'ShowFiltered', true);
% 提取坐标数据  
isolatedCoords1 = points1(isolated1, :);       % 孤立点坐标  
nonIsolatedCoords1 = points1(~isolated1, :);   % 非孤立点坐标 
% 去除离群点
epsilon=50;minPts=10;
[filteredPoints1, numFilteredPoints1] = removeOutliers(nonIsolatedCoords1(:,1), nonIsolatedCoords1(:,2), nonIsolatedCoords1(:,3),  epsilon, minPts);
centerZ = 447;  % 固定Z轴中心坐标
h_values = [100, 120, 140];  % 三组数据对应的高度
L_values = [300, 316, 332];  % 三组数据的长度
threshold = 10;
% 计算x坐标的中位数
x_coords = filteredPoints1(:, 1); 
centerX = median(x_coords);
[collapse3D1, fracture3D1, arcRadiusLeft1, M1, arcAnglemid1, arcRadiusmid1] = process_san1(filteredPoints1, L_values(1), h_values(1), centerX, centerZ);


% (1) 遗传算法参数初始化
pop_size = 100;              % 种群规模
max_generations = 60;        % 进化代数
mutation_prob = 0.2;         % 变异概率
h_range = [5, 50];           % 影响半径范围
elite_ratio = 0.15;          % 精英保留比例
poly_order_fixed = 2;         % 固定多项式阶数(2阶)

% (2) 生成h参数种群
population = h_range(1) + (h_range(2)-h_range(1))*rand(pop_size,1); 
% 初始化进化记录
best_fitness_history = zeros(max_generations,1);

for gen = 1:max_generations
    % (3) 适应度计算
    fitness = zeros(pop_size,1);
    for i = 1:pop_size
        current_h = population(i);
        
        % 计算R²（使用固定参数）
        [~, X_grid, Y_grid, Z_smooth] = mls_fit_function(collapse3D1, current_h,...
                              poly_order_fixed, 50);
        F = scatteredInterpolant(X_grid(:), Y_grid(:), Z_smooth(:), 'natural', 'nearest');
        Z_fit = F(collapse3D1(:,1), collapse3D1(:,2));
        valid = ~isnan(Z_fit);
        SSE = sum((collapse3D1(valid,3) - Z_fit(valid)).^2);
        SST = sum((collapse3D1(valid,3) - mean(collapse3D1(valid,3))).^2);
        fitness(i) = 1 - SSE/SST;
    end
    
    % 记录最佳个体
    [best_fit, best_idx] = max(fitness);
    best_fitness_history(gen) = best_fit;
    best_h = population(best_idx);
    
    % (4) 遗传算子操作
    % 锦标赛选择
    tournament_size = 3;
    selected = zeros(pop_size,1);
    for i = 1:pop_size
        candidates = randperm(pop_size, tournament_size);
        [~, idx] = max(fitness(candidates));
        selected(i) = population(candidates(idx));
    end
    
    % 算术交叉
    offspring = selected;
    for i = 1:2:pop_size-1
        if rand() < 0.7
            alpha = rand();
            parent1 = selected(i);
            parent2 = selected(i+1);
            offspring(i) = alpha*parent1 + (1-alpha)*parent2;
            offspring(i+1) = alpha*parent2 + (1-alpha)*parent1;
            % 边界处理
            offspring(i) = min(max(offspring(i),h_range(1)),h_range(2));
            offspring(i+1) = min(max(offspring(i+1),h_range(1)),h_range(2));
        end
    end
    
    % 高斯变异
    for i = 1:pop_size
        if rand() < mutation_prob
            offspring(i) = offspring(i) + randn()*(h_range(2)-h_range(1))/10;
            offspring(i) = min(max(offspring(i),h_range(1)),h_range(2));
        end
    end
    
    % 精英保留
    [~, elite_idx] = maxk(fitness, round(elite_ratio*pop_size));
    offspring(1:length(elite_idx)) = population(elite_idx);
    
    population = offspring;
    
    % 显示进化信息
    fprintf('Generation %d: R²=%.4f | h=%.1f\n', gen, best_fit, best_h);
    
    % 停止条件
    if gen > 10 && std(best_fitness_history(gen-9:gen)) < 1e-5
        break;
    end
end

% 获取最终优化结果
[~, final_best_idx] = max(fitness);
h_opt = population(final_best_idx);

% 显示优化结果
fprintf('\n优化结果：\n最优h=%.1f\nR²=%.4f\n', h_opt, max(fitness));

% 使用优化参数重新拟合
[~, X_grid, Y_grid, Z_smooth] = mls_fit_function(collapse3D1, h_opt,...
                            poly_order_fixed, 100);

% 步骤4：结果可视化
% ========== 新增：自定义红蓝配色方案 ==========
n = 256; % 色图分辨率
red_component = linspace(1,0,n)';    % 红色分量从1降到0
blue_component = linspace(0,1,n)';   % 蓝色分量从0升到1
green_component = zeros(n,1);        % 绿色分量保持0
redblue = [red_component, green_component, blue_component];
% ==========================================

% 三维曲面图
figure('Name','MLS拟合曲面', 'Position',[100 100 800 600])
surf(X_grid, Y_grid, Z_smooth, 'EdgeColor','none')
hold on
scatter3(collapse3D1(:,1), collapse3D1(:,2), collapse3D1(:,3),...
        20, 'r', 'filled', 'MarkerEdgeColor','k')
title(sprintf('Optimized MLS Fit (h=%.1f, basis=%d)', h_opt));
xlabel('X'); ylabel('Y'); zlabel('Z')
colormap(jet)
colorbar
view(-30,45)
grid on
light('Position',[0 0 100],'Style','local')
lighting gouraud
material dull
legend('拟合曲面', '原始数据')
% 评估指标计算
valid = ~isnan(Z_fit);
SSE = sum((collapse3D1(valid,3) - Z_fit(valid)).^2);
SST = sum((collapse3D1(valid,3) - mean(collapse3D1(valid,3))).^2);
R2 = 1 - SSE/SST;
MAE = mean(abs(residuals(valid)));
RMSE = sqrt(mean(residuals(valid).^2));

fprintf('\n评估指标：\n')
fprintf('决定系数 R²: %.4f\n', R2)
fprintf('平均绝对误差 MAE: %.4f\n', MAE)
fprintf('均方根误差 RMSE: %.4f\n', RMSE)


%融合前期数据对数据进行筛选
%对数据读取和聚类分析，将聚类分析后的数据与前段时间筛选出的数据相结合后再次筛选在继续拟合

xStep=3; yStep =3;
finalPoints2_1 = projectAndDivide1(filteredPoints1, filteredPoints2, xStep, yStep);
% 计算x坐标的中位数
x_coords1 = filteredPoints2(:, 1); 
centerX1 = median(x_coords1);
[collapse3D2, fracture3D2, arcRadiusLeft2, M2, arcAnglemid2, arcRadiusmid2] = process_san1(filteredPoints2, L_values(2), h_values(2), centerX1, centerZ);
