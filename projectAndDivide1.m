function [finalPoints] = projectAndDivide1(filteredPoints1, filteredPoints2, xStep, yStep)
    % 获取XY平面范围 (使用两组点云的并集范围)
    xMin = min([min(filteredPoints1(:,1)); min(filteredPoints2(:,1))]);
    xMax = max([max(filteredPoints1(:,1)); max(filteredPoints2(:,1))]);
    yMin = min([min(filteredPoints1(:,2)); min(filteredPoints2(:,2))]);
    yMax = max([max(filteredPoints1(:,2)); max(filteredPoints2(:,2))]);
    
    % 创建XY网格划分
    xEdges = xMin:xStep:xMax;
    yEdges = yMin:yStep:yMax;
    
    % 初始化网格单元存储
    gridCells1 = cell(length(yEdges)-1, length(xEdges)-1);
    gridCells2 = cell(length(yEdges)-1, length(xEdges)-1);
    
    % 将点云分配到网格单元
    for i = 1:length(xEdges)-1
        for j = 1:length(yEdges)-1
            % 当前网格边界
            xLow = xEdges(i);
            xHigh = xEdges(i+1);
            yLow = yEdges(j);
            yHigh = yEdges(j+1);
            
            % 筛选在网格内的点
            inCell1 = filteredPoints1(:,1) >= xLow & filteredPoints1(:,1) < xHigh & ...
                      filteredPoints1(:,2) >= yLow & filteredPoints1(:,2) < yHigh;
            gridCells1{j,i} = filteredPoints1(inCell1, :);
            
            inCell2 = filteredPoints2(:,1) >= xLow & filteredPoints2(:,1) < xHigh & ...
                      filteredPoints2(:,2) >= yLow & filteredPoints2(:,2) < yHigh;
            gridCells2{j,i} = filteredPoints2(inCell2, :);
        end
    end

    % 初始化结果点集
    finalPoints = [];
    
    % 遍历所有网格单元进行点云融合
    for i = 1:size(gridCells1, 2)
        for j = 1:size(gridCells1, 1)
            cellPoints1 = gridCells1{j,i};
            cellPoints2 = gridCells2{j,i};
            
            % 处理非空网格单元
            if ~isempty(cellPoints1) || ~isempty(cellPoints2)
                % 双点云网格单元处理
                if ~isempty(cellPoints1) && ~isempty(cellPoints2)
                    % 创建当前网格的点云容器
                    mergedPoints = [];
                    
                    % 处理点云1中的点
                    for k = 1:size(cellPoints1, 1)
                        currentPoint = cellPoints1(k, :);
                        
                        % 在点云2中寻找最近邻点
                        dists = vecnorm(cellPoints2(:,1:2) - currentPoint(1:2), 2, 2);
                        [minDist, idx] = min(dists);
                        
                        if minDist < sqrt(xStep^2 + yStep^2)/2  % 匹配点判断
                            if currentPoint(3) > cellPoints2(idx, 3)
                                mergedPoints = [mergedPoints; currentPoint];
                            else
                                mergedPoints = [mergedPoints; cellPoints2(idx, :)];
                            end
                        else  % 无匹配点
                            mergedPoints = [mergedPoints; currentPoint];
                        end
                    end
                    
                    % 处理点云2中的未匹配点
                    for k = 1:size(cellPoints2, 1)
                        currentPoint = cellPoints2(k, :);
                        
                        % 检查是否已被处理
                        dists = vecnorm(mergedPoints(:,1:2) - currentPoint(1:2), 2, 2);
                        if min(dists) > sqrt(xStep^2 + yStep^2)/2  % 未匹配点
                            mergedPoints = [mergedPoints; currentPoint];
                        end
                    end
                    
                    finalPoints = [finalPoints; mergedPoints];
                    
                % 单点云网格单元处理
                else
                    if ~isempty(cellPoints1)
                        finalPoints = [finalPoints; cellPoints1];
                    end
                    if ~isempty(cellPoints2)
                        finalPoints = [finalPoints; cellPoints2];
                    end
                end
            end
        end
    end
    
    % 空间坐标去重 (基于XY坐标)
    [~, uniqueIdx] = unique(finalPoints(:,1:2), 'rows', 'stable');
    finalPoints = finalPoints(uniqueIdx, :);
end