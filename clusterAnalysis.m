% function [isolated, clusterIndices, filteredDataNonIsolated] = clusterAnalysis(filteredData, points, varargin)
% % CLUSTERANALYSIS 三维点集聚类分析及可视化
% % 输入参数：
% %   filteredData1 - 包含时间等信息的表格数据
% %   points        - N×3矩阵，三维坐标点集
% %   varargin      - 可选参数对：
% %       'Threshold'      - 邻接距离阈值（默认13.5）
% %       'ShowAllClusters'- 显示包含孤立点的聚类图（默认true）
% %       'ShowFiltered'   - 显示过滤孤立点后的聚类图（默认true）
% % 输出参数：
% %   isolated              - 孤立点逻辑索引向量
% %   clusterIndices        - 原始聚类标签向量
% %   filteredDataNonIsolated - 过滤孤立点后的数据
% 
% % 解析输入参数
% p = inputParser;
% addParameter(p, 'Threshold', 13.5, @isnumeric);
% addParameter(p, 'ShowAllClusters', true, @islogical);
% addParameter(p, 'ShowFiltered', true, @islogical);
% parse(p, varargin{:});
% threshold = p.Results.Threshold;
% showAll = p.Results.ShowAllClusters;
% showFiltered = p.Results.ShowFiltered;
% 
% % 基础数据准备
% sequenceNumbers = (1:height(filteredData))';
% 
% % 最近邻聚类核心算法
% D = pdist2(points, points); % 欧氏距离矩阵
% A = D <= threshold & ~eye(size(D)); % 邻接矩阵
% G = graph(A);
% clusterIndices = conncomp(G, 'OutputForm', 'vector')';
% isolated = sum(A, 2) == 0; % 孤立点判断
% 
% % % 可视化原始聚类结果
% % if showAll
% %     visualizeClusters(points, clusterIndices, sequenceNumbers, ...
% %         sprintf('原始聚类结果（阈值=%.1f，共%d类）', threshold, max(clusterIndices)), ...
% %         isolated);
% % end
% 
% % 过滤孤立点处理
% filteredDataNonIsolated = filteredData(~isolated, :);
% pointsNonIsolated = points(~isolated, :);
% clusterIndicesFiltered = clusterIndices(~isolated);
% seqNumbersFiltered = sequenceNumbers(~isolated);
% 
% % 可视化过滤结果
% if showFiltered
%     visualizeClusters(pointsNonIsolated, clusterIndicesFiltered, seqNumbersFiltered, ...
%         sprintf('过滤孤立点后结果（共%d类）', max(clusterIndicesFiltered)), ...
%         false(size(pointsNonIsolated,1),1));
% end
% 
% % 输出统计信息
% fprintf('孤立点数量：%d\n', sum(isolated));
% fprintf('有效聚类数量：%d（已排除%d个孤立点）\n',...
%     max(clusterIndicesFiltered), sum(isolated));
% end
% 
% % 可视化子函数
% function visualizeClusters(points, clusters, seqNumbers, titleStr, isIsolated)
%     figure;
%     hold on;
%     axis equal;
%     grid on;
%     view(3);
%     rotate3d on;
%     
%     % 生成颜色映射
%     numClusters = max(clusters);
%     colors = hsv(numClusters);
%     
%     % 绘制每个聚类
%     for k = 1:numClusters
%         mask = (clusters == k);
%         clusterPoints = find(mask);
%         
%         % 绘制点集
%         scatter3(points(clusterPoints,1),...
%                  points(clusterPoints,2),...
%                  points(clusterPoints,3),...
%                  40, colors(k,:), 'filled');
%         
%         % 添加序列号标注
%         for i = 1:length(clusterPoints)
%             idx = clusterPoints(i);
%             text(points(idx,1), points(idx,2), points(idx,3)+1,...
%                 num2str(seqNumbers(idx)),...
%                 'FontSize',7, 'Color','k',...
%                 'HorizontalAlignment','center');
%         end
%     end
%     
%     % 绘制孤立点（如果存在）
%     if any(isIsolated)
%         scatter3(points(isIsolated,1),...
%                 points(isIsolated,2),...
%                 points(isIsolated,3),...
%                 40, [0.5 0.5 0.5], 'filled', 'MarkerEdgeColor','k');
%     end
%     
%     % 添加标注
%     xlabel('X坐标');
%     ylabel('Y坐标');
%     zlabel('Z坐标');
%     title(titleStr);
%     hold off;
% end


function [isolated, clusterIndices, filteredDataNonIsolated] = clusterAnalysis(filteredData, points, varargin)
% CLUSTERANALYSIS 三维点集聚类分析及可视化
% 输入参数：
%   filteredData - 包含时间等信息的表格数据
%   points       - N×3矩阵，三维坐标点集
%   varargin     - 可选参数对：
%       'Threshold'      - 邻接距离阈值（默认13.5）
%       'ShowAllClusters'- 显示包含孤立点的聚类图（默认true）
%       'ShowFiltered'   - 显示过滤孤立点后的聚类图（默认true）
% 输出参数：
%   isolated              - 孤立点逻辑索引向量
%   clusterIndices        - 原始聚类标签向量
%   filteredDataNonIsolated - 过滤孤立点后的数据

% 解析输入参数
p = inputParser;
addParameter(p, 'Threshold', 13.5, @isnumeric);
addParameter(p, 'ShowAllClusters', true, @islogical);
addParameter(p, 'ShowFiltered', true, @islogical);
parse(p, varargin{:});
threshold = p.Results.Threshold;
showAll = p.Results.ShowAllClusters;
showFiltered = p.Results.ShowFiltered;

% 基础数据准备
sequenceNumbers = (1:height(filteredData))';

% 最近邻聚类核心算法
D = pdist2(points, points); % 欧氏距离矩阵
A = D <= threshold & ~eye(size(D)); % 邻接矩阵
G = graph(A);
clusterIndices = conncomp(G, 'OutputForm', 'vector')';

% 计算簇大小并确定孤立点（簇大小<5的）
clusterSizes = accumarray(clusterIndices, 1); 
isolated = clusterSizes(clusterIndices) <3; % 修改后的孤立点判断条件

% % 可视化原始聚类结果
% if showAll
%     visualizeClusters(points, clusterIndices, sequenceNumbers, ...
%         sprintf('原始聚类结果（阈值=%.1f，共%d类）', threshold, max(clusterIndices)), ...
%         isolated);
% end

% 过滤孤立点处理
filteredDataNonIsolated = filteredData(~isolated, :);
pointsNonIsolated = points(~isolated, :);
clusterIndicesFiltered = clusterIndices(~isolated);
seqNumbersFiltered = sequenceNumbers(~isolated);

% % 可视化过滤结果
% if showFiltered
%     visualizeClusters(pointsNonIsolated, clusterIndicesFiltered, seqNumbersFiltered, ...
%         sprintf('过滤孤立点后结果（共%d类）', numel(unique(clusterIndicesFiltered))), ... % 修正统计显示
%         false(size(pointsNonIsolated,1),1));
% end

% 输出统计信息
uniqueClustersFiltered = unique(clusterIndicesFiltered);
numClustersFiltered = numel(uniqueClustersFiltered); % 正确统计有效聚类数量
fprintf('孤立点数量：%d\n', sum(isolated));
fprintf('有效聚类数量：%d（已排除%d个孤立点）\n',...
    numClustersFiltered, sum(isolated));
end

% % 可视化子函数（保持不变）
% function visualizeClusters(points, clusters, seqNumbers, titleStr, isIsolated)
%     figure;
%     hold on;
%     axis equal;
%     grid on;
%     view(3);
%     rotate3d on;
%     
%     % 生成颜色映射
%     numClusters = max(clusters);
%     colors = hsv(numClusters);
%     
%     % 绘制每个聚类
%     for k = 1:numClusters
%         mask = (clusters == k);
%         clusterPoints = find(mask);
%         
%         % 绘制点集
%         scatter3(points(clusterPoints,1),...
%                  points(clusterPoints,2),...
%                  points(clusterPoints,3),...
%                  40, colors(k,:), 'filled');
%         
%         % 添加序列号标注
%         for i = 1:length(clusterPoints)
%             idx = clusterPoints(i);
%             text(points(idx,1), points(idx,2), points(idx,3)+1,...
%                 num2str(seqNumbers(idx)),...
%                 'FontSize',7, 'Color','k',...
%                 'HorizontalAlignment','center');
%         end
%     end
%     
%     % 绘制孤立点（如果存在）
%     if any(isIsolated)
%         scatter3(points(isIsolated,1),...
%                 points(isIsolated,2),...
%                 points(isIsolated,3),...
%                 40, [0.5 0.5 0.5], 'filled', 'MarkerEdgeColor','k');
%     end
%     
%     % 添加标注
%     xlabel('X坐标');
%     ylabel('Y坐标');
%     zlabel('Z坐标');
%     title(titleStr);
%     hold off;
% end