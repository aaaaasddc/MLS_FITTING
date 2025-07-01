%%  原版本
% function [collapse3D, fracture3D, arcRadiusLeft, M1, arcAnglemid, arcRadiusmid] = process_san1(currentPoints, L, h, D_x, D_z)
%     % 处理单组数据的崩落点筛选（三心圆拱约束）
%     % 输入输出参数同原函数
%     
%     % 矩形 AFEG 的顶点坐标（二维XZ平面）
%     A = [D_x - L/2, D_z];      % 左下角
%     F = [D_x + L/2, D_z];      % 右下角
%     E = [D_x + L/2, D_z + h];  % 右上角
%     G = [D_x - L/2, D_z + h];  % 左上角
%     
%     % 计算角平分线交点 M1
%     C = [D_x, D_z + h];  % 上顶点
%     D = [D_x, D_z];      % AF中点
%     
%     % 向量计算
%     vec_CG = G - C;
%     vec_CA = A - C;
%     vec_AG = G - A;
%     vec_AC = C - A;
%     
%     % 单位向量
%     unit_CG = vec_CG / norm(vec_CG);
%     unit_CA = vec_CA / norm(vec_CA);
%     unit_AG = vec_AG / norm(vec_AG);
%     unit_AC = vec_AC / norm(vec_AC);
%     
%     % 角平分线方向
%     dir_C_bisector = unit_CG + unit_CA;  % C点的角平分线方向向量
%     dir_A_bisector = unit_AG + unit_AC;  % A点的角平分线方向向量
%     
%     % 解联立方程求交点 M1
%     % 参数方程: C + t*dir_C_bisector = A + s*dir_A_bisector
%     % 转换为线性方程组求解t和s
%     A_matrix = [dir_C_bisector(1), -dir_A_bisector(1); 
%                 dir_C_bisector(2), -dir_A_bisector(2)];
%     b_vector = [A(1) - C(1); 
%                 A(2) - C(2)];
%     solution = A_matrix \ b_vector;
%     t = solution(1);
%     M1 = C + t * dir_C_bisector;  % 交点坐标
%     
%     % 计算中间圆弧圆心O（过M1作AC的垂线与CD的交点）
%     AC_slope = (C(2) - A(2)) / (C(1) - A(1));
%     if AC_slope == 0
%         % AC水平，垂线为竖直线x=M1(1)
%         O = [D_x, M1(2)];  % CD线为x=D_x
%     else
%         perp_slope = -1 / AC_slope;
%         % 垂线方程: z = perp_slope*(x - M1(1)) + M1(2)
%         % CD线方程: x = D_x
%         O_x = D_x;
%         O_z = perp_slope * (D_x - M1(1)) + M1(2);
%         O = [O_x, O_z];
%     end
%     
%     % 计算左侧圆弧圆心O1（MO与AF的交点）
%     if O(1) == M1(1)
%         O1 = [M1(1), A(2)];  % 垂直线
%     else
%         MO_slope = (O(2) - M1(2)) / (O(1) - M1(1));
%         O1_x = M1(1) + (A(2) - M1(2)) / MO_slope;
%         O1 = [O1_x, A(2)];
%     end
% 
%     % 绘制圆弧 MCX (以 O 为圆心)
%     theta_M = atan2(M1(2)-O(2), M1(1)-O(1));
%     theta_C = atan2(C(2)-O(2), C(1)-O(1));
%     theta = linspace(theta_M, theta_C, 100);
%     arc_MCX_x = O(1) + norm(O-M1)*cos(theta);
%     arc_MCX_z = O(2) + norm(O-M1)*sin(theta);
%     
%     % 绘制圆弧 AM (以 O1 为圆心)
%     theta_A = atan2(A(2)-O1(2), A(1)-O1(1));
%     theta_M_O1 = atan2(M1(2)-O1(2), M1(1)-O1(1));
%     theta_AM = linspace(theta_A, theta_M_O1, 100);
%     arc_AM_x = O1(1) + norm(O1-A)*cos(theta_AM);
%     arc_AM_z = O1(2) + norm(O1-A)*sin(theta_AM);
%     
%     % 镜像右侧
%     arc_MKF_x = 2*D_x - arc_MCX_x;  % 关于 x = D_x 镜像
%     arc_AM_mirror_x = 2*D_x - arc_AM_x;
%     
%     % 右侧对称点
%     O2 = [2*D_x - O1(1), O1(2)];
%     M2 = [2*D_x - M1(1), M1(2)];
%     
%     % 计算圆弧参数
%     arcRadiusLeft = norm(O1 - A);   % 左侧圆弧半径
%     arcRadiusmid = norm(O - C);     % 中间圆弧半径
%     disp(arcRadiusmid);
%     disp(arcRadiusLeft);
%     % 角度范围处理函数
%     adjust_angle = @(theta) mod(theta + 2*pi, 2*pi);
%     
%     % 中间圆弧角度范围
%     theta_M = adjust_angle(atan2(M1(2)-O(2), M1(1)-O(1)));
%     theta_C = adjust_angle(atan2(C(2)-O(2), C(1)-O(1)));
%     theta_M2 = adjust_angle(atan2(M2(2)-O(2), M2(1)-O(1)));
%     arcAnglemid = adjust_angle(theta_C - theta_M);
% 
%     % 初始化输出
%     collapse3D = [];
%     fracture3D = [];
%     
%     for i = 1:size(currentPoints, 1)
%         pt = currentPoints(i, [1,3]); % 提取XZ坐标
%         x = pt(1); z = pt(2);
%         
%         % 1. 排除矩形外点
%         if x < A(1) || x > F(1) || z < D_z || z > D_z + h
%             fracture3D = [fracture3D; currentPoints(i, :)];
%             continue;
%         end
%         
%         % 2. 判断点是否在任一圆弧内
%         inArc = false;
%         
%         % 左侧圆弧判断
%         if check_arc(pt, O1, A, M1, arcRadiusLeft)
%             inArc = true;
%         end      
%         % 中间圆弧判断
%         if check_arc_mid(pt, O, theta_M, theta_C, arcRadiusmid)
%             inArc = true;
%         end
%         if check_arc_mid(pt, O, theta_C, theta_M2, arcRadiusmid)
%             inArc = true;
%         end
%         % 右侧圆弧判断
%         if check_arc(pt, O2, F, M2, arcRadiusLeft)
%             inArc = true;
%         end
%         
%         if inArc
%             collapse3D = [collapse3D; currentPoints(i, :)];
%         else
%             fracture3D = [fracture3D; currentPoints(i, :)];
%         end
%     end
%     
%     % 去除重复点
%     [collapse3D, ~] = unique(collapse3D, 'rows', 'stable');
%     [fracture3D, ~] = unique(fracture3D, 'rows', 'stable');
%     
%     % 在绘制三维图形之前，绘制二维三心圆拱
%     figure;
%     hold on;
%     plot(arc_MCX_x, arc_MCX_z, 'r', 'DisplayName', 'MCX');
%     plot(arc_AM_x, arc_AM_z, 'g', 'DisplayName', 'AM');
%     plot(arc_MKF_x, arc_MCX_z, 'b', 'DisplayName', 'MKF');
%     plot(arc_AM_mirror_x, arc_AM_z, 'm', 'DisplayName', 'AM_mirror');
%     plot([A(1), F(1), F(1), A(1), A(1)], [A(2), A(2), E(2), E(2), A(2)], 'k--', 'DisplayName', 'AFEG');
%     % 绘制二维崩落点
%     scatter(collapse3D(:,1), collapse3D(:,3), 10, 'b', 'filled', 'DisplayName', '崩落点（二维）');
%     % 绘制二维破裂点
%     scatter(fracture3D(:,1), fracture3D(:,3), 10, 'g', 'filled', 'DisplayName', '破裂点（二维）');
%     
%     xlabel('X');
%     xlabel('X');
%     ylabel('Z');
%     legend show;
%     title('三心圆拱（二维视图）');
%     grid on;
%     hold off;
%     
%     %     计算 y 轴范围
%     ymin = min(currentPoints(:, 2));
%     ymax = ymin + 500;  % y轴长度为500
%     
%     % 将二维点扩展为三维点 (沿y轴延伸)
%     y_vals = linspace(ymin, ymax, 100);  % y轴方向的点
%     [X_MCX, Y_MCX] = meshgrid(arc_MCX_x, y_vals);
%     [Z_MCX, ~] = meshgrid(arc_MCX_z, y_vals);
%     
%     [X_AM, Y_AM] = meshgrid(arc_AM_x, y_vals);
%     [Z_AM, ~] = meshgrid(arc_AM_z, y_vals);
%     
%     [X_MKF, Y_MKF] = meshgrid(arc_MKF_x, y_vals);
%     [Z_MKF, ~] = meshgrid(arc_MCX_z, y_vals);
%     
%     [X_AM_mirror, Y_AM_mirror] = meshgrid(arc_AM_mirror_x, y_vals);
%     [Z_AM_mirror, ~] = meshgrid(arc_AM_z, y_vals);
%     % 绘制三维图形
%     figure; 
%     hold on; 
%     % 绘制崩落点
%     scatter3(collapse3D(:,1), collapse3D(:,2), collapse3D(:,3), ...
%         10, 'b', 'filled', 'DisplayName', '崩落点');
%     % 绘制破裂点
%     scatter3(fracture3D(:,1), fracture3D(:,2), fracture3D(:,3), ...
%        10, 'g', 'filled', 'DisplayName', '破裂点');
%     % 绘制三维三心圆拱
%     surf(X_MCX, Y_MCX, Z_MCX, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%     surf(X_AM, Y_AM, Z_AM, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%     surf(X_MKF, Y_MKF, Z_MKF, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%     surf(X_AM_mirror, Y_AM_mirror, Z_AM_mirror, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%     
%     % 图形设置
%     xlabel('X');
%     ylabel('Y');
%     zlabel('Z');
%     legend show; % 显示图例
%     title('三心圆拱筛选结果（三维视图）');
%     grid on;
%     hold off; 
% end
% 
% % 辅助函数：判断点是否在左侧/右侧圆弧内
% function inArc = check_arc(pt, O, P_start, P_end, R)
%     vec_OP = pt - O;
%     dist = norm(vec_OP);
%     if dist > R
%         inArc = false;
%         return;
%     end
%     
%     theta = adjust_angle(atan2(vec_OP(2), vec_OP(1)));
%     theta_start = adjust_angle(atan2(P_start(2)-O(2), P_start(1)-O(1)));
%     theta_end = adjust_angle(atan2(P_end(2)-O(2), P_end(1)-O(1)));
%     
%     % 处理角度跨越2π的情况
%     if theta_end < theta_start
%         inArc = (theta >= theta_start) || (theta <= theta_end);
%     else
%         inArc = (theta >= theta_start) && (theta <= theta_end);
%     end
% end
% 
% % 辅助函数：判断点是否在中间圆弧内
% function inArc = check_arc_mid(pt, O, theta_start, theta_end, R)
%     vec_OP = pt - O;
%     dist = norm(vec_OP);
%     if dist > R
%         inArc = false;
%         return;
%     end
%     
%     theta = adjust_angle(atan2(vec_OP(2), vec_OP(1)));
%     
%     if theta_end < theta_start
%         inArc = (theta >= theta_start) || (theta <= theta_end);
%     else
%         inArc = (theta >= theta_start) && (theta <= theta_end);
%     end
% end
% 
% % 角度调整函数
% function theta = adjust_angle(theta)
%     theta = mod(theta + 2*pi, 2*pi);
% end  

%% 原版合并旋转版本

% function [collapse3D, fracture3D, arcRadiusLeft, M1, arcAnglemid, arcRadiusmid] = process_san1(currentPoints, L, h, D_x, D_z)
%     % 处理单组数据的崩落点筛选（三心圆拱约束）
%     % 输入输出参数同原函数
%     
%     % 矩形 AFEG 的顶点坐标（二维XZ平面）
%     A = [D_x - L/2, D_z];      % 左下角
%     F = [D_x + L/2, D_z];      % 右下角
%     E = [D_x + L/2, D_z + h];  % 右上角
%     G = [D_x - L/2, D_z + h];  % 左上角
%     
%     % 计算角平分线交点 M1
%     C = [D_x, D_z + h];  % 上顶点
%     D = [D_x, D_z];      % AF中点
%     
%     % 向量计算
%     vec_CG = G - C;
%     vec_CA = A - C;
%     vec_AG = G - A;
%     vec_AC = C - A;
%     
%     % 单位向量
%     unit_CG = vec_CG / norm(vec_CG);
%     unit_CA = vec_CA / norm(vec_CA);
%     unit_AG = vec_AG / norm(vec_AG);
%     unit_AC = vec_AC / norm(vec_AC);
%     
%     % 角平分线方向
%     dir_C_bisector = unit_CG + unit_CA;  % C点的角平分线方向向量
%     dir_A_bisector = unit_AG + unit_AC;  % A点的角平分线方向向量
%     
%     % 解联立方程求交点 M1
%     % 参数方程: C + t*dir_C_bisector = A + s*dir_A_bisector
%     % 转换为线性方程组求解t和s
%     A_matrix = [dir_C_bisector(1), -dir_A_bisector(1); 
%                 dir_C_bisector(2), -dir_A_bisector(2)];
%     b_vector = [A(1) - C(1); 
%                 A(2) - C(2)];
%     solution = A_matrix \ b_vector;
%     t = solution(1);
%     M1 = C + t * dir_C_bisector;  % 交点坐标
%     
%     % 计算中间圆弧圆心O（过M1作AC的垂线与CD的交点）
%     AC_slope = (C(2) - A(2)) / (C(1) - A(1));
%     if AC_slope == 0
%         % AC水平，垂线为竖直线x=M1(1)
%         O = [D_x, M1(2)];  % CD线为x=D_x
%     else
%         perp_slope = -1 / AC_slope;
%         % 垂线方程: z = perp_slope*(x - M1(1)) + M1(2)
%         % CD线方程: x = D_x
%         O_x = D_x;
%         O_z = perp_slope * (D_x - M1(1)) + M1(2);
%         O = [O_x, O_z];
%     end
%     
%     % 计算左侧圆弧圆心O1（MO与AF的交点）
%     if O(1) == M1(1)
%         O1 = [M1(1), A(2)];  % 垂直线
%     else
%         MO_slope = (O(2) - M1(2)) / (O(1) - M1(1));
%         O1_x = M1(1) + (A(2) - M1(2)) / MO_slope;
%         O1 = [O1_x, A(2)];
%     end
% 
%     % 绘制圆弧 MCX (以 O 为圆心)
%     theta_M = atan2(M1(2)-O(2), M1(1)-O(1));
%     theta_C = atan2(C(2)-O(2), C(1)-O(1));
%     theta = linspace(theta_M, theta_C, 100);
%     arc_MCX_x = O(1) + norm(O-M1)*cos(theta);
%     arc_MCX_z = O(2) + norm(O-M1)*sin(theta);
%     
%     % 绘制圆弧 AM (以 O1 为圆心)
%     theta_A = atan2(A(2)-O1(2), A(1)-O1(1));
%     theta_M_O1 = atan2(M1(2)-O1(2), M1(1)-O1(1));
%     theta_AM = linspace(theta_A, theta_M_O1, 100);
%     arc_AM_x = O1(1) + norm(O1-A)*cos(theta_AM);
%     arc_AM_z = O1(2) + norm(O1-A)*sin(theta_AM);
%     
%     % 镜像右侧
%     arc_MKF_x = 2*D_x - arc_MCX_x;  % 关于 x = D_x 镜像
%     arc_AM_mirror_x = 2*D_x - arc_AM_x;
%     
%     % 右侧对称点
%     O2 = [2*D_x - O1(1), O1(2)];
%     M2 = [2*D_x - M1(1), M1(2)];
%     
%     % 计算圆弧参数
%     arcRadiusLeft = norm(O1 - A);   % 左侧圆弧半径
%     arcRadiusmid = norm(O - C);     % 中间圆弧半径
%     disp(arcRadiusmid);
%     disp(arcRadiusLeft);
%     % 角度范围处理函数
%     adjust_angle = @(theta) mod(theta + 2*pi, 2*pi);
%     
%     % 中间圆弧角度范围
%     theta_M = adjust_angle(atan2(M1(2)-O(2), M1(1)-O(1)));
%     theta_C = adjust_angle(atan2(C(2)-O(2), C(1)-O(1)));
%     theta_M2 = adjust_angle(atan2(M2(2)-O(2), M2(1)-O(1)));
%     arcAnglemid = adjust_angle(theta_C - theta_M);
%     % 按顺时针顺序合并所有圆弧段（左侧圆弧 -> 中间圆弧 -> 右侧圆弧 -> 镜像左侧圆弧）
%     combined_x = [arc_AM_x, arc_MCX_x, fliplr(arc_MKF_x), fliplr(arc_AM_mirror_x)];
%     combined_z = [arc_AM_z, arc_MCX_z, fliplr(arc_MCX_z), fliplr(arc_AM_z)];
%     
%     % 去重并确保首尾闭合
%     [combined_x, combined_z] = remove_duplicate_points(combined_x, combined_z);
%     if ~isequal([combined_x(1), combined_z(1)], [combined_x(end), combined_z(end)])
%         combined_x(end+1) = combined_x(1);
%         combined_z(end+1) = combined_z(1);
%     end
%     
%     % 参数化闭合曲线（样条插值）
%     t = 1:length(combined_x);
%     tt = linspace(1, length(combined_x), 1000); % 增加插值密度
%     spline_x = spline(t, combined_x, tt);
%     spline_z = spline(t, combined_z, tt);
%     
%     % ========== 三维曲面生成 ==========
%     % 沿Y轴拉伸生成三维曲面
%     ymin = min(currentPoints(:, 2));
%     ymax = ymin + 500;
%     y_vals = linspace(ymin, ymax, 100);
%     
%     [X, Y] = meshgrid(spline_x, y_vals);
%     [Z, ~] = meshgrid(spline_z, y_vals);
%     % 初始化输出
%     collapse3D = [];
%     fracture3D = [];
%     
%     for i = 1:size(currentPoints, 1)
%         pt = currentPoints(i, [1,3]); % 提取XZ坐标
%         x = pt(1); z = pt(2);
%         
%         % 1. 排除矩形外点
%         if x < A(1) || x > F(1) || z < D_z || z > D_z + h
%             fracture3D = [fracture3D; currentPoints(i, :)];
%             continue;
%         end
%         
%         % 2. 判断点是否在任一圆弧内
%         inArc = false;
%         
%         % 左侧圆弧判断
%         if check_arc(pt, O1, A, M1, arcRadiusLeft)
%             inArc = true;
%         end      
%         % 中间圆弧判断
%         if check_arc_mid(pt, O, theta_M, theta_C, arcRadiusmid)
%             inArc = true;
%         end
%         if check_arc_mid(pt, O, theta_C, theta_M2, arcRadiusmid)
%             inArc = true;
%         end
%         % 右侧圆弧判断
%         if check_arc(pt, O2, F, M2, arcRadiusLeft)
%             inArc = true;
%         end
%         
%         if inArc
%             collapse3D = [collapse3D; currentPoints(i, :)];
%         else
%             fracture3D = [fracture3D; currentPoints(i, :)];
%         end
%     end
%     
%     % 去除重复点
%     [collapse3D, ~] = unique(collapse3D, 'rows', 'stable');
%     [fracture3D, ~] = unique(fracture3D, 'rows', 'stable');
%     
%     
%     % ========== 可视化部分 ==========
%     % 二维视图显示合并后的轮廓
%     figure;
%     plot(spline_x, spline_z, 'r-', 'LineWidth', 2, 'DisplayName', '合并轮廓');
%     hold on;
%     scatter(collapse3D(:,1), collapse3D(:,3), 10, 'b', 'filled', 'DisplayName', '崩落点');
%     scatter(fracture3D(:,1), fracture3D(:,3), 10, 'g', 'filled', 'DisplayName', '破裂点');
%     title('合并后的二维三心圆拱'); 
%     axis equal; legend; hold off;
%     
%     % 三维曲面绘制
%     figure;
%     surf(X, Y, Z, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%     hold on;
%     scatter3(collapse3D(:,1), collapse3D(:,2), collapse3D(:,3), 10, 'b', 'filled');
%     scatter3(fracture3D(:,1), fracture3D(:,2), fracture3D(:,3), 10, 'g', 'filled');
%     title('三维崩落分析曲面'); xlabel('X'); ylabel('Y'); zlabel('Z');
%     axis equal; view(3); grid on;
% end
function [collapse3D, fracture3D, arcRadiusLeft, M1, arcAnglemid, arcRadiusmid] = process_san1(currentPoints, L, h, D_x, D_z)
    % 处理单组数据的崩落点筛选（三心圆拱约束）
    % 输入输出参数同原函数
    % 矩形 AFEG 的顶点坐标（二维XZ平面）
    A = [D_x - L/2, D_z];      % 左下角
    F = [D_x + L/2, D_z];      % 右下角
    E = [D_x + L/2, D_z + h];  % 右上角
    G = [D_x - L/2, D_z + h];  % 左上角
    
    % 计算角平分线交点 M1
    C = [D_x, D_z + h];  % 上顶点
    D = [D_x, D_z];      % AF中点
    
    % 向量计算
    vec_CG = G - C;
    vec_CA = A - C;
    vec_AG = G - A;
    vec_AC = C - A;
    
    % 单位向量
    unit_CG = vec_CG / norm(vec_CG);
    unit_CA = vec_CA / norm(vec_CA);
    unit_AG = vec_AG / norm(vec_AG);
    unit_AC = vec_AC / norm(vec_AC);
    
    % 角平分线方向
    dir_C_bisector = unit_CG + unit_CA;  % C点的角平分线方向向量
    dir_A_bisector = unit_AG + unit_AC;  % A点的角平分线方向向量
    
    % 解联立方程求交点 M1
    % 参数方程: C + t*dir_C_bisector = A + s*dir_A_bisector
    % 转换为线性方程组求解t和s
    A_matrix = [dir_C_bisector(1), -dir_A_bisector(1); 
                dir_C_bisector(2), -dir_A_bisector(2)];
    b_vector = [A(1) - C(1); 
                A(2) - C(2)];
    solution = A_matrix \ b_vector;
    t = solution(1);
    M1 = C + t * dir_C_bisector;  % 交点坐标
    
    % 计算中间圆弧圆心O（过M1作AC的垂线与CD的交点）
    AC_slope = (C(2) - A(2)) / (C(1) - A(1));
    if AC_slope == 0
        % AC水平，垂线为竖直线x=M1(1)
        O = [D_x, M1(2)];  % CD线为x=D_x
    else
        perp_slope = -1 / AC_slope;
        % 垂线方程: z = perp_slope*(x - M1(1)) + M1(2)
        % CD线方程: x = D_x
        O_x = D_x;
        O_z = perp_slope * (D_x - M1(1)) + M1(2);
        O = [O_x, O_z];
    end
    
    % 计算左侧圆弧圆心O1（MO与AF的交点）
    if O(1) == M1(1)
        O1 = [M1(1), A(2)];  % 垂直线
    else
        MO_slope = (O(2) - M1(2)) / (O(1) - M1(1));
        O1_x = M1(1) + (A(2) - M1(2)) / MO_slope;
        O1 = [O1_x, A(2)];
    end

    % 绘制圆弧 MCX (以 O 为圆心)
    theta_M = atan2(M1(2)-O(2), M1(1)-O(1));
    theta_C = atan2(C(2)-O(2), C(1)-O(1));
    theta = linspace(theta_M, theta_C, 100);
    arc_MCX_x = O(1) + norm(O-M1)*cos(theta);
    arc_MCX_z = O(2) + norm(O-M1)*sin(theta);
    
    % 绘制圆弧 AM (以 O1 为圆心)
    theta_A = atan2(A(2)-O1(2), A(1)-O1(1));
    theta_M_O1 = atan2(M1(2)-O1(2), M1(1)-O1(1));
    theta_AM = linspace(theta_A, theta_M_O1, 100);
    arc_AM_x = O1(1) + norm(O1-A)*cos(theta_AM);
    arc_AM_z = O1(2) + norm(O1-A)*sin(theta_AM);
    
    % 镜像右侧
    arc_MKF_x = 2*D_x - arc_MCX_x;  % 关于 x = D_x 镜像
    arc_AM_mirror_x = 2*D_x - arc_AM_x;
    
    % 右侧对称点
    O2 = [2*D_x - O1(1), O1(2)];
    M2 = [2*D_x - M1(1), M1(2)];
    
    % 计算圆弧参数
    arcRadiusLeft = norm(O1 - A);   % 左侧圆弧半径
    arcRadiusmid = norm(O - C);     % 中间圆弧半径

    % 角度范围处理函数
    adjust_angle = @(theta) mod(theta + 2*pi, 2*pi);
    
    % 中间圆弧角度范围
    theta_M = adjust_angle(atan2(M1(2)-O(2), M1(1)-O(1)));
    theta_C = adjust_angle(atan2(C(2)-O(2), C(1)-O(1)));
    theta_M2 = adjust_angle(atan2(M2(2)-O(2), M2(1)-O(1)));
    arcAnglemid = adjust_angle(theta_C - theta_M);
    
    
    % 按顺时针顺序合并所有圆弧段（左侧圆弧 -> 中间圆弧 -> 右侧圆弧 -> 镜像左侧圆弧）
    combined_x = [arc_AM_x, arc_MCX_x, fliplr(arc_MKF_x), fliplr(arc_AM_mirror_x)];
    combined_z = [arc_AM_z, arc_MCX_z, fliplr(arc_MCX_z), fliplr(arc_AM_z)];
    
    % 去重并确保首尾闭合
    [combined_x, combined_z] = remove_duplicate_points(combined_x, combined_z);
    if ~isequal([combined_x(1), combined_z(1)], [combined_x(end), combined_z(end)])
        combined_x(end+1) = combined_x(1);
        combined_z(end+1) = combined_z(1);
    end
    
    % 参数化闭合曲线（样条插值）
    t = 1:length(combined_x);
    tt = linspace(1, length(combined_x), 1000); % 增加插值密度
    spline_x = spline(t, combined_x, tt);
    spline_z = spline(t, combined_z, tt);

    % 计算PCA旋转角度
    if size(currentPoints,1) >= 2
        data = currentPoints(:,[1,3]);
        [coeff, ~, ~] = pca(data);
        principal_dir = coeff(:,1);
        theta = atan2(principal_dir(2), principal_dir(1));
    else
        theta = 0;
    end
    
    % 绕中心(D_x,D_z)旋转轮廓
    cx = D_x;
    cz = D_z;
    dx = spline_x - cx;
    dz = spline_z - cz;
    spline_x_rot = cx + dx*cos(theta) - dz*sin(theta);
    spline_z_rot = cz + dx*sin(theta) + dz*cos(theta);

    % 闭合旋转后的轮廓
    if ~isequal([spline_x_rot(end), spline_z_rot(end)], [spline_x_rot(1), spline_z_rot(1)])
        spline_x_rot(end+1) = spline_x_rot(1);
        spline_z_rot(end+1) = spline_z_rot(1);
    end

    % ========== 三维曲面生成 ==========
    % 沿Y轴拉伸生成三维曲面
    ymin = min(currentPoints(:, 2));
    ymax = ymin + 500;
    y_vals = linspace(ymin, ymax, 100);
    
    [X, Y] = meshgrid(spline_x, y_vals);
    [Z, ~] = meshgrid(spline_z, y_vals);
    % 筛选点：判断是否在旋转后的多边形内
    collapse3D = [];
    fracture3D = [];
    for i = 1:size(currentPoints,1)
        x = currentPoints(i,1);
        z = currentPoints(i,3);
        if inpolygon(x, z, spline_x_rot, spline_z_rot)
            collapse3D = [collapse3D; currentPoints(i,:)];
        else
            fracture3D = [fracture3D; currentPoints(i,:)];
        end
    end

    % 去重
    [collapse3D, ~] = unique(collapse3D, 'rows', 'stable');
    [fracture3D, ~] = unique(fracture3D, 'rows', 'stable');

    % 可视化旋转后的轮廓
    figure;
%     plot(spline_x_rot, spline_z_rot, 'r-', 'LineWidth', 2);
    plot(combined_x, combined_z, 'r-', 'LineWidth', 2);
    hold on;
    scatter(collapse3D(:,1), collapse3D(:,3), 10, 'b', 'filled');
    scatter(fracture3D(:,1), fracture3D(:,3), 10, 'g', 'filled');
    title('旋转后的二维三心圆拱'); axis equal; hold off;

%     % 三维可视化（沿Y轴拉伸）
    [X, Y] = meshgrid(spline_x_rot, linspace(min(currentPoints(:,2)), max(currentPoints(:,2)), 100));
    Z = repmat(spline_z_rot, size(X,1),1);

    figure;
    surf(X, Y, Z, 'FaceAlpha',0.3, 'EdgeColor','none');
    hold on;
    scatter3(collapse3D(:,1), collapse3D(:,2), collapse3D(:,3), 10, 'b', 'filled');
    scatter3(fracture3D(:,1), fracture3D(:,2), fracture3D(:,3), 10, 'g', 'filled');
    xlabel("X");
    ylabel("Y");
    zlabel("Z");
    title('三维旋转崩落分析曲面'); axis equal; view(3); grid on;
end 
% 辅助函数：去除相邻重复点
function [x_unique, z_unique] = remove_duplicate_points(x, z)
    tolerance = 1e-6;
    keep_idx = [true, ~(abs(diff(x)) < tolerance & abs(diff(z)) < tolerance)];
    x_unique = x(keep_idx);
    z_unique = z(keep_idx);
end
% 辅助函数：判断点是否在左侧/右侧圆弧内
function inArc = check_arc(pt, O, P_start, P_end, R)
    vec_OP = pt - O;
    dist = norm(vec_OP);
    if dist > R
        inArc = false;
        return;
    end
    
    theta = adjust_angle(atan2(vec_OP(2), vec_OP(1)));
    theta_start = adjust_angle(atan2(P_start(2)-O(2), P_start(1)-O(1)));
    theta_end = adjust_angle(atan2(P_end(2)-O(2), P_end(1)-O(1)));
    
    % 处理角度跨越2π的情况
    if theta_end < theta_start
        inArc = (theta >= theta_start) || (theta <= theta_end);
    else
        inArc = (theta >= theta_start) && (theta <= theta_end);
    end
end

% 辅助函数：判断点是否在中间圆弧内
function inArc = check_arc_mid(pt, O, theta_start, theta_end, R)
    vec_OP = pt - O;
    dist = norm(vec_OP);
    if dist > R
        inArc = false;
        return;
    end
    
    theta = adjust_angle(atan2(vec_OP(2), vec_OP(1)));
    
    if theta_end < theta_start
        inArc = (theta >= theta_start) || (theta <= theta_end);
    else
        inArc = (theta >= theta_start) && (theta <= theta_end);
    end
end

% 角度调整函数
function theta = adjust_angle(theta)
    theta = mod(theta + 2*pi, 2*pi);
end  
%% 旋转版

% function [collapse3D, fracture3D, arcRadiusLeft, M1, arcAnglemid, arcRadiusmid] = process_san1(currentPoints, L, h, D_x, D_z, theta)
%     % 处理单组数据的崩落点筛选（支持三维旋转的三心圆拱约束）
%     % 输入参数:
%     %   currentPoints - 原始点云数据（N×3矩阵）
%     %   L - 矩形底边长度
%     %   h - 矩形高度
%     %   D_x, D_z - 矩形中心坐标
%     %   theta - 绕点 A 的旋转角度（弧度）
%     % 输出参数:
%     %   collapse3D - 崩落点集合
%     %   fracture3D - 破裂点集合
% 
%     if nargin < 6
%         theta = 0; % 默认不旋转
%     end
% 
%     % 检查 theta 是否为标量
%     if ~isscalar(theta)
%         error('theta 必须是一个标量值，表示绕点 A 的旋转角度（弧度）。');
%     end
% 
%     % 打印 theta 的初始值和维度
%     fprintf('初始 theta 的值: %f\n', theta);
%     fprintf('初始 theta 的维度: %s\n', mat2str(size(theta)));
% 
%     % ================== 旋转坐标系（二维部分）================== 
%     rotate_point = @(p, center, theta) [...
%         (p(1)-center(1))*cos(theta) - (p(2)-center(2))*sin(theta) + center(1),...
%         (p(1)-center(1))*sin(theta) + (p(2)-center(2))*cos(theta) + center(2)];
%     
%     center = [D_x, D_z];
%     A_orig = [D_x - L/2, D_z];
%     F_orig = [D_x + L/2, D_z];
%     E_orig = [D_x + L/2, D_z + h];
%     G_orig = [D_x - L/2, D_z + h];
%     
%     % 应用二维旋转
%     A = rotate_point(A_orig, A_orig, theta);
%     F = rotate_point(F_orig, A_orig, theta);
%     E = rotate_point(E_orig, A_orig, theta);
%     G = rotate_point(G_orig, A_orig, theta);
%     C = rotate_point([D_x, D_z + h], A_orig, theta);
% 
%     % ================== 几何计算（基于旋转后坐标）================== 
%     % 计算角平分线交点 M1
%  
%     % 向量计算（使用旋转后的坐标）
%     vec_CG = G - C;
%     vec_CA = A - C;
%     vec_AG = G - A;
%     vec_AC = C - A;
%     
%     % 单位向量
%     unit_CG = vec_CG / norm(vec_CG);
%     unit_CA = vec_CA / norm(vec_CA);
%     unit_AG = vec_AG / norm(vec_AG);
%     unit_AC = vec_AC / norm(vec_AC);
%     
%     % 角平分线方向
%     dir_C_bisector = unit_CG + unit_CA;
%     dir_A_bisector = unit_AG + unit_AC;
%     
%     % 求解交点 M1
%     A_matrix = [dir_C_bisector(1), -dir_A_bisector(1); 
%                 dir_C_bisector(2), -dir_A_bisector(2)];
%     b_vector = [A(1) - C(1); 
%                 A(2) - C(2)];
%     solution = A_matrix \ b_vector;
%     t = solution(1);
%     M1 = C + t * dir_C_bisector;
%     % 计算中间圆弧圆心O（过M1作AC的垂线与CD的交点）
%     AC_slope = (C(2) - A(2)) / (C(1) - A(1));
%     if AC_slope == 0
%         O = [D_x, M1(2)];
%     else
%         perp_slope = -1 / AC_slope;
%         O_x = D_x;
%         O_z = perp_slope * (D_x - M1(1)) + M1(2);
%         O = [O_x, O_z];
%     end
%     
%     % 计算左侧圆弧圆心O1
%     if O(1) == M1(1)
%         O1 = [M1(1), A(2)];
%     else
%         MO_slope = (O(2) - M1(2)) / (O(1) - M1(1));
%         O1_x = M1(1) + (A(2) - M1(2)) / MO_slope;
%         O1 = [O1_x, A(2)];
%     end
% 
%     % ============== 生成圆弧坐标（已包含旋转）============== 
%     % 中间圆弧
%     theta_M = atan2(M1(2)-O(2), M1(1)-O(1));
%     theta_C = atan2(C(2)-O(2), C(1)-O(1));
%     % 使用新变量存储向量
%     theta_arc = linspace(theta_M, theta_C, 100);
%     arc_MCX_x = O(1) + norm(O-M1)*cos(theta_arc);
%     arc_MCX_z = O(2) + norm(O-M1)*sin(theta_arc);
%     
%     % 左侧圆弧
%     theta_A = atan2(A(2)-O1(2), A(1)-O1(1));
%     theta_M_O1 = atan2(M1(2)-O1(2), M1(1)-O1(1));
%     theta_arc_AM = linspace(theta_A, theta_M_O1, 100);
%     arc_AM_x = O1(1) + norm(O1-A)*cos(theta_arc_AM);
%     arc_AM_z = O1(2) + norm(O1-A)*sin(theta_arc_AM);
%     
%     % 镜像右侧
%     arc_MKF_x = 2*D_x - arc_MCX_x;
%     arc_AM_mirror_x = 2*D_x - arc_AM_x;
%     
%     % 右侧对称点
%     O2 = [2*D_x - O1(1), O1(2)];
%     M2 = [2*D_x - M1(1), M1(2)];
%     
%     % 计算圆弧参数
%     arcRadiusLeft = norm(O1 - A);   % 左侧圆弧半径
%     arcRadiusmid = norm(O - C);     % 中间圆弧半径
%     disp(arcRadiusmid);
%     disp(arcRadiusLeft);
%     % 角度范围处理函数
%     adjust_angle = @(theta) mod(theta + 2*pi, 2*pi);
%     
%     % 中间圆弧角度范围
%     theta_M = adjust_angle(atan2(M1(2)-O(2), M1(1)-O(1)));
%     theta_C = adjust_angle(atan2(C(2)-O(2), C(1)-O(1)));
%     theta_M2 = adjust_angle(atan2(M2(2)-O(2), M2(1)-O(1)));
%     arcAnglemid = adjust_angle(theta_C - theta_M);
%     % 初始化输出
%     collapse3D = [];
%     fracture3D = [];
% 
%     % 计算 y 轴范围
%     ymin = min(currentPoints(:, 2));
%     ymax = ymin + 500;  % y轴长度为500
%     
%     % 将二维点扩展为三维点 (沿y轴延伸)
%     y_vals = linspace(ymin, ymax, 100);  % y轴方向的点
%     [X_MCX, Y_MCX] = meshgrid(arc_MCX_x, y_vals);
%     [Z_MCX, ~] = meshgrid(arc_MCX_z, y_vals);
%     
%     [X_AM, Y_AM] = meshgrid(arc_AM_x, y_vals);
%     [Z_AM, ~] = meshgrid(arc_AM_z, y_vals);
%     
%     [X_MKF, Y_MKF] = meshgrid(arc_MKF_x, y_vals);
%     [Z_MKF, ~] = meshgrid(arc_MCX_z, y_vals);
%     
%     [X_AM_mirror, Y_AM_mirror] = meshgrid(arc_AM_mirror_x, y_vals);
%     [Z_AM_mirror, ~] = meshgrid(arc_AM_z, y_vals);
%     
%     % 打印 theta 在构建旋转矩阵前的值和维度
%     fprintf('构建旋转矩阵前 theta 的值: %f\n', theta);
%     fprintf('构建旋转矩阵前 theta 的维度: %s\n', mat2str(size(theta)));
% 
%     % 将二维圆弧点扩展到三维
%     points_MCX = [X_MCX(:), Y_MCX(:), Z_MCX(:)];
%     points_AM = [X_AM(:), Y_AM(:), Z_AM(:)];
%     points_MKF = [X_MKF(:), Y_MKF(:), Z_MKF(:)];
%     points_AM_mirror = [X_AM_mirror(:), Y_AM_mirror(:), Z_AM_mirror(:)];
% 
%     % 构建围绕点 A 的旋转矩阵
%     rotation_matrix = [cos(theta), 0, sin(theta);
%                        0, 1, 0;
%                       -sin(theta), 0, cos(theta)];
% 
%     % 围绕点 A 旋转
%     A_3D = [A_orig(1), 0, A_orig(2)]; % 将 A 扩展到三维
%     points_MCX_rot = rotate_around_point(points_MCX, A_3D, rotation_matrix);
%     points_AM_rot = rotate_around_point(points_AM, A_3D, rotation_matrix);
%     points_MKF_rot = rotate_around_point(points_MKF, A_3D, rotation_matrix);
%     points_AM_mirror_rot = rotate_around_point(points_AM_mirror, A_3D, rotation_matrix);
%     
%     % 确保旋转后的点的顺序正确（这里简单按行排序示例，具体根据实际情况调整）
%     [~, idx_MCX] = sortrows(points_MCX_rot);
%     points_MCX_rot = points_MCX_rot(idx_MCX, :);
%     [~, idx_AM] = sortrows(points_AM_rot);
%     points_AM_rot = points_AM_rot(idx_AM, :);
%     [~, idx_MKF] = sortrows(points_MKF_rot);
%     points_MKF_rot = points_MKF_rot(idx_MKF, :);
%     [~, idx_AM_mirror] = sortrows(points_AM_mirror_rot);
%     points_AM_mirror_rot = points_AM_mirror_rot(idx_AM_mirror, :);
%     
%     % 重塑旋转后的点
%     X_MCX_rot = reshape(points_MCX_rot(:,1), size(X_MCX));
%     Y_MCX_rot = reshape(points_MCX_rot(:,2), size(Y_MCX));
%     Z_MCX_rot = reshape(points_MCX_rot(:,3), size(Z_MCX));
%     
%     X_AM_rot = reshape(points_AM_rot(:,1), size(X_AM));
%     Y_AM_rot = reshape(points_AM_rot(:,2), size(Y_AM));
%     Z_AM_rot = reshape(points_AM_rot(:,3), size(Z_AM));
%     
%     X_MKF_rot = reshape(points_MKF_rot(:,1), size(X_MKF));
%     Y_MKF_rot = reshape(points_MKF_rot(:,2), size(Y_MKF));
%     Z_MKF_rot = reshape(points_MKF_rot(:,3), size(Z_MKF));
%     
%     X_AM_mirror_rot = reshape(points_AM_mirror_rot(:,1), size(X_AM_mirror));
%     Y_AM_mirror_rot = reshape(points_AM_mirror_rot(:,2), size(Y_AM_mirror));
%     Z_AM_mirror_rot = reshape(points_AM_mirror_rot(:,3), size(Z_AM_mirror));
% 
%     for i = 1:size(currentPoints, 1)
%         pt = currentPoints(i, [1,3]); % 提取XZ坐标
%         x = pt(1); z = pt(2);
%         
%         % 1. 排除矩形外点
%         if x < A_orig(1) || x > F_orig(1) || z < D_z || z > D_z + h
%             fracture3D = [fracture3D; currentPoints(i, :)];
%             continue;
%         end
%         
%         % 2. 判断点是否在任一圆弧内
%         inArc = false;
%         
%         % 转换旋转后的圆心和端点到二维
%         O1_2D = points_MCX_rot(1, [1,3]);
%         A_2D = A_orig;
%         M1_2D = points_MCX_rot(1, [1,3]);
%         O_2D = points_MCX_rot(1, [1,3]);
%         M2_2D = points_MCX_rot(1, [1,3]);
%         F_2D = F_orig;
% 
%         % 左侧圆弧判断
%         if check_arc(pt, O1_2D, A_2D, M1_2D, arcRadiusLeft)
%             inArc = true;
%         end      
%         % 中间圆弧判断
%         if check_arc_mid(pt, O_2D, theta_M, theta_C, arcRadiusmid)
%             inArc = true;
%         end
%         if check_arc_mid(pt, O_2D, theta_C, theta_M2, arcRadiusmid)
%             inArc = true;
%         end
%         % 右侧圆弧判断
%         if check_arc(pt, O1_2D, F_2D, M2_2D, arcRadiusLeft)
%             inArc = true;
%         end
%         
%         if inArc
%             collapse3D = [collapse3D; currentPoints(i, :)];
%         else
%             fracture3D = [fracture3D; currentPoints(i, :)];
%         end
%     end
%     % 绘制二维旋转后的三心圆拱
%     figure;
%     hold on;
%     plot(arc_MCX_x, arc_MCX_z, 'b', 'DisplayName', '中间圆弧');
%     plot(arc_AM_x, arc_AM_z, 'r', 'DisplayName', '左侧圆弧');
%     plot(arc_MKF_x, arc_MCX_z, 'g', 'DisplayName', '右侧中间圆弧镜像');
%     plot(arc_AM_mirror_x, arc_AM_z, 'm', 'DisplayName', '右侧左侧圆弧镜像');
%     plot([A(1), F(1), F(1), A(1), A(1)], [A(2), F(2), E(2), G(2), A(2)], 'k--', 'DisplayName', '矩形边界');
%     scatter([A(1), F(1), E(1), G(1), C(1), M1(1), O(1), O1(1), O2(1), M2(1)], [A(2), F(2), E(2), G(2), C(2), M1(2), O(2), O1(2), O2(2), M2(2)], 'ko', 'DisplayName', '关键点');
%     legend show;
%     xlabel('X');
%     ylabel('Z');
%     title('二维旋转后的三心圆拱');
%     grid on;
%     hold off;
%     
%     % 绘制三维图形
%     figure; 
%     hold on; 
%     % 绘制崩落点
%     if ~isempty(collapse3D)
%         scatter3(collapse3D(:,1), collapse3D(:,2), collapse3D(:,3), ...
%             10, 'b', 'filled', 'DisplayName', '崩落点');
%     end
%     % 绘制破裂点
%     if ~isempty(fracture3D)
%         scatter3(fracture3D(:,1), fracture3D(:,2), fracture3D(:,3), ...
%            10, 'g', 'filled', 'DisplayName', '破裂点');
%     end
%     % 绘制旋转后的三维三心圆拱
%     surf(X_MCX_rot, Y_MCX_rot, Z_MCX_rot, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%     surf(X_AM_rot, Y_AM_rot, Z_AM_rot, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%     surf(X_MKF_rot, Y_MKF_rot, Z_MKF_rot, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%     surf(X_AM_mirror_rot, Y_AM_mirror_rot, Z_AM_mirror_rot, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%     
%     % 图形设置
%     xlabel('X');
%     ylabel('Y');
%     zlabel('Z');
%     legend show; % 显示图例
%     title('三心圆拱筛选结果（三维视图）');
%     grid on;
%     hold off; 
% end
% 
% % 围绕给定点旋转点集的函数
% function rotated_points = rotate_around_point(points, center, rotation_matrix)
%     % 平移到原点
%     translated_points = points - repmat(center, size(points, 1), 1);
%     % 旋转
%     rotated_points = (rotation_matrix * translated_points')';
%     % 平移回去
%     rotated_points = rotated_points + repmat(center, size(rotated_points, 1), 1);
% end
%  
% % 辅助函数：判断点是否在左侧/右侧圆弧内
% function inArc = check_arc(pt, O, P_start, P_end, R)
%     vec_OP = pt - O;
%     dist = norm(vec_OP);
%     if dist > R
%         inArc = false;
%         return;
%     end
%     
%     theta = adjust_angle(atan2(vec_OP(2), vec_OP(1)));
%     theta_start = adjust_angle(atan2(P_start(2)-O(2), P_start(1)-O(1)));
%     theta_end = adjust_angle(atan2(P_end(2)-O(2), P_end(1)-O(1)));
%     
%     % 处理角度跨越2π的情况
%     if theta_end < theta_start
%         inArc = (theta >= theta_start) || (theta <= theta_end);
%     else
%         inArc = (theta >= theta_start) && (theta <= theta_end);
%     end
% end
% 
% % 辅助函数：判断点是否在中间圆弧内
% function inArc = check_arc_mid(pt, O, theta_start, theta_end, R)
%     vec_OP = pt - O;
%     dist = norm(vec_OP);
%     if dist > R
%         inArc = false;
%         return;
%     end
%     
%     theta = adjust_angle(atan2(vec_OP(2), vec_OP(1)));
%     
%     if theta_end < theta_start
%         inArc = (theta >= theta_start) || (theta <= theta_end);
%     else
%         inArc = (theta >= theta_start) && (theta <= theta_end);
%     end
% end
% 
% % 角度调整函数
% function theta = adjust_angle(theta)
%     theta = mod(theta + 2*pi, 2*pi);
% end    