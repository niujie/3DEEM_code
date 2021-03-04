function [removeRay_data, data] = RemoveRay(ex_band, em_band, spec_data, scale)
% Remove Rayleigh scattering, change the value to NaN
% Mainly for EX = 500nm beyond 
for i = 1 : 4    % EM, 650-680nm
    for j = 11 : 18    % EX, start from 510nm 
        if spec_data(i, j) - spec_data(i, 10) > 5
            spec_data(i, j) = NaN;    % column 510nm－490nm, if greater than 5, change to NaN
        end
    end
end
% Remove Ramann scattering
% mainly for EX = 490nm, EM = 750nm beyond 
for i = 8 : 11    % EM
    for j = 9 : 16    % EX, start from EX = 480nm
        if spec_data(i, j) - spec_data(i, 8) > 1    % EM >= 720nm，EX >= 490nm
            spec_data(i,j) = NaN;    % if -EX480 >1, chang to NaN
        end
    end
end
 
% 插值去除NaN（散射）
[X, Y] = meshgrid(ex_band, em_band);    % 构建二维矩阵
idxgood = ~isnan(spec_data);    % 非NaN
intp = scatteredInterpolant(X(idxgood), Y(idxgood), spec_data(idxgood));
removeRay_data = intp(X, Y);
 
% 函数绘制光谱图
% 输入：  ex_band = 激发波段
%        em_band = 发射波段
%      spec_band = 波谱数据
%          scale = 插值精度，默认为1，即不插值

if nargin < 4
    % 如果未指定，则默认为1
    scale = 1;
end
d_ex = ex_band(2) - ex_band(1);     % EX波段间隔                    
d_em = em_band(2) - em_band(1);     % EM波段间隔
ex_intp = min(ex_band) : d_ex/scale : max(ex_band);     % EX波段插值
em_intp = min(em_band) : d_em/scale : max(em_band);     % EM波段插值

% 构建二维矩阵：绘图和插值（原始大小）
[EX_ori, EM_ori] = meshgrid(ex_band, em_band);
% 构建二维矩阵：绘图和插值（插值大小）
[EX_intp, EM_intp] = meshgrid(ex_intp, em_intp);
if scale ~=1
    % 二维插值
    data = interp2(EX_ori, EM_ori,removeRay_data,EX_intp, EM_intp,'spline');
else
    % 若scale == 1， 则不改变数据
    data = removeRay_data;
end


% figure                                      % 新建图像
% subplot(1, 2, 1)                            % 子图1
% % imagesc(ex_intp, em_intp, data)             % 绘制二维栅格图
% contour(ex_intp, em_intp, data, 40)         % 绘制20条等高线
% shading flat                                % 去除网格线
% set(gca, 'YDir', 'reverse')                 % 倒转Y轴，使显示与数据文件一致
% xlabel('EX')                                % X轴标签
% ylabel('EM')                                % Y轴标签
% set(gca, 'FontSize', 18)                    % 增大字体
% 
% subplot(1, 2, 2)                            % 子图2
% surf(EX_intp, EM_intp, data+mean(data(:)))  % 三维图，注意图像在Z轴方向抬升
% shading flat                                % 同上
% hold on                                     % 同一张图上同时绘制二维和三维
% contour(ex_intp, em_intp, data, 40)             % 绘制等高线
% xlabel('EX')                                % X轴标签
% ylabel('EM')                                % Y轴标签
% set(gca, 'FontSize', 18)                    % 增大字体
% colorbar                                    % 颜色图例

