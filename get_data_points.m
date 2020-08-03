function [ex_band, em_band, spec_data] = get_data_points(data_pts)
% 从Data Points数组提取数据
% 输入： data_pts = 光谱数据文件中的Data Points, file_info.Data_Points
% 输出：  ex_band = 激发波段
%        em_band = 发射波段
%      spec_band = 波谱数据

% 提取激发波段
ex_band = data_pts(1, :);
% 去除NaN
ex_band(isnan(ex_band)) = [];

% 提取发射波段
em_band = data_pts(2:end, 1);
em_band = em_band';     % 转置

% 提取波谱数据
spec_data = data_pts(2:end, 2:end);