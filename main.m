% 主程序
clear
clc
% 1. 修改以下文件路径及要打开的文件名（最好使用全路径）
% 文件路径
% file_path = 'C:\Users\jay_n\Documents\学生\苟偲钰\interval10\non-toxic\Chlorophyta\P.subcordiformis\';
% filename = [file_path, 'PS2.TXT'];
file_path = 'C:\Users\jay_n\Documents\学生\曹洁茹\光照60\14-4-15\';
filename = [file_path, '311.TXT'];

% 结构体变量，用于存储文件中的数据
% 包括
% Sample
% File name
% ... ...
% Data Points
file_info = read_3DEEM_file(filename);

% 提取激发波段
ex_band = file_info.Data_Points(1, :);
% 去除NaN
ex_band(isnan(ex_band)) = [];

% 提取发射波段
em_band = file_info.Data_Points(2:end, 1);
em_band = em_band';     % 转置

% 提取波谱数据
spec_data = file_info.Data_Points(2:end, 2:end);