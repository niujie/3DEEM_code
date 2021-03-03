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

% 2. 从Data Points数组提取数据
%   ex_band = 激发波段
%   em_band = 发射波段
% spec_band = 波谱数据
[ex_band, em_band, spec_data] = get_data_points(file_info.Data_Points);

% 3. 绘制三维光谱图
scale = 10;     % 插值精度，默认为1，即不插值
draw_3D_spectrum(ex_band, em_band, spec_data, scale)

% 4. 绘制"所谓"的一维光谱图
draw_1D_spectrum(spec_data)
