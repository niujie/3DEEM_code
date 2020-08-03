% 主程序
clear
clc
% 1. 修改以下文件路径及要打开的文件名（最好使用全路径）
% 文件路径
% file_path = 'C:\Users\jay_n\Documents\学生\苟偲钰\interval10\non-toxic\Chlorophyta\P.subcordiformis\';
% filename = [file_path, 'PS2.TXT'];
file_path = 'C:\Users\jay_n\Documents\学生\曹洁茹\光照60\14-4-15\';
filename = [file_path, '311.TXT'];

data = read_3DEEM_file(filename);