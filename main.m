% ������
clear
clc
% 1. �޸������ļ�·����Ҫ�򿪵��ļ��������ʹ��ȫ·����
% �ļ�·��
% file_path = 'C:\Users\jay_n\Documents\ѧ��\������\interval10\non-toxic\Chlorophyta\P.subcordiformis\';
% filename = [file_path, 'PS2.TXT'];
file_path = 'C:\Users\jay_n\Documents\ѧ��\�ܽ���\����60\14-4-15\';
filename = [file_path, '311.TXT'];

% �ṹ����������ڴ洢�ļ��е�����
% ����
% Sample
% File name
% ... ...
% Data Points
file_info = read_3DEEM_file(filename);

% 2. ��Data Points������ȡ����
%   ex_band = ��������
%   em_band = ���䲨��
% spec_band = ��������
[ex_band, em_band, spec_data] = get_data_points(file_info.Data_Points);

% 3. ������ά����ͼ
scale = 10;     % ��ֵ���ȣ�Ĭ��Ϊ1��������ֵ
draw_3D_spectrum(ex_band, em_band, spec_data, scale)

% 4. ����"��ν"��һά����ͼ
draw_1D_spectrum(spec_data)
