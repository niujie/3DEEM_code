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

% ��ȡ��������
ex_band = file_info.Data_Points(1, :);
% ȥ��NaN
ex_band(isnan(ex_band)) = [];

% ��ȡ���䲨��
em_band = file_info.Data_Points(2:end, 1);
em_band = em_band';     % ת��

% ��ȡ��������
spec_data = file_info.Data_Points(2:end, 2:end);