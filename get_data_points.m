function [ex_band, em_band, spec_data] = get_data_points(data_pts)
% ��Data Points������ȡ����
% ���룺 data_pts = ���������ļ��е�Data Points, file_info.Data_Points
% �����  ex_band = ��������
%        em_band = ���䲨��
%      spec_band = ��������

% ��ȡ��������
ex_band = data_pts(1, :);
% ȥ��NaN
ex_band(isnan(ex_band)) = [];

% ��ȡ���䲨��
em_band = data_pts(2:end, 1);
em_band = em_band';     % ת��

% ��ȡ��������
spec_data = data_pts(2:end, 2:end);
