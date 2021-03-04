function [] = draw_3D_spectrum(ex_band, em_band, spec_data, scale)
% �������ƹ���ͼ
% ���룺  ex_band = ��������
%        em_band = ���䲨��
%      spec_band = ��������
%          scale = ��ֵ���ȣ�Ĭ��Ϊ1��������ֵ

if nargin < 4
    % ���δָ������Ĭ��Ϊ1
    scale = 1;
end

d_ex = ex_band(2) - ex_band(1);     % EX���μ��                    
d_em = em_band(2) - em_band(1);     % EM���μ��
ex_intp = min(ex_band) : d_ex/scale : max(ex_band);     % EX���β�ֵ
em_intp = min(em_band) : d_em/scale : max(em_band);     % EM���β�ֵ

% ������ά���󣺻�ͼ�Ͳ�ֵ��ԭʼ��С��
[EX_ori, EM_ori] = meshgrid(ex_band, em_band);
% ������ά���󣺻�ͼ�Ͳ�ֵ��ԭʼ��С��
[EX_intp, EM_intp] = meshgrid(ex_intp, em_intp);
if scale ~=1
    % ��ά��ֵ
    data = interp2(EX_ori, EM_ori, spec_data, EX_intp, EM_intp, 'spline');
else
    % ��scale == 1�� �򲻸ı�����
    data = spec_data;
end


figure                                      % �½�ͼ��
subplot(1, 2, 1)                            % ��ͼ1
% imagesc(ex_intp, em_intp, data)             % ���ƶ�άդ��ͼ
contour(ex_intp, em_intp, data, 20)         % ����20���ȸ���
shading flat                                % ȥ��������
set(gca, 'YDir', 'reverse')                 % ��תY�ᣬʹ��ʾ�������ļ�һ��
xlabel('EX')                                % X���ǩ
ylabel('EM')                                % Y���ǩ
set(gca, 'FontSize', 18)                    % ��������

subplot(1, 2, 2)                            % ��ͼ2
surf(EX_intp, EM_intp, data+mean(data(:)))  % ��άͼ��ע��ͼ����Z�᷽��̧��
shading flat                                % ͬ��
hold on                                     % ͬһ��ͼ��ͬʱ���ƶ�ά����ά
contour(ex_intp, em_intp, data)             % ���Ƶȸ���
xlabel('EX')                                % X���ǩ
ylabel('EM')                                % Y���ǩ
set(gca, 'FontSize', 18)                    % ��������
colorbar                                    % ��ɫͼ��
