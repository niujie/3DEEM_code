function [removeRay_data, data] = RemoveRay(ex_band, em_band, spec_data, scale)
% Remove Rayleigh scattering, change the value to NaN
% Mainly for EX = 500nm beyond 
for i = 1 : 4    % EM, 650-680nm
    for j = 11 : 18    % EX, start from 510nm 
        if spec_data(i, j) - spec_data(i, 10) > 5
            spec_data(i, j) = NaN;    % column 510nm��490nm, if greater than 5, change to NaN
        end
    end
end
% Remove Ramann scattering
% mainly for EX = 490nm, EM = 750nm beyond 
for i = 8 : 11    % EM
    for j = 9 : 16    % EX, start from EX = 480nm
        if spec_data(i, j) - spec_data(i, 8) > 1    % EM >= 720nm��EX >= 490nm
            spec_data(i,j) = NaN;    % if -EX480 >1, chang to NaN
        end
    end
end
 
% ��ֵȥ��NaN��ɢ�䣩
[X, Y] = meshgrid(ex_band, em_band);    % ������ά����
idxgood = ~isnan(spec_data);    % ��NaN
intp = scatteredInterpolant(X(idxgood), Y(idxgood), spec_data(idxgood));
removeRay_data = intp(X, Y);
 
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
% ������ά���󣺻�ͼ�Ͳ�ֵ����ֵ��С��
[EX_intp, EM_intp] = meshgrid(ex_intp, em_intp);
if scale ~=1
    % ��ά��ֵ
    data = interp2(EX_ori, EM_ori,removeRay_data,EX_intp, EM_intp,'spline');
else
    % ��scale == 1�� �򲻸ı�����
    data = removeRay_data;
end


% figure                                      % �½�ͼ��
% subplot(1, 2, 1)                            % ��ͼ1
% % imagesc(ex_intp, em_intp, data)             % ���ƶ�άդ��ͼ
% contour(ex_intp, em_intp, data, 40)         % ����20���ȸ���
% shading flat                                % ȥ��������
% set(gca, 'YDir', 'reverse')                 % ��תY�ᣬʹ��ʾ�������ļ�һ��
% xlabel('EX')                                % X���ǩ
% ylabel('EM')                                % Y���ǩ
% set(gca, 'FontSize', 18)                    % ��������
% 
% subplot(1, 2, 2)                            % ��ͼ2
% surf(EX_intp, EM_intp, data+mean(data(:)))  % ��άͼ��ע��ͼ����Z�᷽��̧��
% shading flat                                % ͬ��
% hold on                                     % ͬһ��ͼ��ͬʱ���ƶ�ά����ά
% contour(ex_intp, em_intp, data, 40)             % ���Ƶȸ���
% xlabel('EX')                                % X���ǩ
% ylabel('EM')                                % Y���ǩ
% set(gca, 'FontSize', 18)                    % ��������
% colorbar                                    % ��ɫͼ��

