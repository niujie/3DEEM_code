function [] = draw_1D_spectrum(spec_data)
% ����"��ν"��һά����ͼ
% ���룺spec_data = ����ǿ������

data_em = reshape(spec_data, [1, numel(spec_data)]);    % �ȱ�EM�����EX
data_ex = reshape(spec_data', [1, numel(spec_data)]);   % �ȱ�EX�����EM

figure;                     % �½�ͼ��
subplot(1, 2, 1)            % ��ͼ1
plot(data_ex)               % ����EX����
set(gca, 'FontSize', 18)    % ��������
xlabel('���ݵ�')             % X���ǩ
ylabel('���ӫ��ǿ��')        % Y���ǩ
title('ԭʼӫ���ף�EX��')      % ����

subplot(1, 2, 2)            % ��ͼ2
plot(data_em)               % ����EM����
set(gca, 'FontSize', 18)    % ��������
xlabel('���ݵ�')             % X���ǩ
ylabel('���ӫ��ǿ��')        % Y���ǩ
title('ԭʼӫ���ף�EM��')      % ����
