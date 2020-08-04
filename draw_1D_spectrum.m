function [] = draw_1D_spectrum(spec_data)
% 绘制"所谓"的一维光谱图
% 输入：spec_data = 光谱强度数据

data_em = reshape(spec_data, [1, numel(spec_data)]);    % 先变EM，后变EX
data_ex = reshape(spec_data', [1, numel(spec_data)]);   % 先变EX，后变EM

figure;                     % 新建图像
subplot(1, 2, 1)            % 子图1
plot(data_ex)               % 绘制EX波谱
set(gca, 'FontSize', 18)    % 增大字体
xlabel('数据点')             % X轴标签
ylabel('相对荧光强度')        % Y轴标签
title('原始荧光谱（EX）')      % 标题

subplot(1, 2, 2)            % 子图2
plot(data_em)               % 绘制EM波谱
set(gca, 'FontSize', 18)    % 增大字体
xlabel('数据点')             % X轴标签
ylabel('相对荧光强度')        % Y轴标签
title('原始荧光谱（EM）')      % 标题
