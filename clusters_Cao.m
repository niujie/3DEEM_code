%% read in non-toxic data
[non_toxic_data, non_toxic_files] = read_data('../3DEEM_DATA/曹分类后原数据/non_toxic/');
non_toxic_data = remove_scatter(non_toxic_data);
% plot_verification(non_toxic_data, non_toxic_files);

%% read in toxic data
[toxic_data, toxic_data_files] = read_data('../3DEEM_DATA/曹分类后原数据/toxic/');
toxic_data = remove_scatter(toxic_data);
% plot_verification(toxic_data, toxic_data_files);

%% merge all data
data = cat(3, non_toxic_data, toxic_data);

% normalize
for i = 1 : size(data, 3)
    temp = data(:, :, i);
    data(:, :, i) = (temp - min(temp(:))) / (max(temp(:)) - min(temp(:)));
end

% reshape to 2D array
data2D = reshape(data, [size(data, 1) * size(data, 2), size(data, 3)]);
% now row for samples
data2D = data2D';


%%
function [data, files] = read_data(path)
files = dir(path);
files(1 : 2) = [];
data = zeros(11, 21, length(files));
for i = 1 : length(files)
    data(:, :, i) = importdata([files(i).folder, '/', files(i).name]);
end
end

%%
function plot_verification(data, files)
[X, Y] = meshgrid(1 : size(data, 2), 1 : size(data, 1));
i = 1;
while i <= size(data, 3)
    ind = mod(i, 25);
    if ind == 1, figure; end
    if ind == 0, ind = 25; end
    subplot(5, 5, ind);
    pcolor(X, Y, reshape(data(:, :, i), [size(data, 1), size(data, 2)]))
    shading flat
    xlabel("x")
    ylabel("y")
    title(files(i).name)
    i = i + 1;
end
end

%% Remove Rayleigh Scattering
function data = remove_scatter(data)
% 主要改EX=500nm以后的数据
for n = 1 : size(data, 3)
    for i = 1 : 4    % em, 650-680nm
        for j = 11 : 18    % ex, 从510nm开始 
            if data(i, j, n) - data(i, 10, n) > 5
                data(i, j, n) = NaN;    % >= 510nm－490nm的列, 若数值大于5, 则输出数据为NaN
            end
        end
    end
end
% 去拉曼散射
% 主要改EX=490nm，EM=750nm以后的数据
for n = 1 : size(data, 3)
    for i = 8 : 11    % EM
        for j = 9 : 16    % EX, 从ex=480nm
            if data(i, j, n) - data(i, 8, n) > 1    % em >= 720nm, ex>=490nm
                data(i, j, n) = NaN;    % 若-ex 480大于1, 输出数据为NaN
            end
        end
    end
end
% Interpolation to remove Nan
[X, Y] = meshgrid(1 : size(data, 2), 1 : size(data, 1));
for i = 1 : size(data, 3)
    temp = data(:, :, i);
    idxgood = ~isnan(temp);  % 输出非0 index
    intp = scatteredInterpolant(X(idxgood), Y(idxgood), temp(idxgood));
    data(:, :, i) = intp(X, Y);
end
end