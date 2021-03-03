%% non-toxic data
[non_toxic_ex_bands, non_toxic_em_bands, ...
    non_toxic_data, non_toxic_files] = read_data('../3DEEM_DATA/苟原数据/non_toxic/');
% plot_verification(non_toxic_ex_bands, non_toxic_em_bands, non_toxic_data, non_toxic_files);

n_non_toxic_samples = size(non_toxic_data, 3);  % number of non-toxic data samples

%% toxic data
[toxic_ex_bands, toxic_em_bands, ...
    toxic_data, toxic_files] = read_data('../3DEEM_DATA/苟原数据/toxic/');
% plot_verification(toxic_ex_bands, toxic_em_bands, toxic_data, toxic_files);

n_toxic_samples = size(toxic_data, 3);  % number of toxic data samples

%% merge all data
data = cat(3, non_toxic_data, toxic_data);
ex_bands = [non_toxic_ex_bands; toxic_ex_bands];
em_bands = [non_toxic_em_bands; toxic_em_bands];
n_samples = size(data, 3);  % number of data samples
assert(n_samples == n_non_toxic_samples + n_toxic_samples)  % assert if the number of samples is correct

%% normalize
for i = 1 : n_samples
    temp = data(:, :, i);
    data(:, :, i) = (temp - min(temp(:))) / (max(temp(:)) - min(temp(:)));
end

%% reshape to 2D array
data2D = reshape(data, [size(data, 1) * size(data, 2), n_samples]);
% now row for samples
data2D = data2D';

%% create the one-hot target
% non-toxic = 0
% targets = zeros(size(data2D, 1), 2);
% % toxic = 1
% targets(size(non_toxic_data, 3) + 1 : end, 2) = 1;
targets = zeros(n_samples, 1);
targets(n_non_toxic_samples + 1 : end) = 1;

%% random sample
ind = randsample(1 : n_samples, n_samples);
data = data(:, :, ind);
data2D = data2D(ind, :);
targets = targets(ind, :);


%% neural network
% create 4D array for image network training
dataset = zeros(size(data,1), size(data,2), 1, n_samples);
for i = 1 : n_samples, dataset(:,:,1,i) = data(:,:,i); end

% create categorical for network training
Y = categorical(targets);

train_ind = 1 : round(0.6 * n_samples);   % training data indices
valid_ind = train_ind(end) + 1 : n_samples; % validation data indices
trainX = dataset(:,:,:,train_ind);
trainY = Y(train_ind);
validX = dataset(:,:,:,valid_ind);
validY = Y(valid_ind);

%% train network
layers = [
    imageInputLayer([11 21 1], "Name", "imageinput")
    batchNormalizationLayer("Name", "batchnorm")
    fullyConnectedLayer(20, "Name", "fc_1")
    dropoutLayer(0.5, "Name", "dropout")
    fullyConnectedLayer(2, "Name", "fc_2")
    softmaxLayer("Name", "softmax")
    classificationLayer("Name", "classoutput")];

options = trainingOptions('adam', ...
    'ExecutionEnvironment', 'cpu', ...
    'MaxEpochs', 100, ...
    'MiniBatchSize', 10, ...
    'ValidationData', {validX, validY}, ...
    'GradientThreshold', 2, ...
    'Shuffle', 'every-epoch', ...
    'Verbose', true, ...
    'Plots', 'training-progress');

net = trainNetwork(trainX, trainY, layers, options);
predictedY = classify(net, dataset);
cm = confusionchart(Y, predictedY);


%%
function [ex_bands, em_bands, data, files] = read_data(path)
files = dir(path);
files(1 : 2) = [];
data = zeros(11, 21, length(files));
ex_bands = zeros(length(files), 21);
em_bands = zeros(length(files), 11);
for i = 1 : length(files)
    file_info = read_3DEEM_file([files(i).folder, '/', files(i).name]);
    [ex_band, em_band, spec_data] = get_data_points(file_info.Data_Points);
    [removeRay_data, ~] = removeRay(ex_band, em_band, spec_data);
    data(:, :, i) = removeRay_data;
    ex_bands(i, :) = ex_band;
    em_bands(i, :) = em_band;
end
end

%%
function plot_verification(ex_bands, em_bands, data, files)
i = 1;
while i <= size(data, 3)
    ind = mod(i, 25);
    if ind == 1, figure; end
    if ind == 0, ind = 25; end
    subplot(5, 5, ind);
    [X, Y] = meshgrid(ex_bands(i, :), em_bands(i, :));
    pcolor(X, Y, reshape(data(:, :, i), [size(data, 1), size(data, 2)]))
    shading flat
    xlabel("x")
    ylabel("y")
    title(files(i).name)
    i = i + 1;
end
end