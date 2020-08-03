function data_struct = read_3DEEM_file(filename)
% 函数用于读取3D荧光谱数据文件
% 输入：filename = 文件名
% 输出：data_struct = 包含整个文件中数据的数据结构
[fid, errmsg] = fopen(filename);
% 文件不存在，报错，结束
if fid < 0
    error(errmsg)
end
% 空的结构体，用于存储文件中的数据
% 包括
% Sample
% File name
% ... ...
% 等信息
data_struct = struct();
% 按行读取文件，直到光谱矩阵数据之前
line = [];
sub_field = false;  % 是否为子结构
while ~feof(fid)
    line = fgetl(fid);
    % 删除前后空格
    line = strtrim(line);
    % 跳过空行
    if isempty(line)
        sub_field = false;  % 空行表示一个子字段结构结束
        continue
    end
    
    % 处理Peaks数据
    if strcmpi(line, 'Peaks')
        ApexEX = [];
        ApexEM = [];
        Height = [];
        n = 0;    % 数据个数
        % 去掉该行 ‘No.	ApexEX/EM(nm)	Height(Data)’
        fgetl(fid);
        while true
            n = n + 1;
            line = fgetl(fid);
            if line == -1
                % 文件为空或者读到文件尾巴
                break
            end
            % 直到读到空行或者Data Points字段
            if isempty(line) || strcmpi(line, 'Data Points')
                break
            end
            line = strtrim(line);
            [A, ~, ~] = sscanf(line, '%d %f / %f %f');
            ApexEX(n, 1) = A(2);
            ApexEM(n, 1) = A(3);
            Height(n, 1) = A(4);
        end
        data_struct.Peaks = struct();
        data_struct.Peaks.ApexEX = ApexEX;
        data_struct.Peaks.ApexEM = ApexEM;
        data_struct.Peaks.Height = Height;
        continue
    end
    
    % 处理Data Points数据
    if strcmpi(line, 'Data Points')
        data_pts = [];
        % 处理第一行，激发谱波段值
        n = 1;
        line = fgetl(fid);
        data = sscanf(line, '%f')';     % 列向量转置
        data_pts(n, :) = [NaN, data];
        % 读到文件尾或者空行停止
        while ~feof(fid) || ~isempty(line)
            line = fgetl(fid);
            if line == -1
                % 文件为空或者读到文件尾巴
                break
            end
            n = n + 1;
            data = sscanf(line, '%f')';   % 将列向量转置
            data_pts(n, :) = sscanf(line, '%f')';
        end
        data_struct.Data_Points = data_pts;
        continue
    end
    
    k = strfind(line, ':');
    % 跳过非字段信息
    if isempty(k)
        % 子结构
        sub_field = true;
        struct_name = ['data_struct.', line];
        % 将空格替换为'_'
        struct_name = strrep(struct_name, ' ', '_');
        eval([struct_name, ' = struct();']);
        continue
    end
    % 字段名
    field_name = line(1:k-1);
    % 将空格替换为'_'
    field_name = strrep(field_name, ' ', '_');
    % 字段内容
    info = line(k+1:end);
    % 删除前后空格
    info = strtrim(info);
    % 将字段名和信息添加进data_struct
    if ~isempty(info)
        if sub_field
            % 如果是子结构
            eval([struct_name, '.', field_name, ' = info;']);
        else
            eval(['data_struct.', field_name, ' = info;']);
        end
    else
        % 处理空内容
        if sub_field
            % 子字段
            eval([struct_name, '.', field_name, ' = [];']);
        else
            eval(['data_struct.', field_name, ' = [];']);
        end
    end
    
end

% 关闭文件
fclose(fid);

end