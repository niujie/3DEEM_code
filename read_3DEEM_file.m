function file_info = read_3DEEM_file(filename)
% �������ڶ�ȡ3Dӫ���������ļ�
% ���룺filename = �ļ���
% �����file_info = ���������ļ������ݵ����ݽṹ
[fid, errmsg] = fopen(filename);
% �ļ������ڣ���������
if fid < 0
    error(errmsg)
end
% �յĽṹ�壬���ڴ洢�ļ��е�����
% ����
% Sample
% File name
% ... ...
% ����Ϣ
file_info = struct();
% ���ж�ȡ�ļ���ֱ�����׾�������֮ǰ
line = [];
sub_field = false;  % �Ƿ�Ϊ�ӽṹ
while ~feof(fid)
    line = fgetl(fid);
    % ɾ��ǰ��ո�
    line = strtrim(line);
    % ��������
    if isempty(line)
        sub_field = false;  % ���б�ʾһ�����ֶνṹ����
        continue
    end
    
    % ����Peaks����
    if strcmpi(line, 'Peaks')
        ApexEX = [];
        ApexEM = [];
        Height = [];
        n = 0;    % ���ݸ���
        % ȥ������ ��No.	ApexEX/EM(nm)	Height(Data)��
        fgetl(fid);
        while true
            n = n + 1;
            line = fgetl(fid);
            if line == -1
                % �ļ�Ϊ�ջ��߶����ļ�β��
                break
            end
            % ֱ���������л���Data Points�ֶ�
            if isempty(line) || strcmpi(line, 'Data Points')
                break
            end
            line = strtrim(line);
            [A, ~, ~] = sscanf(line, '%d %f / %f %f');
            ApexEX(n, 1) = A(2);
            ApexEM(n, 1) = A(3);
            Height(n, 1) = A(4);
        end
        file_info.Peaks = struct();
        file_info.Peaks.ApexEX = ApexEX;
        file_info.Peaks.ApexEM = ApexEM;
        file_info.Peaks.Height = Height;
        continue
    end
    
    % ����Data Points����
    if strcmpi(line, 'Data Points')
        data_pts = [];
        % �����һ�У������ײ���ֵ
        n = 1;
        line = fgetl(fid);
        data = sscanf(line, '%f')';     % ������ת��
        data_pts(n, :) = [NaN, data];
        % �����ļ�β���߿���ֹͣ
        while ~feof(fid) || ~isempty(line)
            line = fgetl(fid);
            if line == -1
                % �ļ�Ϊ�ջ��߶����ļ�β��
                break
            end
            n = n + 1;
            data = sscanf(line, '%f')';   % ��������ת��
            data_pts(n, :) = sscanf(line, '%f')';
        end
        file_info.Data_Points = data_pts;
        continue
    end
    
    k = strfind(line, ':');
    % �������ֶ���Ϣ
    if isempty(k)
        % �ӽṹ
        sub_field = true;
        struct_name = ['file_info.', line];
        % ���ո��滻Ϊ'_'
        struct_name = strrep(struct_name, ' ', '_');
        eval([struct_name, ' = struct();']);
        continue
    end
    % �ֶ���
    field_name = line(1:k-1);
    % ���ո��滻Ϊ'_'
    field_name = strrep(field_name, ' ', '_');
    % �ֶ�����
    info = line(k+1:end);
    % ɾ��ǰ��ո�
    info = strtrim(info);
    % ���ֶ�������Ϣ��ӽ�data_struct
    if ~isempty(info)
        if sub_field
            % ������ӽṹ
            eval([struct_name, '.', field_name, ' = info;']);
        else
            eval(['file_info.', field_name, ' = info;']);
        end
    else
        % ���������
        if sub_field
            % ���ֶ�
            eval([struct_name, '.', field_name, ' = [];']);
        else
            eval(['file_info.', field_name, ' = [];']);
        end
    end
    
end

% �ر��ļ�
fclose(fid);

end
