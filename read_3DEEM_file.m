function data_struct = read_3DEEM_file(filename)
% �������ڶ�ȡ3Dӫ���������ļ�
% ���룺filename = �ļ���
% �����data_struct = ���������ļ������ݵ����ݽṹ
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
data_struct = struct();
% ���ж�ȡ�ļ���ֱ�����׾�������֮ǰ
line = [];
sub_field = false;  % �Ƿ�Ϊ�ӽṹ
while ~strcmpi(line, 'Data Points') && ~strcmpi(line, 'Peaks')
    line = fgetl(fid);
    % ɾ��ǰ��ո�
    line = strtrim(line);
    % ��������
    if isempty(line)
        sub_field = false;  % ���б�ʾһ�����ֶνṹ����
        continue
    end
    k = strfind(line, ':');
    % �������ֶ���Ϣ
    if isempty(k)
        % �ӽṹ
        sub_field = true;
        struct_name = ['data_struct.', line];
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
            eval(['data_struct.', field_name, ' = info;']);
        end
    else
        % ���������
        if sub_field
            % ���ֶ�
            eval([struct_name, '.', field_name, ' = [];']);
        else
            eval(['data_struct.', field_name, ' = [];']);
        end
    end
end

% ����Peaks����
ApexEX = [];
ApexEM = [];
Height = [];
n = 0;    % ���ݸ���
if strcmp(line, 'Peaks')
    % ȥ������ ��No.	ApexEX/EM(nm)	Height(Data)��
    fgetl(fid);
    while true
        n = n + 1;
        line = fgetl(fid);
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
end
data_struct.Peaks = struct();
data_struct.Peaks.ApexEX = ApexEX;
data_struct.Peaks.ApexEM = ApexEM;
data_struct.Peaks.Height = Height;

% �ر��ļ�
fclose(fid);

end