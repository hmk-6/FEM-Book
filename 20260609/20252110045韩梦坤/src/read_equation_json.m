function data = read_equation_json(filename)
%READ_EQUATION_JSON 读取作业建议格式的JSON方程组文件。
    if ~exist(filename,'file')
        error('JSON文件不存在：%s',filename);
    end
    data=jsondecode(fileread(filename));
    if isfield(data,'K_FF')
        data.K_FF=double(data.K_FF);
        data.rhs=double(data.rhs(:));
    elseif isfield(data,'K')
        data.K=double(data.K);
        data.R=double(data.R(:));
    else
        error('JSON中未找到K_FF/rhs或K/R字段。');
    end
end
