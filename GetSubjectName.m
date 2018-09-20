function Names = GetSubjectName(base_dir)

addpath('C:\Users\ZhangX1\Documents\MATLAB\cviParser\');
name_glob = glob(cat(2, base_dir, '*\'));
Names = cell(length(name_glob), 1);

for i = 1:length(name_glob)
    string = split(name_glob{i}, '\');
    Names(i) = string(end-1);
end
end