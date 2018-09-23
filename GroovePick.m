clear all;
close all;

base_dir = 'C:\Users\ZhangX1\Desktop\contour_exporting_Guan\';
out_dir = 'C:\Users\ZhangX1\Documents\MATLAB\masked\';
sequence_label = {'LGE', 'T1'};
anatomy_label = {'Heart', 'Myocardium', 'excludeContour', 'MyoReference'};
anatomy = anatomy_label{1};
overwrite_label = 0;

name_glob = glob(cat(2, base_dir, '/*_*'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'\');
    name = strings{end-1};
    Names{i} = name;
end

for ll = 1:length(sequence_label)
    label = sequence_label{ll};
    GroovePickNCheckFunc(Names, label, anatomy, base_dir, out_dir, overwrite_label);
end

%% double check saved coordinates
doublecheck_label = 1;
label = sequence_label{2};
GroovePickNCheckFunc(Names, label, anaomy, base_dir, out_dir, overwrite_label, doublecheck_label);

%% write coordinates to txt file
% [pathstr, name, ext] = fileparts(im_path);
% outputFileName = [out_dir, label, '_coords.txt'];
% [fid, msg] = fopen(outputFileName, 'w');
% for i = 1:length(x)
%     fprintf(fid, '%.1f\n', x(i));
% end
% for i = 1:length(y)
%     fprintf(fid, '%.1f\n', y(i));
% end
% fclose(fid);





