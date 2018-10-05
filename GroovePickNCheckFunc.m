function T = GroovePickNCheckFunc(Names, label, anatomy, base_dir, out_dir, overwrite_label, doublecheck_label)

if nargin == 5
    overwrite_label = 1;
    doublecheck_label = 0;
elseif nargin == 6
    doublecheck_label = 0;
end

addpath('C:\Users\ZhangX1\Documents\MATLAB\cviParser\');
p = 1; % number of landmarks to pick
x = zeros(length(Names), p); % Rows
y = zeros(length(Names), p); % Columns
x_centroid = zeros(length(Names), 1);
y_centroid = zeros(length(Names), 1);
outputFileName = [out_dir, label, '_coords.csv'];

if ~(exist(outputFileName, 'file') && overwrite_label == 0)
    for n = 1:length(Names)
        name = Names{n};
        mat_glob = glob([base_dir, name, '\', label, '\*\*']);
        mat_file = GetFullPath(cat(2, mat_glob{1}, 'VOLUME_IMAGE.mat'));
        load(mat_file, 'volume_image');
        slice_num = ceil(size(volume_image, 3) / 2);
        
        % Load Myocardium mask
        mask_path = cat(2, out_dir, name, '\', label, '\', anatomy, '\', 'masked_heart', num2str(slice_num), '.mat');
        if exist(mask_path, 'file')
            load(mask_path, 'mask_heart');
        else
            mask_path = cat(2, out_dir, name, '\', label, '\', anatomy, '\', 'masked_heart', num2str(slice_num+1), '.mat');
            load(mask_path, 'mask_heart')
        end
        [x_heart, y_heart] = find(mask_heart ~= 0);
        x_centroid(n) = round(mean(x_heart),1);
        y_centroid(n) = round(mean(y_heart),1);
        
        figure();
        imagesc(volume_image(:,:,slice_num))
        truesize([3*size(volume_image,1), 3*size(volume_image,2)]);
        axis off;
        
        for i = 1:p
            [y(n,i), x(n,i)] = ginput(1);
            y(n,i) = round(y(n,i), 1);
            x(n,i) = round(x(n,i), 1);
            h1 = text(y(n,i), x(n,i), '*', 'HorizontalAlignment', 'center', 'Color', [0 0 0], 'FontSize', 16);
            h2 = text(y(n,i), x(n,i), num2str(i), 'HorizontalAlignment', 'center', 'Color', [1 0 0], 'FontSize', 14);
        end
    end
    T = table(Names, x, y, x_centroid, y_centroid);
    writetable(T, outputFileName);
else
    fprintf('The file already exists for %s\n', label);
end

if doublecheck_label == 1
    % double check saved coordinates
    [num,txt,raw] = xlsread(outputFileName);
    for i = 1:length(Names)
        name = Names{i};
        mat_glob = glob([base_dir, name, '\', label, '\*\*']);
        mat_file = GetFullPath(cat(2, mat_glob{1}, 'VOLUME_IMAGE.mat'));
        load(mat_file, 'volume_image');
        slice_num = ceil(size(volume_image, 3) / 2);
        mask_path = cat(2, out_dir, name, '\', label, '\', anatomy, '\', 'masked_heart', num2str(slice_num), '.mat');
        if ~exist(mask_path, 'file')
            slice_num = slice_num + 1;
        end
        
        figure(), imagesc(volume_image(:,:,slice_num)), truesize([3*size(volume_image,1), 3*size(volume_image,2)]);
        axis off;
        hold on; plot(num(i,2), num(i,1),'r*', 'MarkerSize', 12)
    end
end

end