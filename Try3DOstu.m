clear all;
close all;

addpath('C:/Users/ZhangX1/Documents/MATLAB/cviParser/');
sequence_label = {'LGE', 'T1'};
anatomy_label = {'Heart', 'Myocardium'};
label = char(sequence_label(1));
anatomy = anatomy_label{2};

name_glob = glob(cat(2, 'C:/Users/ZhangX1/Documents/MATLAB/masked/*/'));
mat_glob = glob(cat(2, 'C:/Users/ZhangX1/Documents/MATLAB/masked/CHOI_DAE_SUK/', label, '/', anatomy, '/masked_*.mat'));

% Reordering
idx_array = zeros(length(mat_glob), 1);
for i = 1:length(mat_glob)
    B = regexp(mat_glob(i),'\d*','Match');
    
    for ii= 1:length(B)
        if ~isempty(B{ii})
            Num(ii,1)=str2double(B{ii}(end));
        else
            Num(ii,1)=NaN;
        end
    end
    fprintf('%d\n', Num)
    idx_array(i) = Num;
end

% Display masks
gap_map = zeros(length(mat_glob), 1);
min_idx = min(idx_array);
gap = min_idx - 1;

for i = 1:length(mat_glob)
    idx = find(idx_array == i + gap, 1);
    if isempty(idx)
        gap = gap + 1;
    end
    gap_map(i) = gap;
end

idx = find(idx_array == 4 + gap_map(4));
readout = load(mat_glob{idx});
I = readout.mask_myocardium;
I(I == 0) = NaN;

%% 2D Ostu's
uniqLevels = unique(I(:));
disp(['Number of unique levels = ' int2str(length(uniqLevels))]);
Nvals = [1 2 3 4 5 6 7 8 9 10];
for i = 1:length(Nvals)
    [thresh, metric] = multithresh(I, Nvals(i) );
    disp(['N = ' int2str(Nvals(i)) '  |  metric = ' num2str(metric)]);
end

[thresh, metric] = multithresh(I, 2);
I(isnan(I)) = 0;
seg_Neq2 = imquantize(I, thresh);

binary_mask = I > 0;
seg_Neq2 = binary_mask .* seg_Neq2;

figure();
imagesc(seg_Neq2)
axis equal

%% Try 3D

for i = 1:length(mat_glob)
    idx = find(idx_array == i + gap_map(i));
    readout = load(mat_glob{idx});
    if strcmp(anatomy, 'Myocardium')
        I = readout.mask_myocardium;
    end
    
    volume_img(:,:,i) = I;
end
volume_img(volume_img == 0) = NaN;
uniqLevels = unique(volume_img(:));
disp(['Number of unique levels = ' int2str(length(uniqLevels))]);
Nvals = [1 2 3 4 5 6 7 8 9 10];
for i = 1:length(Nvals)
    [thresh, metric] = multithresh(volume_img, Nvals(i) );
    disp(['N = ' int2str(Nvals(i)) '  |  metric = ' num2str(metric)]);
end

[thresh, metric] = multithresh(volume_img, 3);
volume_img(isnan(volume_img)) = 0;
seg_Neq2 = imquantize(volume_img, thresh);

binary_mask = volume_img > 0;
seg_Neq2 = binary_mask .* seg_Neq2;

n = ceil(sqrt(length(mat_glob)));
figure();
for i = 1:length(mat_glob)
    subplot(n, n, i)
    imagesc(seg_Neq2(:,:,i))
    axis equal
end

%% T1 Mapping

label = char(sequence_label(2));
anatomy = anatomy_label{2};

name_glob = glob(cat(2, 'C:/Users/ZhangX1/Documents/MATLAB/masked/*/'));
mat_glob = glob(cat(2, 'C:/Users/ZhangX1/Documents/MATLAB/masked/LEE_KWAN_JOON/', label, '/', anatomy, '/masked_*.mat'));

% Reordering
idx_array = zeros(length(mat_glob), 1);
for i = 1:length(mat_glob)
    B = regexp(mat_glob(i),'\d*','Match');
    
    for ii= 1:length(B)
        if ~isempty(B{ii})
            Num(ii,1)=str2double(B{ii}(end));
        else
            Num(ii,1)=NaN;
        end
    end
    fprintf('%d\n', Num)
    idx_array(i) = Num;
end

% Display masks
gap_map = zeros(length(mat_glob), 1);
min_idx = min(idx_array);
gap = min_idx - 1;

for i = 1:length(mat_glob)
    idx = find(idx_array == i + gap, 1);
    if isempty(idx)
        gap = gap + 1;
    end
    gap_map(i) = gap;
end

idx = find(idx_array == 4 + gap_map(4));
readout = load(mat_glob{idx});
I = readout.mask_myocardium;
I(I == 0) = NaN;

%% 2D Ostu's to T1 map
uniqLevels = unique(I(:));
disp(['Number of unique levels = ' int2str(length(uniqLevels))]);
Nvals = [1 2 3 4 5 6 7 8 9 10];
for i = 1:length(Nvals)
    [thresh, metric] = multithresh(I, Nvals(i) );
    disp(['N = ' int2str(Nvals(i)) '  |  metric = ' num2str(metric)]);
end

[thresh, metric] = multithresh(I, 4);
I(isnan(I)) = 0;
seg_Neq2 = imquantize(I, thresh);

binary_mask = I > 0;
seg_Neq2 = binary_mask .* seg_Neq2;

figure();
imagesc(seg_Neq2)
axis equal

%% Try 3D T1 Mapping
clear volume_img
for i = 1:length(mat_glob)
    idx = find(idx_array == i + gap_map(i));
    readout = load(mat_glob{idx});
    if strcmp(anatomy, 'Myocardium')
        I = readout.mask_myocardium;
    end
    
    volume_img(:,:,i) = I;
end
volume_img(volume_img == 0) = NaN;
uniqLevels = unique(volume_img(:));
disp(['Number of unique levels = ' int2str(length(uniqLevels))]);
Nvals = [1 2 3 4 5 6 7 8 9 10];
for i = 1:length(Nvals)
    [thresh, metric] = multithresh(volume_img, Nvals(i) );
    disp(['N = ' int2str(Nvals(i)) '  |  metric = ' num2str(metric)]);
end

[thresh, metric] = multithresh(volume_img, 1);
volume_img(isnan(volume_img)) = 0;
seg_Neq2 = imquantize(volume_img, thresh);

binary_mask = volume_img > 0;
seg_Neq2 = binary_mask .* seg_Neq2;

n = ceil(sqrt(length(mat_glob)));
figure();
% imagesc(seg_Neq2(:,:,4));
% axis equal

for i = 1:length(mat_glob)
    subplot(n, n, i)
    imagesc(seg_Neq2(:,:,i))
    axis equal
end

%% Read and Display 08/23/2018
name_glob = glob(cat(2, 'C:/Users/ZhangX1/Documents/MATLAB/masked/*/'));
name = 'KIM_BONG_KI';
label = 'T1';
anatomy = 'Myocardium';
[mask, dicom_idx] = ReadMatFile3D(name, label, anatomy);
load(cat(2, 'C:/Users/ZhangX1/Documents/MATLAB/masked/', name, '/', label, '/', anatomy_label{3}, '/T1_excludeContour.mat'));
fname = fieldnames(excludeContour);
mask(mask == 0) = NaN;
Nvals = [1 2 3 4];
inf_size = zeros(length(Nvals), 1);

for j = 1:length(Nvals)
    figure();
    [thresh, metric] = multithresh(mask, Nvals(j) );
    mask(isnan(mask)) = 0;
    seg_Neq = imquantize(mask, thresh);
    binary_mask = mask > 0;
    seg_Neq = binary_mask .* seg_Neq;
    n = ceil(sqrt(size(seg_Neq, 3)));
%     for i = 1:size(seg_Neq, 3)
%         subplot(n,n,i)
%         imagesc(seg_Neq(:,:,i))
%         axis equal
%     end
    if ~ isempty(fname)
        exctr_idx = excludeContour.(fname{1}){2};
        for k = 1: length(dicom_idx)
            if ~isempty(excludeContour.(fname{1}))
                ind = find(dicom_idx(k) == exctr_idx);
                if ~isempty(ind)
                    counterMask = ~ excludeContour.(fname{1}){1}(:,:,ind);
                    seg_Neq(:,:,k) = counterMask .* seg_Neq(:,:,k);
                end
            end
        end
    end
    for i = 1:size(seg_Neq, 3)
        subplot(n,n,i)
        imagesc(seg_Neq(:,:,i))
        axis equal
    end
    inf_size(j) = sum(seg_Neq(:) == Nvals(j)+1);
    fprintf('infarct size: %d #ofThresh: %d \n', sum(seg_Neq(:) == Nvals(j)+1), j);
    
end

lv_size(i) = sum(binary_mask(:) == 1);
infarct_size(i, :) = {name, inf_size};

