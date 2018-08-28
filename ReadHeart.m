% Read .mat files
clear all;
close all;


addpath('C:/Users/ZhangX1/Documents/MATLAB/cviParser/');
sequence_label = {'LGE', 'T1'};
anatomy_label = {'Heart', 'Myocardium'};
label = char(sequence_label(1));
anatomy = anatomy_label{2};

name_glob = glob(cat(2, 'C:/Users/ZhangX1/Documents/MATLAB/masked/*/'));

mat_glob = glob(cat(2, 'C:/Users/ZhangX1/Documents/MATLAB/masked/LEE_KWAN_JOON/', label, '/', anatomy, '/masked_*.mat'));

%% Reordering
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

%% Display masks
gap_map = zeros(length(mat_glob), 1);
min_idx = min(idx_array);
gap = min_idx - 1;

for i = 1:length(mat_glob)
    idx = find(idx_array == i + gap);
    if isempty(idx)
        gap = gap + 1;
    end
    gap_map(i) = gap;
end

for i = 1:length(mat_glob)
    figure();
    idx = find(idx_array == i + gap_map(i));
    
    readout = load(char(mat_glob(idx)));
    if strcmp(anatomy, 'Heart')
        imagesc(readout.mask_heart);
    elseif strcmp(anatomy, 'Myocardium')
        imagesc(readout.mask_myocardium)
    end
    axis equal;
    
end

