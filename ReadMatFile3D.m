function [mask, dicom_idx] = ReadMatFile3D(name, modality, anatomy, base_dir)

if nargin == 3
    base_dir = 'C:/Users/ZhangX1/Documents/MATLAB/masked/';
end

% Too include glob function
addpath('C:/Users/ZhangX1/Documents/MATLAB/cviParser/');

mat_glob = glob(cat(2, base_dir, name, '/', modality, '/', anatomy, '/masked_*.mat'));

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
    % fprintf('%d\n', Num)
    idx_array(i) = Num;
end

% Generate gap map
gap_map = zeros(length(mat_glob), 1);
min_idx = min(idx_array);
gap = min_idx - 1;

for i = 1:length(mat_glob)
    idx = find(idx_array == i + gap, 1);
    while isempty(idx)
        gap = gap + 1;
        idx = find(idx_array == i + gap, 1);
    end
    gap_map(i) = gap;
end

for i = 1:length(mat_glob)

    idx = find(idx_array == i + gap_map(i));
    
    readout = load(char(mat_glob(idx)));
    if strcmp(anatomy, 'Heart')
        mask(:,:,i) = readout.mask_heart;
    elseif strcmp(anatomy, 'Myocardium')
        mask(:,:,i) = readout.mask_myocardium;
    end

end

dicom_idx = sort(idx_array);

end