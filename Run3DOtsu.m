clear all;
close all;

sequence_label = {'LGE', 'T1'};
anatomy_label = {'Heart', 'Myocardium', 'excludeContour'};
label = char(sequence_label(2));
anatomy = anatomy_label{2};

name_glob = glob(cat(2, 'C:/Users/ZhangX1/Documents/MATLAB/masked/*/'));

infarct_size = cell(length(name_glob), 2);
lv_size = zeros(length(name_glob), 1);
for i = 1:length(name_glob)
   strings = strsplit(name_glob{i},'\'); 
   name = strings{end-1};
   disp(name)
   if length(ls(cat(2, 'C:/Users/ZhangX1/Documents/MATLAB/masked/', name, '/', label, '/', anatomy))) > 2
        [mask, dicom_idx] = ReadMatFile3D(name, label, anatomy);
        load(cat(2, 'C:/Users/ZhangX1/Documents/MATLAB/masked/', name, '/', label, '/', anatomy_label{3}, '/', label, '_excludeContour.mat'));
        fname = fieldnames(excludeContour);
        
        Nvals = [1 2 3 4];
        inf_size = zeros(length(Nvals), 1);
        
        for j = 1:length(Nvals)
            mask(mask == 0) = NaN;
            [thresh, metric] = multithresh(mask, Nvals(j) );
            mask(isnan(mask)) = 0;
            seg_Neq = imquantize(mask, thresh);
            
            binary_mask = mask > 0;
            seg_Neq = binary_mask .* seg_Neq;
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
            inf_size(j) = sum(seg_Neq(:) == Nvals(j)+1);
            fprintf('infarct size: %d #ofThresh: %d \n', sum(seg_Neq(:) == Nvals(j)+1), j);
            
        end
        lv_size(i) = sum(binary_mask(:) == 1);
        infarct_size(i, :) = {name, inf_size};
        
   else
       disp('Folder is empty!');
       infarct_size(i, :) = {name, []};
   end
end

%% Convert to matrix
clear infarct_matrix lv
gap = 0;
for i = 1: length(name_glob)
    if ~isempty(infarct_size{i, 2})
        infarct_matrix(:, i-gap) = infarct_size{i, 2};
        lv(i-gap) = lv_size(i);
    else
        gap = gap + 1;
    end
end
infarct_matrix = infarct_matrix* 1.1875^2 / (10^2);
true_lv = lv * 1.1875^2 / (10^2);
infarct_mean = mean(infarct_matrix, 2);
infarct_std = std(infarct_matrix, 0, 2);

inf_fraction = infarct_matrix ./ true_lv;
inf_fraction_perc = inf_fraction * 100;

%% Error bar
figure();
errorbar(Nvals,infarct_mean,infarct_std,'-s','MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red')
xlabel('Number of threshold'); ylabel('Infarct Area (cm^2)')
title('3D Otsu on LGE');
grid on;

save('T1OtsuMat.mat', 'inf_fraction');

% %% T1 map
% sequence_label = {'LGE', 'T1'};
% anatomy_label = {'Heart', 'Myocardium'};
% label = char(sequence_label(2));
% anatomy = anatomy_label{2};
% 
% name_glob = glob(cat(2, 'C:/Users/ZhangX1/Documents/MATLAB/masked/*/'));
% 
% infarct_size = cell(length(name_glob), 2);
% lv_size = zeros(length(name_glob), 1);
% for i = 1:length(name_glob)
%    strings = strsplit(name_glob{i},'\'); 
%    name = strings{end-1};
%    disp(name)
%    if length(ls(cat(2, 'C:/Users/ZhangX1/Documents/MATLAB/masked/', name, '/', label, '/', anatomy))) > 2
%         mask = ReadMatFile3D(name, label, anatomy);
%         mask(mask == 0) = NaN;
%         Nvals = [1 2 3 4];
%         inf_size = zeros(length(Nvals), 1);
% 
%         for j = 1:length(Nvals)
%             [thresh, metric] = multithresh(mask, Nvals(j) );
%             mask(isnan(mask)) = 0;
%             seg_Neq = imquantize(mask, thresh);
%             
%             binary_mask = mask > 0;
%             seg_Neq = binary_mask .* seg_Neq;
%             inf_size(j) = sum(seg_Neq(:) == Nvals(j)+1);
%             fprintf('infarct size: %d #ofThresh: %d \n', sum(seg_Neq(:) == Nvals(j)+1), j);
%         end
%         lv_size(i) = sum(binary_mask(:) == 1);
%         infarct_size(i, :) = {name, inf_size};
%         
%    else
%        disp('Folder is empty!');
%        infarct_size(i, :) = {name, []};
%    end
% end
% 
% % Convert to matrix
% clear infarct_matrix lv
% gap = 0;
% for i = 1: length(name_glob)
%     if ~isempty(infarct_size{i, 2})
%         infarct_matrix(:, i-gap) = infarct_size{i, 2};
%         lv(i-gap) = lv_size(i);
%     else
%         gap = gap + 1;
%     end
% end
% infarct_matrix = infarct_matrix* 1.484375^2 / (10^2);
% true_lv = lv * 1.484375^2 / (10^2);
% infarct_mean = mean(infarct_matrix, 2);
% infarct_std = std(infarct_matrix, 0, 2);
% 
% inf_fraction = infarct_matrix ./ true_lv;
% inf_fraction_perc = inf_fraction * 100;
% 
% % Error bar
% figure();
% errorbar(Nvals,infarct_mean,infarct_std,'-s','MarkerSize',10,...
%     'MarkerEdgeColor','red','MarkerFaceColor','red')
% xlabel('Number of threshold'); ylabel('Infarct Area (cm^2)')
% title('3D Otsu on T1 Map');
% grid on;
% save('T1OtsuMat.mat', 'infarct_matrix');


%%
% clear volume_img
% for i = 1:length(mat_glob)
%     idx = find(idx_array == i + gap_map(i));
%     readout = load(mat_glob{idx});
%     if strcmp(anatomy, 'Myocardium')
%         I = readout.mask_myocardium;
%     end
%     
%     volume_img(:,:,i) = I;
% end

%%
sequence_label = {'LGE', 'T1'};
anatomy_label = {'Heart', 'Myocardium', 'excludeContour'};
label = char(sequence_label(2));
anatomy = anatomy_label{2};

name_glob = glob(cat(2, 'C:/Users/ZhangX1/Documents/MATLAB/masked/*/'));

infarct_size = cell(length(name_glob), 2);
lv_size = zeros(length(name_glob), 1);
for i = 1:length(name_glob)
   strings = strsplit(name_glob{i},'\'); 
   name = strings{end-1};
   disp(name)
   if length(ls(cat(2, 'C:/Users/ZhangX1/Documents/MATLAB/masked/', name, '/', label, '/', anatomy))) > 2
        [mask, dicom_idx] = ReadMatFile3D(name, label, anatomy);
        load(cat(2, 'C:/Users/ZhangX1/Documents/MATLAB/masked/', name, '/', label, '/', anatomy_label{3}, '/', label, '_excludeContour.mat'));
        fname = fieldnames(excludeContour);
        
        Nvals = [1 2 3 4];
        inf_size = zeros(length(Nvals), 1);
        
        for j = 1:length(Nvals)
            if ~ isempty(fname)
                exctr_idx = excludeContour.(fname{1}){2};
                for k = 1: length(dicom_idx)
                    if ~isempty(excludeContour.(fname{1}))
                        ind = find(dicom_idx(k) == exctr_idx);
                        if ~isempty(ind)
                            counterMask = ~ excludeContour.(fname{1}){1}(:,:,ind);
                            mask(:,:,k) = counterMask .* mask(:,:,k);
                        end
                    end
                end
            end
            mask(mask == 0) = NaN;
            [thresh, metric] = multithresh(mask, Nvals(j) );
            mask(isnan(mask)) = 0;
            seg_Neq = imquantize(mask, thresh);
            
            binary_mask = mask > 0;
            seg_Neq = binary_mask .* seg_Neq;

            inf_size(j) = sum(seg_Neq(:) == Nvals(j)+1);
            fprintf('infarct size: %d #ofThresh: %d \n', sum(seg_Neq(:) == Nvals(j)+1), j);
            
        end
        lv_size(i) = sum(binary_mask(:) == 1);
        infarct_size(i, :) = {name, inf_size};
        
   else
       disp('Folder is empty!');
       infarct_size(i, :) = {name, []};
   end
end

%% Convert to matrix
clear infarct_matrix lv
gap = 0;
for i = 1: length(name_glob)
    if ~isempty(infarct_size{i, 2})
        infarct_matrix(:, i-gap) = infarct_size{i, 2};
        lv(i-gap) = lv_size(i);
    else
        gap = gap + 1;
    end
end
infarct_matrix = infarct_matrix* 1.1875^2 / (10^2);
true_lv = lv * 1.1875^2 / (10^2);
infarct_mean = mean(infarct_matrix, 2);
infarct_std = std(infarct_matrix, 0, 2);

inf_fraction = infarct_matrix ./ true_lv;
inf_fraction_perc = inf_fraction * 100;