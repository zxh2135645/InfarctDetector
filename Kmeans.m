clear all;
close all;

base_dir = 'C:\Users\ZhangX1\Documents\MATLAB\masked\';
addpath('C:\Users\ZhangX1\Documents\MATLAB\cviParser\');
sequence_label = {'LGE', 'T1'};
label = sequence_label{2};
anatomys = {'Heart', 'Myocardium', 'excludeContour', 'MyoReference', 'MI'};
anatomy = anatomys{2};


name_glob = glob(cat(2, base_dir, '*\'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'\');
    name = strings{end-1};
    Names{i} = name;
end

RuleOutLabel = NameRuleOutFunc(Names);
Names = Names(RuleOutLabel == 0);

K = 2:5;
infarct_perc_array = zeros(length(Names), length(K));
optK = zeros(length(Names),1);

for j = 1:length(Names)
    exzmple = Names{j};
    disp(exzmple);
    Myo3D = [];
    excludeContour = [];
    MyoPath = cat(2, base_dir ,exzmple, '/', label, '/', anatomy, '/');
    Myo = imageDatastore(MyoPath, 'FileExtensions', '.mat');
    [Myo3D, dicom_idx] = ReadMatFile3D(exzmple, label, anatomy);
    
    if isempty(Myo3D)
        disp('break')
        infarct_perc_array(j,:) = nan(1, length(K));
    else
        Myo3D(Myo3D == 0) = nan;
        clust = zeros(size(Myo3D(:),1),length(K));
        C = cell(1, length(K));
        for k = 1:length(K)
            [clust(:,k), C{k}] = kmeans(Myo3D(:), K(k), 'replicate', 5);
        end
        eva = evalclusters(Myo3D(:),clust,'CalinskiHarabasz');
        optK(j) = eva.OptimalK;
        
        clustered = zeros(size(Myo3D,1), size(Myo3D,2), size(Myo3D,3), length(K));
        for k = 1:length(K)
            clustered(:,:,:,k) = reshape(clust(:,k), size(Myo3D));
        end
        clustered(isnan(clustered)) = 0;
        clustered_sorted = zeros(size(clustered));
        infarct_Kmeans = zeros(size(clustered_sorted));
        
        for k = 1:length(K)
            [B, I] = sort(C{k});
            for i = 1:(k+1)
                clustered_sorted(:,:,:,k) = clustered_sorted(:,:,:,k) + (clustered(:,:,:,k) == I(i))* i; 
            end
            infarct_temp = infarct_Kmeans(:,:,:,k);
            infarct_temp(clustered_sorted(:,:,:,k) == K(k)) = 1;
            infarct_Kmeans(:,:,:,k) = infarct_temp;
        end
        
        excludeCtrPath = cat(2, base_dir ,exzmple, '/', label, '/', anatomys{3}, '/');
        if length(ls(excludeCtrPath)) > 2
            load(cat(2, excludeCtrPath, label, '_excludeContour.mat'))
            [excludeContourMat, exctr_idx] = ReadExContour(excludeContour);
            for i = 1:length(dicom_idx)
                ex_ind = find(dicom_idx(i) == exctr_idx);
                if ~isempty(ex_ind)
                    for k = 1:length(K)
                        infarct_Kmeans(:,:,i,k) = infarct_Kmeans(:,:,i,k) .* ~excludeContourMat(:,:,ex_ind);
                    end
                end
            end
        end
        
        myo_size = sum(Myo3D(:) > 0);

        for k = 1:length(K)
            infarct = ImPostProc(infarct_Kmeans(:,:,:,k));
            infarct_ex_masked = infarct;
            infarct_out = cat(2, base_dir, exzmple, '\', label, '\MI_Kmeans\');
            
            if ~ exist(infarct_out, 'dir')
                mkdir(infarct_out);
            end
            
            % Always overwrite
            infarct_out_path = cat(2, infarct_out, 'MyoInfarct', num2str(k), '.mat');
            save(infarct_out_path, 'infarct_ex_masked');
            
            infarct_size = sum(infarct(:) > 0);
            infarct_perc = infarct_size / myo_size * 100;
            disp(infarct_perc);
            infarct_perc_array(j, k) = infarct_perc;
        end
    end
end

K2 = infarct_perc_array(:,1);
K3 = infarct_perc_array(:,2);
K4 = infarct_perc_array(:,3);
K5 = infarct_perc_array(:,4);
T = table(Names, K2, K3, K4, K5);
out_dir = 'C:/Users/ZhangX1/Documents/MATLAB/masked/';
writetable(T, cat(2, out_dir, 'Kmeans.csv'));

% 
% for i = 1:size(clustered, 3)
%     figure();
%     imagesc(infarct_Kmeans(:,:,i));
%     axis equal
% end
% 
% for i = 1:size(clustered_sorted, 3)
%     figure();
%     imagesc(clustered_sorted(:,:,i));
%     axis equal
% end

