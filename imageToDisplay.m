clear all;
close all;
%%
%name = 'OH_SANG_MOO';
%name = 'CHOI_DAE_SUK';
%name = 'HAN_BONG_SANG';
% name = 'HWANG_IN_YONG';
% name = 'KIM_KEUM_HUN';
base_dir = 'C:\Users\ZhangX1\Documents\MATLAB\masked\';
addpath('C:\Users\ZhangX1\Documents\MATLAB\cviParser\');
sequence_label = {'LGE', 'T1'};
anatomys = {'Heart', 'Myocardium', 'excludeContour', 'MyoReference', 'MI'};

dicom_glob = glob(cat(2, 'C:/Users/ZhangX1/Desktop/contour_exporting_Guan/', name, '/', sequence_label{1}, '/*/*'));
dicom = dicom_glob{1};
[volume_image, slice_data, image_meta_data] = dicom23D(dicom);
[mask, dicom_idx] = ReadMatFile3D(name, sequence_label{1}, anatomys{1});

% for i = 1:size(volume_image, 3)
%     figure();
%     imagesc(volume_image(:,:,i))
%     axis equal
%     colormap gray
% end
% 
% for i = 1:length(dicom_idx)
%     figure();
%     imagesc(mask(:,:,i))
%     axis equal
%     colormap gray
% end

load(cat(2, base_dir, name, '/', sequence_label{1}, '/MI/MyoInfarct.mat'))
for i = 1:length(dicom_idx)
  C = imfuse(infarct_ex_masked(:,:,i), volume_image(:,:,dicom_idx(i)), 'blend');
  figure();
  imagesc(C)
  colormap gray
  axis equal
  axis off
end
%%
% figure();
% imagesc(volume_image(:,:,dicom_idx(4)));
% colormap gray;
% axis equal;
% axis off;

%%
load(cat(2, base_dir, name, '/', sequence_label{2}, '/MI/MyoInfarct.mat'))
dicom_glob = glob(cat(2, 'C:/Users/ZhangX1/Desktop/contour_exporting_Guan/', name, '/', sequence_label{2}, '/*/*'));
dicom = dicom_glob{1};
[volume_image, slice_data, image_meta_data] = dicom23D(dicom);
[mask, dicom_idx] = ReadMatFile3D(name, sequence_label{2}, anatomys{1});

for i = 1:length(dicom_idx)
  C = imfuse(infarct_ex_masked(:,:,i), volume_image(:,:,dicom_idx(i)), 'blend');
  figure();
  imagesc(C)
  colormap gray
  axis equal
  axis off
end

%%
Names = GetSubjectName(base_dir);
RuleOutLabel = NameRuleOutFunc(Names);
Names = Names(RuleOutLabel == 0);
idx = find(strcmp(name, Names));
optOtsuIndex = xlsread(cat(2, base_dir, 'optOtsuIndex.csv'));
optKmeansIndex = xlsread(cat(2, base_dir, 'optKmeansIndex.csv'));
optGMMIndex = xlsread(cat(2, base_dir, 'optGMMIndex.csv'));
otsuind = optOtsuIndex(idx, 2);
kmeansind = optKmeansIndex(idx, 2);
gmmind = optGMMIndex(idx, 2);

%%
alg = {'Mean5SD', 'Otsu', 'Kmeans', 'GMM'};
alg_ind = [1, 2, 2, 1];
for alg_idx = 1:length(alg)
    alg_label = alg{alg_idx};
    if strcmp(alg_label, 'Mean5SD')
        T1_MI_Path = cat(2, base_dir, name, '\', sequence_label{2}, '\MI\MyoInfarct.mat');
    else
        T1_MI_Path = cat(2, base_dir, name, '\', sequence_label{2}, '\MI_', alg_label, '\MyoInfarct*.mat');
    end
    alg_glob = glob(T1_MI_Path);
    
    load(alg_glob{alg_ind(alg_idx)})
    C = imfuse(infarct_ex_masked(:,:,4), volume_image(:,:,dicom_idx(4)), 'blend');
    figure();
    imagesc(C)
    colormap gray;
    axis off;
    
    figure();
    imagesc(volume_image(:,:,dicom_idx(4)));
    colormap gray;
    axis equal;
    axis off;
     if ~strcmp(alg_label, 'Mean5SD')
         outpath = cat(2, base_dir, name, '\', sequence_label{2}, '\MI_', alg_label, '\MyoInfarct4.tif');
         imwrite(infarct_ex_masked(:,:,4),outpath)
     end
end

%%
figure();
imagesc(mask_myocardium)
%%
alg = {'Mean5SD', 'Otsu', 'Kmeans', 'GMM'};
for alg_idx = 1:length(alg)
    alg_label = alg{alg_idx};
    if strcmp(alg_label, 'Mean5SD')
        T1_MI_Path = cat(2, base_dir, name, '\', sequence_label{2}, '\MI\MyoInfarct.mat');
    else
        T1_MI_Path = cat(2, base_dir, name, '\', sequence_label{2}, '\MI_', alg_label, '\MyoInfarct*.mat');
    end
    alg_glob = glob(T1_MI_Path);
    
    for i = 1:length(alg_glob)
        load(alg_glob{i})
        C = imfuse(infarct_ex_masked(:,:,4), volume_image(:,:,dicom_idx(4)), 'blend');
        figure();
        imagesc(C)
    end
end


%% Bullseye
alg = {'Mean5SD', 'Otsu', 'Kmeans', 'GMM'};
% name = 'OH_SANG_MOO';
% name = 'CHOI_DAE_SUK';
name = 'HAN_BONG_SANG';
% name = 'KIM_KEUM_HUN';
 GetAHABullsEye(name);
for alg_idx = 1:length(alg)
    alg_label = alg{alg_idx};
    GetAHABullsEye(name, alg_label)
end

% A = infarct_ex_masked(:,:,5);
% B = volume_image(:,:,dicom_idx(5));
% alpha = 0.999;
% C = alpha * A + (1 - alpha) * B;
% figure();
% imagesc(C)
% %%
% figure1 = figure;
% ax1 = axes('Parent', figure1);
% ax2 = axes('Parent',figure1);
% set(ax1,'Visible','off');
% set(ax2,'Visible','off');
% [a,map,alpha] = imread('C:\Users\ZhangX1\Documents\MATLAB\masked\OH_SANG_MOO\LGE\MI\MyoInfarct6.tif');
% I = imagesc(a,'Parent',ax2);
% set(I,'AlphaData',alpha);
% imagesc(volume_image(:,:,6),'Parent',ax1);