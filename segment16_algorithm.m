%% This version is initially used but now aborted
% Prefered to use the segment16_volume_algorithm.m to volume-by-volume
% segment
clear all;
close all;

base_dir = 'C:\Users\ZhangX1\Documents\MATLAB\masked\';
addpath('C:\Users\ZhangX1\Documents\MATLAB\cviParser\');
sequence_label = {'LGE', 'T1'};
anatomys = {'Heart', 'Myocardium', 'excludeContour', 'MyoReference', 'MI'};
anatomy = anatomys{2};

alg = {'Mean5SD', 'Otsu', 'Kmeans', 'GMM'};

CoordsFileName_LGE = [base_dir, sequence_label{1}, '_coords.csv'];
[num_lge,txt,raw_lge] = xlsread(CoordsFileName_LGE);
CoordsNames_LGE = txt(2:end,1);

CoordsFileName_T1 = [base_dir, sequence_label{2}, '_coords.csv'];
[num_t1,txt,raw_t1] = xlsread(CoordsFileName_T1);
CoordsNames_T1 = txt(2:end,1);

if length(CoordsNames_LGE) ~= length(CoordsNames_T1)
    error('Number of LGE and T1 coords are different, please check!');
end

Names = GetSubjectName(base_dir);
RuleOutLabel = NameRuleOutFunc(Names);
Names = Names(RuleOutLabel == 0);
LGE_aha = zeros(length(Names) * 16, 1);
T1_aha = zeros(length(Names) * 16, length(alg));

Coords_idx = zeros(length(CoordsNames_LGE), 1);
for i = 1:length(CoordsNames_LGE)
    if ~isempty(find(strcmp(CoordsNames_LGE{i}, Names), 1))
        Coords_idx(i) = find(strcmp(CoordsNames_LGE{i}, Names));
    end
end

for alg_iter = 1:length(alg)
    alg_label = alg{alg_iter};
    TP = 0;
    TN = 0;
    FP = 0;
    FN = 0;
    if ~strcmp(alg_label, 'Mean5SD')
        Index = xlsread(cat(2, base_dir, 'opt', alg_label, 'Index.csv'));
        Index = Index(:,2);
    end
    
    for name_idx = 1:length(Names)
        name = Names{name_idx};
        LGE_MI_Path = cat(2, base_dir, name, '\', sequence_label{1}, '\MI\MyoInfarct.mat');
        load(LGE_MI_Path);
        
        load(cat(2, base_dir, 'SliceLocDetailLGE.mat'));
        load(cat(2, base_dir, 'SliceLocDetailT1.mat'));
        sliceLocLGE = SliceLocDetailLGE.(name);
        sliceLocT1 = SliceLocDetailT1.(name);
        
        [LGE_Myo3D, LGE_dicom_idx] = ReadMatFile3D(name, sequence_label{1}, anatomy);
        [T1_Myo3D, T1_dicom_idx] = ReadMatFile3D(name, sequence_label{2}, anatomy);
        
        LGE_original_idx = 1:1:length(sliceLocLGE);
        T1_original_idx = 1:1:length(sliceLocT1);
        sliceLocLGE_contoured = sliceLocLGE(LGE_dicom_idx);
        sliceLocT1_contoured = sliceLocT1(T1_dicom_idx);
        
                
        if min(abs(sliceLocLGE_contoured)) == abs(sliceLocLGE_contoured(end))
            LGE_Myo3D = LGE_Myo3D(:,:,end:-1:1);
            infarct_ex_masked = infarct_ex_masked(:,:,end:-1:1);
        end
        
        % Read the reference coordinates
        refP_idx = find(name_idx == Coords_idx);
        x = num_lge(refP_idx, 1);
        y = num_lge(refP_idx, 2);
        x_centroid = num_lge(refP_idx, 3);
        y_centroid = num_lge(refP_idx, 4);
        BaseGroove = atan2(x - x_centroid, y - y_centroid) * 180 / pi;
        
        clear LGE_idx_mapped T1_idx_mapped
        for i = 1:length(sliceLocT1_contoured)
            idx = find(abs(sliceLocT1_contoured(i) - sliceLocLGE_contoured) <= 4); % Slice Thickness is 8mm
            if ~isempty(idx)
                LGE_idx_mapped(i) = LGE_dicom_idx(min(idx));
                T1_idx_mapped(i) = T1_dicom_idx(i);
            end
        end
        
        LGE_idx_mapped = nonzeros(LGE_idx_mapped);
        T1_idx_mapped = nonzeros(T1_idx_mapped);
        
        if length(LGE_idx_mapped) == length(T1_idx_mapped)
            n = length(LGE_idx_mapped);
            mode = mod(n,3);
            if n > 3
                switch mode
                    case {0, 1}
                        integ = fix(n/3);
                        aha_slice = [2; 2+integ; 2+integ*2];
                    case {2}
                        integ = fix(n/3) + 1;
                        aha_slice = [1; 1+integ; 1+integ*2];
                end
            elseif n ==3
                aha_slice = [1 2 3];
            else
                error("Available slice numbers are smaller than 3.");
            end
        end
        
        LGE_aha_seg_slice = zeros(length(aha_slice), 1);
        T1_aha_seg_slice = zeros(length(aha_slice), 1);
        for i = 1:length(aha_slice)
            LGE_aha_seg_slice(i) = find(LGE_dicom_idx == LGE_idx_mapped(aha_slice(i)));
            T1_aha_seg_slice(i) = find(T1_dicom_idx == T1_idx_mapped(aha_slice(i)));
        end
        
        
        LGE_MyoMask = LGE_Myo3D > 0;
        LGE_SegPixCount = cell(3, 1);
        LGE_SegTotalPixCount = cell(3, 1);
        
        for i = 1:length(aha_slice)
            if i == 1 || i == 2
                LocPixCount = zeros(6, 1);
                SegTotalPixCount = zeros(6, 1);
                Groove = BaseGroove + 60;
                [Segmentpix, stats] = AHASegmentation(infarct_ex_masked(:,:,LGE_aha_seg_slice(i)), LGE_MyoMask(:,:,LGE_aha_seg_slice(i)),6, Groove);
            elseif i == 3
                LocPixCount = zeros(4, 1);
                SegTotalPixCount = zeros(4, 1);
                Groove = BaseGroove + 75;
                [Segmentpix, stats] = AHASegmentation(infarct_ex_masked(:,:,LGE_aha_seg_slice(i)), LGE_MyoMask(:,:,LGE_aha_seg_slice(i)),4, Groove);
            end
            
            for j = 1:length(Segmentpix)
                LocPixCount(j) = sum(Segmentpix{j});
                SegTotalPixCount(j) = length(Segmentpix{j});
            end
            LGE_SegPixCount{i} = LocPixCount;
            LGE_SegTotalPixCount{i} = SegTotalPixCount;
        end
        
        if strcmp(alg_label, 'Mean5SD')
            T1_MI_Path = cat(2, base_dir, name, '\', sequence_label{2}, '\MI\MyoInfarct.mat');
        else
            T1_MI_Path = cat(2, base_dir, name, '\', sequence_label{2}, '\MI_', alg_label, '\MyoInfarct', num2str(Index(name_idx)), '.mat');
        end
        load(T1_MI_Path);
        
        if min(abs(sliceLocT1_contoured)) == abs(sliceLocT1_contoured(end))
            T1_Myo3D = T1_Myo3D(:,:,end:-1:1);
            infarct_ex_masked = infarct_ex_masked(:,:,end:-1:1);
        end
        
        T1_MyoMask = T1_Myo3D > 0;
        T1_SegPixCount = cell(3, 1);
        T1_SegTotalPixCount = cell(3, 1);
             
        % Read the reference coordinates
        x = num_t1(refP_idx, 1);
        y = num_t1(refP_idx, 2);
        x_centroid = num_t1(refP_idx, 3);
        y_centroid = num_t1(refP_idx, 4);
        BaseGroove = atan2(x - x_centroid, y - y_centroid) * 180 / pi;
        
        for i = 1:length(aha_slice)
            if i == 1 || i == 2
                LocPixCount = zeros(6, 1);
                SegTotalPixCount = zeros(6, 1);
                Groove = BaseGroove + 60;
                [Segmentpix, stats] = AHASegmentation(infarct_ex_masked(:,:,T1_aha_seg_slice(i)), T1_MyoMask(:,:,T1_aha_seg_slice(i)),6,Groove);
            elseif i == 3
                LocPixCount = zeros(4, 1);
                SegTotalPixCount = zeros(4, 1);
                Groove = BaseGroove + 75;
                [Segmentpix, stats] = AHASegmentation(infarct_ex_masked(:,:,T1_aha_seg_slice(i)), T1_MyoMask(:,:,T1_aha_seg_slice(i)),4,Groove);
            end
            
            for j = 1:length(Segmentpix)
                LocPixCount(j) = sum(Segmentpix{j});
                SegTotalPixCount(j) = length(Segmentpix{j});
            end
            T1_SegPixCount{i} = LocPixCount;
            T1_SegTotalPixCount{i} = SegTotalPixCount;
        end
        
        LGE_SegPixCal = zeros(16, 1);
        T1_SegPixCal = zeros(16, 1);
        LGE_SegPixPerc = zeros(16, 1);
        T1_SegPixPerc = zeros(16, 1);
        for i = 1:length(aha_slice)
            if i == 1
                LGE_SegPixCal(1:6) = LGE_SegPixCount{i} ./ LGE_SegTotalPixCount{i} * 100 > 1;
                T1_SegPixCal(1:6) = T1_SegPixCount{i} ./ T1_SegTotalPixCount{i} * 100 > 1;
                LGE_SegPixPerc(1:6) = LGE_SegPixCount{i} ./ LGE_SegTotalPixCount{i} * 100;
                T1_SegPixPerc(1:6) = T1_SegPixCount{i} ./ T1_SegTotalPixCount{i} * 100;
            elseif i == 2
                LGE_SegPixCal(7:12) = LGE_SegPixCount{i} ./ LGE_SegTotalPixCount{i} * 100 > 1;
                T1_SegPixCal(7:12) = T1_SegPixCount{i} ./ T1_SegTotalPixCount{i} * 100 > 1;
                LGE_SegPixPerc(7:12) = LGE_SegPixCount{i} ./ LGE_SegTotalPixCount{i} * 100;
                T1_SegPixPerc(7:12) = T1_SegPixCount{i} ./ T1_SegTotalPixCount{i} * 100;
            elseif i == 3
                LGE_SegPixCal(13:16) = LGE_SegPixCount{i} ./ LGE_SegTotalPixCount{i} * 100 > 1;
                T1_SegPixCal(13:16) = T1_SegPixCount{i} ./ T1_SegTotalPixCount{i} * 100 > 1;
                LGE_SegPixPerc(13:16) = LGE_SegPixCount{i} ./ LGE_SegTotalPixCount{i} * 100;
                T1_SegPixPerc(13:16) = T1_SegPixCount{i} ./ T1_SegTotalPixCount{i} * 100;
            end
        end
        
        pos = find(T1_SegPixCal == 1);
        TP = TP + sum(LGE_SegPixCal(pos));
        FP = FP + sum(~LGE_SegPixCal(pos));
        neg = find(T1_SegPixCal == 0);
        TN = TN + sum(~LGE_SegPixCal(neg));
        FN = FN + sum(LGE_SegPixCal(neg));
        
        if alg_iter == 1
            LGE_aha((16*name_idx - 15):(16*name_idx), alg_iter) = LGE_SegPixPerc;
        end
        T1_aha((16*name_idx - 15):(16*name_idx), alg_iter) = T1_SegPixPerc;
    end
    
    Sensitivity = TP / (TP + FN);
    Specificity = TN / (TN + FP);
    Accuracy = (TP + TN) / (TP + TN + FN + FP);
    Precision = TP / (TP + FP);
    
    T = table(alg_label, Sensitivity, Specificity, Accuracy, Precision, TP, TN, FP, FN);
    out_dir = 'C:/Users/ZhangX1/Documents/MATLAB/masked/';
    writetable(T, cat(2, out_dir, 'roc', alg_label,'.csv'));
end

%%
alg = {'Mean5SD', 'Otsu', 'Kmeans', 'GMM'};
LGE_Mean5SD = LGE_aha;
T1_Mean5SD = T1_aha(:,1);
T1_Otsu = T1_aha(:,2);
T1_Kmeans = T1_aha(:,3);
T1_GMM = T1_aha(:,4);

T2 = table(LGE_Mean5SD, T1_Mean5SD, T1_Otsu, T1_Kmeans, T1_GMM);
writetable(T2, cat(2, out_dir, 'rocForReal.csv'));