clear all;
close all;

base_dir = 'C:\Users\ZhangX1\Documents\MATLAB\masked\';
addpath('C:\Users\ZhangX1\Documents\MATLAB\cviParser\');
sequence_label = {'LGE', 'T1'};
label = sequence_label{1};
anatomys = {'Heart', 'Myocardium', 'excludeContour', 'MyoReference', 'MI'};
anatomy = anatomys{2};

Names = GetSubjectName(base_dir);
RuleOutLabel = NameRuleOutFunc(Names);
Names = Names(RuleOutLabel == 0);
Accuracy = zeros(length(Names), 1);
Precision = zeros(length(Names), 1);
Sensitivity = zeros(length(Names), 1);
Specificity = zeros(length(Names), 1);


for name_idx = 1:length(Names)
    name = Names{name_idx};
    LGE_MI_Path = cat(2, base_dir, name, '\', label, '\MI\MyoInfarct.mat');
    load(LGE_MI_Path);
    
    load(cat(2, base_dir, 'SliceLocDetailLGE.mat'));
    load(cat(2, base_dir, 'SliceLocDetailT1.mat'));
    sliceLocLGE = SliceLocDetailLGE.(name);
    sliceLocT1 = SliceLocDetailT1.(name);
    [LGE_Myo3D, LGE_dicom_idx] = ReadMatFile3D(name, sequence_label{1}, anatomy);
    [T1_Myo3D, T1_dicom_idx] = ReadMatFile3D(name, sequence_label{2}, anatomy);
    sliceLocLGE_contoured = sliceLocLGE(LGE_dicom_idx);
    sliceLocT1_contoured = sliceLocT1(T1_dicom_idx);
        
    clear LGE_idx_mapped T1_idx_mapped
    for i = 1:length(sliceLocT1_contoured)
        idx = find(abs(sliceLocT1_contoured(i) - sliceLocLGE_contoured) <= 4); % Slice Thickness is 8mm
        if ~isempty(idx)
            LGE_idx_mapped(i) = min(idx);
            T1_idx_mapped(i) = T1_dicom_idx(i);
        end
    end
    
    LGE_idx_mapped = nonzeros(LGE_idx_mapped);
    T1_idx_mapped = nonzeros(T1_idx_mapped);
    
    if length(LGE_idx_mapped) == length(T1_idx_mapped)
        n = length(LGE_idx_mapped);
        mode = mod(n,3);
        switch mode
            case {0, 1}
                integ = fix(n/3);
                aha_slice = [2; 2+integ; 2+integ*2];
            case {2}
                integ = fix(n/3) + 1;
                aha_slice = [1; 1+integ; 1+integ*2];
        end
    end
    
    LGE_MyoMask = LGE_Myo3D > 0;
    LGE_SegPixCount = cell(3, 1);
    
    for i = 1:length(aha_slice)
        if i == 1 || i == 2
            LocPixCount = zeros(6, 1);
            [Segmentpix, stats] = AHASegmentation(infarct_ex_masked(:,:,aha_slice(i)), LGE_MyoMask(:,:,aha_slice(i)),6,0);
        elseif i == 3
            LocPixCount = zeros(4, 1);
            [Segmentpix, stats] = AHASegmentation(infarct_ex_masked(:,:,aha_slice(i)), LGE_MyoMask(:,:,aha_slice(i)),4,45);
        end
        
        for j = 1:length(Segmentpix)
            LocPixCount(j) = sum(Segmentpix{j});
        end
        LGE_SegPixCount{i} = LocPixCount;
    end
    
    T1_MyoMask = T1_Myo3D > 0;
    T1_SegPixCount = cell(3, 1);
    T1_MI_Path = cat(2, base_dir, name, '\', sequence_label{2}, '\MI\MyoInfarct.mat');
    load(T1_MI_Path);
    for i = 1:length(aha_slice)
        if i == 1 || i == 2
            LocPixCount = zeros(6, 1);
            [Segmentpix, stats] = AHASegmentation(infarct_ex_masked(:,:,aha_slice(i)), T1_MyoMask(:,:,aha_slice(i)),6,0);
        elseif i == 3
            LocPixCount = zeros(4, 1);
            [Segmentpix, stats] = AHASegmentation(infarct_ex_masked(:,:,aha_slice(i)), T1_MyoMask(:,:,aha_slice(i)),4,45);
        end
        
        for j = 1:length(Segmentpix)
            LocPixCount(j) = sum(Segmentpix{j});
        end
        T1_SegPixCount{i} = LocPixCount;
    end
    
    LGE_SegPixCal = zeros(16, 1);
    T1_SegPixCal = zeros(16, 1);
    for i = 1:length(aha_slice)
        if i == 1
            LGE_SegPixCal(1:6) = LGE_SegPixCount{i} > 0;
            T1_SegPixCal(1:6) = T1_SegPixCount{i} > 0;
        elseif i == 2
            LGE_SegPixCal(7:12) = LGE_SegPixCount{i} > 0;
            T1_SegPixCal(7:12) = T1_SegPixCount{i} > 0;
        elseif i == 3
            LGE_SegPixCal(13:16) = LGE_SegPixCount{i} > 0;
            T1_SegPixCal(13:16) = T1_SegPixCount{i} > 0;
        end
    end
    
    pos = find(T1_SegPixCal == 1);
    TP = sum(LGE_SegPixCal(pos));
    FP = sum(~LGE_SegPixCal(pos));
    neg = find(T1_SegPixCal == 0);
    TN = sum(~LGE_SegPixCal(neg));
    FN = sum(LGE_SegPixCal(neg));
    
  
    Sensitivity(name_idx) = TP / (TP + FN);
    Specificity(name_idx) = TN / (TN + FP);
    Accuracy(name_idx) = (TP + TN) / (TP + TN + FN + FP);
    Precision(name_idx) = TP / (TP + FP);
end

T = table(Names, Sensitivity, Specificity, Accuracy, Precision);
out_dir = 'C:/Users/ZhangX1/Documents/MATLAB/masked/';
%writetable(T, cat(2, out_dir, 'rocMean5SD.csv'));
writetable(T, cat(2, out_dir, 'rocOtsu.csv'));