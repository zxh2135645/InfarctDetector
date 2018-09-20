clear all;
close all;
sequence_label = {'LGE', 'T1'};
base_dir = 'C:/Users/ZhangX1/Desktop/contour_exporting_Guan/*/';
name_glob = glob(base_dir);

%% Get the name
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'\');
    name = strings(end-1);
    Names(i) = name;
    disp(name)
end

RuleOutLabel = NameRuleOutFunc(Names);
Names = Names(RuleOutLabel == 0);

%% Get trigger time and slice location
T1_TriggerTime = zeros(length(Names), 1);
LGE_TriggerTime = zeros(length(Names), 1);
T1_SliceLoc = zeros(length(Names), 1);
LGE_SliceLoc = zeros(length(Names), 1);
T1_dicom_glob = glob(cat(2, base_dir, sequence_label{2}, '/*/*'));
LGE_dicom_glob = glob(cat(2, base_dir, sequence_label{1}, '/*/*'));
T1_dicom_glob = T1_dicom_glob(RuleOutLabel == 0);
LGE_dicom_glob = LGE_dicom_glob(RuleOutLabel == 0);
SliceLocDetailLGE = struct;
SliceLocDetailT1 = struct;

dicom_fields = {... 
        'Filename',...
        'Height', ...
        'Width', ...            
        'Rows',...
        'Columns', ...
        'PixelSpacing',...
        'SliceThickness',...
        'SliceLocation',...
        'SpacingBetweenSlices'...
        'ImagePositionPatient',...
        'ImageOrientationPatient',...
        'MediaStorageSOPInstanceUID',...
        'TriggerTime',...
        };
    
for j = 1:length(Names)
    
    [T1_volume_image, T1_slice_data, T1_image_meta_data] = dicom23D(T1_dicom_glob{j}, dicom_fields);
    [LGE_volume_image, LGE_slice_data, LGE_image_meta_data] = dicom23D(LGE_dicom_glob{j}, dicom_fields);
    T1_tt = zeros(1,length(T1_slice_data));
    T1_SLDetail = zeros(1,length(T1_slice_data));
    LGE_tt = zeros(1,length(LGE_slice_data));
    LGE_SLDetail = zeros(1,length(LGE_slice_data));
    
    for i = 1:length(LGE_slice_data)
        LGE_tt(i) = LGE_slice_data(i).TriggerTime;
        LGE_SLDetail(i) = LGE_slice_data(i).SliceLocation;
    end
    
    for i = 1:length(T1_slice_data)
        T1_tt(i) = T1_slice_data(i).TriggerTime;
        T1_SLDetail(i) = T1_slice_data(i).SliceLocation;
    end
    
    LGE_SliceLoc(j) = LGE_slice_data(1).SliceLocation;
    T1_SliceLoc(j) = T1_slice_data(1).SliceLocation;
    % fprintf('T1 tt: %.4f \n', mean(T1_tt))
    % fprintf('LGE tt: %.4f \n', mean(LGE_tt))
    T1_TriggerTime(j) = mean(T1_tt);
    LGE_TriggerTime(j) = mean(LGE_tt);
    SliceLocDetailLGE.(Names{j}) = LGE_SLDetail;
    SliceLocDetailT1.(Names{j}) = T1_SLDetail;
   
end

%% Set thresholds and labels
diff_TriggerTime_SmallThan200 = abs(LGE_TriggerTime - T1_TriggerTime) < 200;
diff_SliceLoc = round(abs(LGE_SliceLoc - T1_SliceLoc)) == 0;
T = table(Names, LGE_TriggerTime, T1_TriggerTime, diff_TriggerTime_SmallThan200, ...
    LGE_SliceLoc, T1_SliceLoc, diff_SliceLoc);
out_dir = 'C:/Users/ZhangX1/Documents/MATLAB/masked/';
writetable(T, cat(2, out_dir, 'SliceLocNTriggerTime.csv'));
save(cat(2, out_dir, 'SliceLocDetailLGE.mat'), 'SliceLocDetailLGE');
save(cat(2, out_dir, 'SliceLocDetailT1.mat'), 'SliceLocDetailT1');