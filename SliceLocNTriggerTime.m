clear all;
close all;
sequence_label = {'LGE', 'T1'};
label = char(sequence_label(1));

base_dir = 'C:/Users/ZhangX1/Desktop/contour_exporting_Guan/*/';
dicom_glob = glob(cat(2, base_dir, label, '/*/*'));

%% Get the name
name_cell = cell(1, length(dicom_glob));
for i = 1:length(dicom_glob)
    strings = strsplit(dicom_glob{i},'\');
    name = strings(end-4);
    name_cell(i) = name;
    disp(name)
end

%% Get trigger time and slice location
T1_TriggerTime = zeros(length(name_cell), 1);
LGE_TriggerTime = zeros(length(name_cell), 1);
T1_SliceLoc = zeros(length(name_cell), 1);
LGE_SliceLoc = zeros(length(name_cell), 1);
T1_dicom_glob = glob(cat(2, base_dir, sequence_label{1}, '/*/*'));
LGE_dicom_glob = glob(cat(2, base_dir, sequence_label{2}, '/*/*'));

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
    
for j = 1:length(name_cell)
    display(name_cell{i});
    
    [T1_volume_image, T1_slice_data, T1_image_meta_data] = dicom23D(T1_dicom_glob{j}, dicom_fields);
    [LGE_volume_image, LGE_slice_data, LGE_image_meta_data] = dicom23D(LGE_dicom_glob{j}, dicom_fields);
    T1_tt = zeros(1,length(T1_slice_data));
    LGE_tt = zeros(1,length(LGE_slice_data));
    for i = 1:length(LGE_slice_data)
        LGE_tt(i) = LGE_slice_data(i).TriggerTime;
    end
    
    for i = 1:length(T1_slice_data)
        T1_tt(i) = T1_slice_data(i).TriggerTime;
    end
    LGE_SliceLoc(j) = LGE_slice_data(1).SliceLocation;
    T1_SliceLoc(j) = T1_slice_data(1).SliceLocation;
    % fprintf('T1 tt: %.4f \n', mean(T1_tt))
    % fprintf('LGE tt: %.4f \n', mean(LGE_tt))
    T1_TriggerTime(j) = mean(T1_tt);
    LGE_TriggerTime(j) = mean(LGE_tt);
end

%% Set thresholds and labels
diff_TriggerTime_SmallThan100 = abs(LGE_TriggerTime - T1_TriggerTime) < 100;
diff_SliceLoc = round(abs(LGE_SliceLoc - T1_SliceLoc)) == 0;
Names = name_cell';
T = table(Names, LGE_TriggerTime, T1_TriggerTime, diff_TriggerTime_SmallThan100, ...
    LGE_SliceLoc, T1_SliceLoc, diff_SliceLoc);
out_dir = 'C:/Users/ZhangX1/Documents/MATLAB/masked/';
writetable(T, cat(2, out_dir, 'SliceLocNTriggerTime.csv'));
