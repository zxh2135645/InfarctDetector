%% Try gap statistic on number of K clusters

clear all
close all

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
        gapstats = evalclusters(Myo3D(:), 'kmeans', 'Gap', 'KList', [2:6]);
        optK(j) = gapstats.OptimalK;
    end
end

% Seems not to work