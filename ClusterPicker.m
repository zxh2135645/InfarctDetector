%% Pick the optimal cluster for Otsu, Kmeans and GMM
clear all;
close all;

base_dir = 'C:\Users\ZhangX1\Documents\MATLAB\masked\';
addpath('C:\Users\ZhangX1\Documents\MATLAB\cviParser\');
sequence_label = {'LGE', 'T1'};
label = sequence_label{2};
anatomys = {'Heart', 'Myocardium', 'excludeContour', 'MyoReference', 'MI'};
algorithms = {'Otsu', 'Kmeans', 'GMM'};
anatomy = anatomys{2};
infarct_t1 = 1470;


name_glob = glob(cat(2, base_dir, '*\'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'\');
    name = strings{end-1};
    Names{i} = name;
end

RuleOutLabel = NameRuleOutFunc(Names);
Names = Names(RuleOutLabel == 0);
[OtsuIndex, ~, ~] = xlsread(fullfile(base_dir, cat(2, 'opt', algorithms{1}, 'Index.csv')));
[KmeansIndex, ~, ~] = xlsread(fullfile(base_dir, cat(2, 'opt', algorithms{2}, 'Index.csv')));
[GMMIndex, ~, ~] = xlsread(fullfile(base_dir, cat(2, 'opt', algorithms{3}, 'Index.csv')));
OtsuPicker = zeros(length(Names), 1);
KmeansPicker = zeros(length(Names), 1);
GMMPicker = zeros(length(Names), 1);

for i = 1:length(Names)
    name = Names{i};
    disp(name);
    
    myo_glob = glob(fullfile(base_dir, name, label, anatomys{2}, '*'));
    for k = 1:length(myo_glob)
        load(myo_glob{k});
        if k > 1
            mask_3d = cat(3, mask_3d, mask_myocardium);
        elseif k == 1
            mask_3d = mask_myocardium;    
        end
    end
    
    for j = 1:length(algorithms)
        alg = algorithms{j};
        disp(alg);
        mi_glob = glob(fullfile(base_dir, name, label, cat(2,anatomys{5}, '_', alg), '*'));
        
        infarct_mean_array = zeros(length(mi_glob), 1);
        for l = 1:length(mi_glob)
            load(mi_glob{l});
            mapped_infarct = mask_3d .* infarct_ex_masked;
            infarct_mean = mean(nonzeros(mapped_infarct(:)));
            disp(infarct_mean)
            infarct_mean_array(l) = infarct_mean;
%             figure();
%             n = ceil(sqrt(size(mapped_infarct, 3)));
%             for m = 1:size(mapped_infarct, 3)
%                 subplot(n,n,m)
%                 imagesc(mapped_infarct(:,:,m))
%             end
        end
        [M, I] = min(abs(infarct_mean_array - infarct_t1));
        switch alg
            case algorithms{1}
                % disp(infarct_mean_array(OtsuIndex(i, 2)))
                OtsuPicker(i) = I;
            case algorithms{2}
                % disp(infarct_mean_array(KmeansIndex(i, 2)))
                KmeansPicker(i) = I;
            case algorithms{3}
                % disp(infarct_mean_array(GMMIndex(i, 2)))
                GMMPicker(i) = I;
        end
    end
end
T = table(Names, OtsuPicker, KmeansPicker, GMMPicker);
writetable(T, fullfile(base_dir, 'ClusterPicker1.csv'));