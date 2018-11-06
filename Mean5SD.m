%% To get Mean + 5SD as ground truth 
% This version is for clustering methods
clear all;
close all;

addpath('C:\Users\ZhangX1\Documents\MATLAB\cviParser\');

sequence_label = {'LGE', 'T1'};
anatomy_label = {'Heart', 'Myocardium', 'excludeContour', 'MyoReference'};
anatomy = anatomy_label{2};
base_dir = 'C:/Users/ZhangX1/Documents/MATLAB/masked/';

name_glob = glob(cat(2, base_dir, '/*_*'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'\');
    name = strings{end-1};
    Names{i} = name;
end

RuleOutLabel = NameRuleOutFunc(Names);
Names = Names(RuleOutLabel == 0);

for la = 1:length(sequence_label)
    infarct_perc_struct = struct;
    label = sequence_label{la};
    infarct_perc_array = zeros(size(Names));
    for i = 1:length(Names)
        name = Names{i};
        disp(name)
        
        mask = [];
        excludeContour = [];
        myoRefCell = [];
        
        if length(ls(cat(2, base_dir, name, '/', label, '/', anatomy))) > 2
            [mask, dicom_idx] = ReadMatFile3D(name, label, anatomy);
        end
        if length(ls(cat(2, base_dir, name, '/', label, '/', anatomy_label{3}))) > 2
            load(cat(2, base_dir, name, '/', label, '/', anatomy_label{3}, '/', label, '_excludeContour.mat'));
        end
        if length(ls(cat(2, base_dir, name, '/', label, '/', anatomy_label{4}))) > 2
            load(cat(2, base_dir, name, '/', label, '/', anatomy_label{4}, '/', label, '_myoRef.mat'));
        end
        if isempty(mask) || isempty(myoRefCell)
            disp('break')
            % infarct_perc_struct.(name) = [];
            infarct_perc_array(i) = nan;
        else
            binarize_mask = mask > 0;
            fname = fieldnames(excludeContour);
            if ~isempty(struct2table(excludeContour))
                excludeContourMat = excludeContour.(fname{1}){1};
                exctr_idx = excludeContour.(fname{1}){2};
            else
                exctr_idx = [];
            end
            myoRefMat = myoRefCell{1};
            refctr_idx = myoRefCell{2};
            compositeIm = zeros(size(mask));
            
            for j = 1:length(dicom_idx)
                
                ex_ind = find(dicom_idx(j) == exctr_idx);
                ref_ind = find(dicom_idx(j) == refctr_idx);
                if ~isempty(ref_ind)
                    compositeIm(:,:,j) = binarize_mask(:,:,j) + myoRefMat(:,:,ref_ind);
                else
                    compositeIm(:,:,j) = binarize_mask(:,:,j);
                end
                
                if ~isempty(ex_ind)
                    compositeIm(:,:,j) = compositeIm(:,:,j) + excludeContourMat(:,:,ex_ind) * (-3);
                end
            end
            %        % Display Images
            %                    figure();
            %                    n = ceil(sqrt(size(mask, 3)));
            %                    for k = 1: size(mask, 3)
            %                        subplot(n,n,k)
            %                        imagesc(compositeIm(:,:,k))
            %                        axis equal
            %                    end
            Outpath = cat(2, base_dir, name, '/', label, '/', 'compositeMat.mat');
            if ~exist(Outpath, 'file')
                save(Outpath, 'compositeIm');
            end
            
            interp_idx = zeros(size(dicom_idx));
            for j = 1:length(dicom_idx)
                [minDistance, indexOfMin] = min(abs(refctr_idx - dicom_idx(j)));
                interp_idx(j) = refctr_idx(indexOfMin);
            end
            
            infarct_ex_masked = zeros(size(mask));
            for j = 1:length(interp_idx)
                refIm = (compositeIm(:,:,j) == -1) | (compositeIm(:,:,j) == 2);
                exIm = compositeIm(:,:,j) < 0;
                ref_ind = find(interp_idx(j) == dicom_idx);
                
                if any(refIm(:))
                    refMasked = refIm .* mask(:,:,j);
                else
                    refIm = (compositeIm(:,:,ref_ind) == -1) | (compositeIm(:,:,ref_ind) == 2);
                    refMasked = refIm .* mask(:,:,ref_ind);
                end
                
                ref_mean = mean(nonzeros(refMasked(:)));
                ref_sd = std(nonzeros(refMasked(:)));
                thresh = ref_mean + 5*ref_sd;
                infarct_raw = mask(:,:,j) > thresh;
                infarct_ex_masked(:,:,j) = infarct_raw .* (~exIm);
            end
            myo_size = sum(mask(:) > 0);
            
            infarct_ex_masked = ImPostProc(infarct_ex_masked);
            infarct_size = sum(infarct_ex_masked(:) > 0);
            infarct_perc = infarct_size / myo_size * 100;
            disp(infarct_perc);
            infarct_perc_array(i) = infarct_perc;
            
            infarct_out = cat(2, base_dir, name, '/', label, '/MI/');
            if ~ exist(infarct_out, 'dir')
                mkdir(infarct_out);
            end
         
            % Always overwrite
            infarct_out_path = cat(2, infarct_out, 'MyoInfarct.mat');
            save(infarct_out_path, 'infarct_ex_masked');
            
            % Always overwrite
            for j = 1:length(interp_idx)
                infarct_out_path = cat(2, infarct_out, 'MyoInfarct', num2str(j),'.tif');
                imwrite(infarct_ex_masked(:,:,j), infarct_out_path);
            end
        end
    end
    if strcmp(label, 'LGE')
        LGE_InfarctPerc = infarct_perc_array;
    elseif strcmp(label, 'T1')
        T1_InfarctPerc = infarct_perc_array;
    end
end

T = table(Names, LGE_InfarctPerc, T1_InfarctPerc);
writetable(T, cat(2, base_dir, 'InfarctPercMean5SD.csv'));