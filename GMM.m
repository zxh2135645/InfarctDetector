clear all;
close all;

% GMModel = fitgmdist(X,2);
base_dir = 'C:\Users\ZhangX1\Documents\MATLAB\masked\';
addpath('C:\Users\ZhangX1\Documents\MATLAB\cviParser\');
sequence_label = {'LGE', 'T1'};
label = sequence_label{2};
anatomys = {'Heart', 'Myocardium', 'excludeContour', 'MyoReference', 'MI'};
anatomy = anatomys{2};

Names = GetSubjectName(base_dir);
RuleOutLabel = NameRuleOutFunc(Names);
Names = Names(RuleOutLabel == 0);

K = 2:5;
infarct_perc_array = zeros(length(Names), length(K));
optK = zeros(length(Names),1);

for n = 1:length(Names)
    exzmple = Names{n};
    [Myo3D, dicom_idx] = ReadMatFile3D(exzmple, label, anatomy);
    Myo_nnz = Myo3D(:);
    nnz_idx = Myo_nnz ~= 0;
    Myo_nnz = Myo_nnz(nnz_idx);
    infarct_pick = zeros([size(Myo3D), length(K)]);
    
    for clus = 1:length(K)
        GMModel = fitgmdist(Myo_nnz,K(clus));
        idx = cluster(GMModel,Myo_nnz);
        mu = GMModel.mu;
        [sorted_mu, I] = sort(mu);
        clust = zeros(length(idx), K(clus));
        for k = 1:K(clus)
            clust(:, k) = (idx == k);
        end
        
        clustered_gmm = zeros(length(nnz_idx), 1);
        count = 1;
        for i = 1:length(nnz_idx)
            if nnz_idx(i) == 1
                ind = find(clust(count, :) == 1);
                clustered_gmm(i) = ind;
                count = count + 1;
            end
        end
        
        clustered_recon = reshape(clustered_gmm, size(Myo3D));
        clustered_recon_sorted = zeros(size(clustered_recon));
        for i = 1:length(sorted_mu)
            clustered_recon_sorted(clustered_recon == I(i)) = i;
        end
        excludeCtrPath = cat(2, base_dir ,exzmple, '/', label, '/', anatomys{3}, '/');
        if length(ls(excludeCtrPath)) > 2
            load(cat(2, excludeCtrPath, label, '_excludeContour.mat'))
            [excludeContourMat, exctr_idx] = ReadExContour(excludeContour);
            for i = 1:length(dicom_idx)
                ex_ind = find(dicom_idx(i) == exctr_idx);
                if ~isempty(ex_ind)
                    clustered_recon_sorted(:,:,i) = clustered_recon_sorted(:,:,i) .* ~excludeContourMat(:,:,ex_ind);
                end
            end
        end
        infarct_pick(:,:,:,clus) = clustered_recon_sorted == K(clus);
        infarct = ImPostProc(infarct_pick(:,:,:,clus));
        myo_size = sum(Myo3D(:) > 0);
        infarct_perc_array(n, clus) = sum(infarct(:)>0) / myo_size * 100;
        
        infarct_out = cat(2, base_dir, exzmple, '\', label, '\MI_GMM\');
        infarct_ex_masked = infarct;
        if ~ exist(infarct_out, 'dir')
            mkdir(infarct_out);
        end
        
        % Always overwrite
        infarct_out_path = cat(2, infarct_out, 'MyoInfarct', num2str(clus), '.mat');
        save(infarct_out_path, 'infarct_ex_masked');
    end
end

C2 = infarct_perc_array(:,1);
C3 = infarct_perc_array(:,2);
C4 = infarct_perc_array(:,3);
C5 = infarct_perc_array(:,4);

T = table(Names, C2, C3, C4, C5);
out_dir = 'C:/Users/ZhangX1/Documents/MATLAB/masked/';
writetable(T, cat(2, out_dir, 'GMM.csv'));


