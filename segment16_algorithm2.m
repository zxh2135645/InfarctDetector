clear all;
close all;

base_dir = 'C:\Users\ZhangX1\Documents\MATLAB\masked\';
addpath('C:\Users\ZhangX1\Documents\MATLAB\cviParser\');
sequence_label = {'LGE', 'T1'};
anatomys = {'Heart', 'Myocardium', 'excludeContour', 'MyoReference', 'MI'};
anatomy = anatomys{2};

alg = {'Mean5SD', 'Otsu', 'Kmeans', 'GMM'};

CoordsFileName_LGE = [base_dir, sequence_label{1}, '_coords2.mat'];
LGE_coords = load(CoordsFileName_LGE);

CoordsFileName_T1 = [base_dir, sequence_label{2}, '_coords2.mat'];
T1_coords = load(CoordsFileName_T1);

Names = GetSubjectName(base_dir);
RuleOutLabel = NameRuleOutFunc(Names);
Names = Names(RuleOutLabel == 0);

CoordsFileName_LGE = [base_dir, sequence_label{1}, '_coords.csv'];
[num_lge,txt,raw_lge] = xlsread(CoordsFileName_LGE);
CoordsNames_LGE = txt(2:end,1);

Coords_idx = zeros(length(CoordsNames_LGE), 1);
for i = 1:length(CoordsNames_LGE)
    if ~isempty(find(strcmp(CoordsNames_LGE{i}, Names), 1))
        Coords_idx(i) = find(strcmp(CoordsNames_LGE{i}, Names));
    end
end

