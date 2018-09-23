clear all;
close all;

base_dir = 'C:\Users\ZhangX1\Documents\MATLAB\masked\';
Names = GetSubjectName(base_dir);
RuleOutLabel = zeros(length(Names), 1);

for i = 1:length(Names)
    switch Names{i}
        case {'BOO_SEUNG_HYUN', 'LIM_DAE_YUNG', 'KIM_YUNG', 'CHUNG_SAM_YUN', 'LI_SHU_YONG'}
            disp(Names{i});
            RuleOutLabel(i) = 1;
        otherwise 
    end
end

T = table(Names, RuleOutLabel);
out_dir = base_dir;
writetable(T, cat(2, out_dir, 'NamesToRuleOut.csv'));
