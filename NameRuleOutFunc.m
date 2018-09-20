function RuleOutLabel = NameRuleOutFunc(Names)
RuleOutLabel = zeros(length(Names), 1);

for i = 1:length(Names)
    switch Names{i}
        case {'BOO_SEUNG_HYUN', 'LIM_DAE_YUNG', 'KIM_YUNG', 'CHUNG_SAM_YUN', 'LI_SHU_YONG'}
            RuleOutLabel(i) = 1;
        otherwise
            RuleOutLabel(i) = 0;
    end
end
end