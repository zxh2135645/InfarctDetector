function [excludeContourMat, exctr_idx] = ReadExContour(excludeContour)
    % excludeContour is a struct
fname = fieldnames(excludeContour);
if ~isempty(struct2table(excludeContour))
    excludeContourMat = excludeContour.(fname{1}){1};
    exctr_idx = excludeContour.(fname{1}){2};
else
    excludeContourMat = [];
    exctr_idx = [];
    
end
end