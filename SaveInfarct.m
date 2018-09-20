function SaveInfarct(infarct_out, infarct_ex_masked)
interp_idx = size(infarct_ex_masked);
if ~ exist(infarct_out, 'dir')
    mkdir(infarct_out);
end
if length(ls(infarct_out)) == 2
    % Not always overwrite
    infarct_out_path = cat(2, infarct_out, 'MyoInfarct.mat');
    save(infarct_out_path, 'infarct_ex_masked');
end
if length(ls(infarct_out)) > 2
    % Not always overwrite
    for j = 1:length(interp_idx)
        infarct_out_path = cat(2, infarct_out, 'MyoInfarct', num2str(j),'.tif');
        imwrite(infarct_ex_masked(:,:,j), infarct_out_path);
    end
end
end