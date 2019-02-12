function mloc = mask2mloc(mask)
% the mloc returned has mask as its 1st col

dims = size(mask);
if iscolumn(mask), dims = dims(1); end
ndim = numel(dims);

cSub = ctrSub(dims);
inds_c = cell(1, ndim);

[inds_c{:}] = ind2sub(dims, (1:numel(mask)).');
loc = bsxfun(@minus, cell2mat(inds_c), cSub(:)');
mloc = [~~mask(:), loc];

end

