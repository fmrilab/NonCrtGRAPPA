function mask = normBallMask(bSize, p)
if ~exist('p', 'var'), p = 2; end
if isinf(p), mask = true(bSize); return; end

ndim = numel(bSize);
cSub = shiftdim(ctrSub(bSize), -ndim+1);

coord_c = cellfun(@(x)(1:x).',num2cell(bSize), 'UniformOutput',false);
grid_c  = cell(1,numel(bSize));

[grid_c{:}] = ndgrid(coord_c{:});
grids = bsxfun(@rdivide, cat(ndim+1, grid_c{:}), cSub) - 1;

mask = sum(abs(grids).^p, ndim+1).^(1/p) <= (0.99+eps);

end

