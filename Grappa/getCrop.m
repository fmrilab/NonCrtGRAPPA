function crop_fn = getCrop(xSize, pSize)
%% sub-function return a crop_fn function handle for GRAPPA
% For GRAPPA, kcalib is the ACS, when forming the calibration matrix, boundaries
% need to be cutted, since it doesn't have a "full" patch

cSub  = ctrSub(pSize);
npdim = numel(pSize);

cropSubsc = cell(1, npdim);
for ii = 1:npdim % xSize must be larger than pSize along all dim
  cropSubsc{ii} = cSub(ii)+(0:(xSize(ii)-pSize(ii)));
end

crop_fn = @(x)x(cropSubsc{:}, :,:,:,:);

end


