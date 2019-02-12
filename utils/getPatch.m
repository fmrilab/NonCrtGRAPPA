function [patch_c] = getPatch(sMask_sTraj, pSize, p, method)
% kd-tree finding ngbs, then group into patches
%INPUTS
% - sMask_sTraj
% - pSize
% - p (1,): minkowski dist parameter, 1-cityblock, 2-euclidean, inf-chebychev
% - method str: 'kdTree', 'grid', etc.
%OUTPUTS
% - patch_c (#_patch_types, 2) cell: take iith row as an example
%     (ii,1):
%       sMask case: (#ngb, 1),   mask identifying sampled ngbs;
%       sTraj case: (#ngb,ndim), ngb-ctr shifts, rowsort()'ed along dimensions
%     (ii,2):
%       array composed as [ctrInds, ngbInds];

%% Prep

if islogical(sMask_sTraj)||isinteger(sMask_sTraj)
  sMask = sMask_sTraj;
  if ~exist('method','var'), method = 'grid'; end
  switch lower(method)
    case 'grid', patch_c = patch_grid(sMask, pSize, p);
    otherwise, error('not supported for sMask');
  end
elseif isfloat(sMask_sTraj)
  sTraj = sMask_sTraj;
  if ~exist('method','var'), method = 'kdtree'; end
  switch lower(method)
    case 'kdtree', patch_c = patch_kdTree(sTraj, pSize, p);
    otherwise, error('not supported for sTraj');
  end
else, error('not supported');
end

end

%% spcialized for sMask
function [patch_c] = patch_grid(sMask, pSize, p)
pMask = normBallMask(pSize, p); % patch mask, containing both ctr and ngb
pMask(ctrInd(pSize)) = false;
% pcSub = ctrSub(pSize);

nindsMask = embed((1:nnz(sMask)).', ~~sMask);

ninds = patch2row(nindsMask, pMask, 'pad0', ctrSub(pSize));
ninds(~~nindsMask(:), :) = []; % remove where ctr is sampled

pattern = ~~ninds; % logical
inds = [(1:size(ninds,1)).', ninds]; % add patch center indices into inds

inds_c = cellfun(@(x)x(~~x), num2cell(inds, 2),'UniformOutput',false);

patch_c = [num2cell(pattern.',1).', inds_c];
end

%% specialized for sTraj
function [patch_c] = patch_kdTree(sTraj, pSize, p)
trajMask = ~~sTraj(:,1);
[traj0, traj1] = deal(sTraj(~trajMask,2:end), sTraj(trajMask,2:end));

inds_c = inds_kdTree(traj0, traj1, pSize, p);

% row sort shifts, so duplicate patterns can be identified easily later
[n0s, ndim] = size(traj0);
shifts_rc = cell(n0s, 1); % _c: cell, _r: raw
colOrder = ndim:-1:1;
parfor ii = 1:n0s
  % calc shifts w/in a patch and keep them in order
  nInds_i = inds_c{ii}(2:end); % 1st entry is ctr, not ngb
  if isempty(nInds_i), continue; end
  shifts_r_i = bsxfun(@minus, traj1(nInds_i,:), traj0(ii,:)); % _i: iter
  
  % sortrows orders ngbs, thus easier to be processed beyond this function
  shifts_i_int = round(shifts_r_i);
  [~, inds_s]  = sortrows(shifts_i_int, colOrder);
  % unique does sorting automatically
  [shifts_r_i, nInds_i]= deal(shifts_r_i(inds_s,:), nInds_i(inds_s));
  % -----shifts_i & NInd_c is now sorted-----
  [shifts_rc{ii}, inds_c{ii}(2:end)] = deal(single(shifts_r_i), nInds_i);
end
patch_c = [shifts_rc, inds_c];
end

