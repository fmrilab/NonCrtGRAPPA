function [P, PS, kPSo, argo] = GrappaNaive(kPSi, gMDL, varargin)
%Inputs
% - kPS,
%    case sMask: (nx, ny, nz, nc, nf)
%    case sTraj:
%     kPSi (nk, nc, nf), unsampled data are set to 0. nf for n-frames
% OPTIONS

%%
argo = struct;
[argo, ~, argoXtra] = attrParser(argo, varargin); % for future usage
argo = mrgfield(argoXtra, argo);
%% Prep
patch_c = gMDL.patch_c;
npatch = size(patch_c,1);

%% form raw k-data vector, unify dimension into (nk, nc, nf)
if isattr(gMDL,'sMask')
  if 0
    kMask_d = any(any(kPSi(:,:,:,:,1), 4), 5); % _d: data; 4 nc, 5 nf
    % The following if-else is useful for compare rec quality across ACS
    % size, it allows gMDL.sMask to change while keeps kMask_d the same
    if ~isequal(kMask_d, kMask) %%% sMask MIS-MATCH %%%
      warning('data sMask mismatches calib sMask');
      [patch_cd] = getPatch(kMask_d, gMDL.pSize, gMDL.pDist, 'grid');
      [type, type_d] = deal(cat(patch_c{:,1},1), cat(patch_cd{:,1},1));
      
      [~, ind, ind_d] = intersect(type, type_d, 'rows', 'stable');
      patch_cd = patch_cd(ind_d); % reorder patch_cd to match calibed order
      
      % replace calib patch w/ reordered data patch, and coeff_c correspondinly
      [gMDL.coeff_c, patch_c] = deal(gMDL.coeff_c(ind), patch_cd);
    end
  end
  
  kMask = ~~gMDL.sMask;
  if isequal(size(kPSi(:,:,:,1)), size(gMDL.sMask))
    kPSi = reshape(kPSi, [], size(kPSi,4), size(kPSi,5));
  end
  nkDim = size(kMask);
  if ismatrix(kMask), nkDim = [nkDim, 1]; end
elseif isattr(gMDL,'sTraj')
  kMask = ~~gMDL.sTraj(:,1);
  nkDim = numel(kMask);
  kPSi = squeeze(kPSi); % ensure kPSi is (nk, nc, nf)
else,   error('unsupported');
end

[~, nc, nf] = size(kPSi);
[nk1, nk] = deal(nnz(kMask), numel(kMask));
nk0 = nk - nk1;
if size(kPSi, 1) ~= nk1, kPSi = kPSi(kMask, :,:); end % now kPSi is kPS1

%% recon
if npatch
  kPS0_c = cell(npatch, 1);
  if ~isa(gMDL, 'matlab.io.MatFile')
    kPS0_c = doRecon(patch_c(:,2), kPSi, gMDL.coeff_c);
  else % gMDL is a matfile handle, assume the coeffs set is too large to load
    [nP, ndnc] = deal(size(patch_c,1), 10); % 10 was chosen arbitrarily
    ie = min(nP, floor(nP/ndnc*(1:ndnc)));
    ib = [1, ie(1:end-1)+1];
    
    for ii = 1:numel(ib)
      coeff_c = gMDL.coeff_c(ib(ii):ie(ii),1);
      kPS0_c(ib(ii):ie(ii)) = doRecon(patch_c(ib(ii):ie(ii),2), kPSi,coeff_c);
    end
  end
  
  % embedding or reshaping and k -> img
  cInd = cell2mat( cellfun(@(x)x(:,1), patch_c(:,2), 'UniformOutput',false) );
  kPS0 = zeros(nk0,size(kPSi,2), size(kPSi,3));
%   kPS0 = cat(1, kPS0_c{:});
  kPS0(cInd,:,:) = cat(1, kPS0_c{:}); % restore sTraj order
else
  kPS0 = 0;
end

if gMDL.doGrid
  kMask = mloc2mask([~gMDL.sTraj(:,1), gMDL.sTraj(:,2:end)], gMDL.imSize);
  kPS0 = reshape(kPS0, [size(kPS0,1), ones(1,3-ndims(kMask)), nc, nf]);
  kPSo = embed(kPS0, kMask);
else
  kPSo = zeros([prod(nkDim), nc, nf]);
  kPSo(~~kMask, :,:) = kPSi;
  kPSo( ~kMask, :,:) = kPS0;
  kPSo = reshape(kPSo, [nkDim, nc, nf]);
end

k2imgPara = gMDL.k2imgPara;

PS = k2imgPara.kPS2PS_fn(kPSo);
P  = k2imgPara.PS2P_fn(PS);

end

%% sub-function that does the recon
function kPS0_c = doRecon(inds_c, kPS1, coeff_c)
%Inputs:
% - kPS1 (nk, nc, nf)
%Outputs:
% - kPS0 (numel(patch_c),): cells of (nCtr, nc, nf)
[~, nc, nf] = size(kPS1);
nType = numel(inds_c);
kPS0_c = cell(nType,1);

parfor ii = 1 : nType % parfor
  [ctrInds, ngbInds] = deal(inds_c{ii}(:,1), inds_c{ii}(:,2:end));
  if isempty(ngbInds)
    kPS0_i = zeros(numel(ctrInds), nc, nf); % when no ngb was captured
  else
    kNgbs = permute(kPS1(ngbInds, :,:), [3,1,2]);  % -> (nf, nNgb*nCtr, nc)
    kNgbs = reshape(kNgbs, nf*numel(ctrInds), []); % -> (nCtr*nf, nc*nNgb)
    
    kPS0_i = permute(reshape((kNgbs*coeff_c{ii}).', nc,nf,[]), [3,1,2]);
  end
  kPS0_c{ii} = kPS0_i;
end
end

