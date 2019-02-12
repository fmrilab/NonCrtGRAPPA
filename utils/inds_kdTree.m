function inds_c = inds_kdTree(loc0, loc1_kdtMdl, pSize, p)

if isa(loc1_kdtMdl, 'KDTreeSearcher')
  kdtMdl = loc1_kdtMdl;
else
  loc1_n = bsxfun(@rdivide, loc1_kdtMdl, (pSize-1)/2);
  kdtMdl = KDTreeSearcher(loc1_n, 'Distance', 'minkowski', 'P', p);
end
assert(strcmp(kdtMdl.Distance, 'minkowski'));
if p ~= kdtMdl.DistParameter, kdtMdl.DistParameter = p; end

% normalize sTraj coordinates wrt pSize, since matlab kd-tree doesn't support
% varying range along different direction.
loc0_n = bsxfun(@rdivide, loc0, (pSize-1)/2);
r = 1+1e-6; % just in case rangesearch use '<' instead of '<='

%% kdTree get indices
% [divide] & conqure
n0s = size(loc0_n,1 );
pool = gcp; % praise matlab parfor
nWorker = pool.NumWorkers;

n0s_dnc  = floor(n0s/nWorker); % _dnc for divide and conquer
n0s_xtra = mod(n0s, nWorker);
dnc_cnt  = [(n0s_dnc+1)*ones(n0s_xtra,1); n0s_dnc*ones(nWorker-n0s_xtra,1)];
loc0_n_dnc = mat2cell(loc0_n, dnc_cnt);

% divide & [conquer]
nInds_cc = cell(nWorker,1);
parfor ii = 1:nWorker
  nInds_cc{ii} = rangesearch(kdtMdl, loc0_n_dnc{ii}, r);
end % notice this is index for loc1
cInds_c = num2cell((1:n0s).');
nInds_c = cat(1, nInds_cc{:});

% inds_c = cellfun(@(c,n)uint64([c,n(:).']), cInds_c,nInds_c, 'Uni', 0);
inds_c = cellfun(@(c,n)[c,n(:).'], cInds_c,nInds_c, 'UniformOutput', 0);

end