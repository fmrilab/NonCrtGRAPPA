function gMDL = grappaPrep(sMask_sTraj, kcalib, pSize, varargin)
%function gMDL = grappaPrep(sMask_sTraj, kcalib, pSize, varargin)
%
% This function contains non-Cartesian GRAPPA calibration algorithm, see,
% doi: [10.1002/mrm.27801](https://doi.org/10.1002/mrm.27801)
% The function finishes calibration and prepares model, `gMDL`, to be used in
% reconstruction, `grappa.m`.
%
%GENERAL COMMENTS
% In this file, patch and constellation are synonyms.
% `_c` as a prefix indicates a variable being a cell array.
% Let np be the number of distinct constellations/patches.
%
%Input
% - sMask_sTraj, sampling info
%   sMask (nx, (ny, (nz,))), Cartesian readout.
%     integer or logical: 0, not sampled; 1, sampled.
%   sTraj (nk, 1+ndim), General readout, each row is a k-space location
%     straj(:,1): 0, not sampled; 1, sampled.
%     sTraj(:,2:end): k-space coordinate, cycle/fov,
% - kcalib (nkx, nky, nkz, nc), ACS on a cycle/FOV Cartesian grid
% - pSize  (ndim,), cycle/FOV, patch size along each dimension
% - imSize (ndim,), image size, required input only for sTraj
%OPTIONAL
% - pDist (1,), as p-norm ball, demarcates patch shape, i.e., a sampled point
%   that locates w/in a p-norm ball centered at a ctr can be used as its ngb.
% - sysEqPara struct, describe system eq parameters, see also sysEqPara_()
% - doGrid     [t/F], for sTraj, if true, directly reconstruct a grid k-space.
% - doSave     [t/F], save gMDL as a matfile, and return its handle.
% (.calibPara wraps below into output)
% - doSift     [T/f], for sTraj, if true, sifts too crowded patch ngbs
% - Tik (1,), Tikhonov regularizor lambda
% - calibSize (ndim,): resize kcalib for ram safety when input size is too large
% (.k2imgPara wraps below into output)
% - kPS2PS_fn fn_handle, reconstruct coil images from full coil data
% - PS2P_fn   fn_handle, combine coil images into a single one
%Output
% - gMDL struct, besides described in Input/OPTIONAL section, adding fields:
%   .imSize (ndim,); size of recon'ed image
%   .pSize  (ndim,); size of recon patch
%   .pMask  (px, py, pz); patch ball inclusion mask expressed by Cartesian grid
%   .patch_c (np, 2) cell, record ngb indices (in sTraj) of each patch
%   .coeff_c (np, 1) cell, GRAPPA weights of each patch

%% prep
%%% gMDL parsing

% ACS normalization, for consistent Tik behavior.
kcalib = single(kcalib);
kcalib = kcalib/mean(sqrt(sum(abs(reshape(kcalib,[],size(kcalib,4))).^2,1)),2);

gMDL.pDist = inf; % inf-norm ball sets constellation region to be a square/cube.
gMDL.sysEqPara = [];

if isinteger(sMask_sTraj) || islogical(sMask_sTraj) % Cartesian Readout
  imSize = size(sMask_sTraj);
  
  [gMDL, extra] = attrParser(gMDL, varargin);
  gMDL.sMask = sMask_sTraj;
  
  isFullySampled = ~nnz(~gMDL.sMask);
elseif isfloat(sMask_sTraj)                         % General Readout.
  if max(max(abs(sMask_sTraj(:,2:end)))) < 0.6, error('sTraj cycle/FOV'); end
  imSize = varargin{1};
  
  gMDL.doGrid = false;
  [gMDL, extra] = attrParser(gMDL, varargin(2:end));
  gMDL.sTraj = sMask_sTraj;
  
  isFullySampled = ~nnz(~gMDL.sTraj(:,1));
else, error('not supported');
end

opt.doSave = false;
[opt, extra] = attrParser(opt, extra);

gMDL.pMask = normBallMask(pSize, gMDL.pDist);

if isattr(gMDL, 'sTraj') && gMDL.doGrid
  assert(all(gMDL.sTraj(:,1)), 'For doGrid: input only sampled locations')
  % Determine grid loc's that can possibly be recon'ed, append to the sTraj.
  crtMask = imdilate(mloc2mask(gMDL.sTraj, imSize), gMDL.pMask);
  sTraj0 = mask2mloc(crtMask); % add in grid locations into sTraj
  gMDL.sTraj = [gMDL.sTraj; [zeros(nnz(crtMask),1), sTraj0(crtMask, 2:end)]];
end

if isempty(gMDL.sysEqPara)
  kfull = sMask_sTraj;
  if isattr(gMDL, 'sMask') || gMDL.doGrid, kfull = true(imSize);
  else, kfull(:,1) = 1;
  end
  gMDL.sysEqPara = sysEqPara_(kfull,true(imSize));
end
[gMDL.imSize, gMDL.pSize] = deal(imSize, pSize);

%%% calibPara parsing
calibSize = size(kcalib); % [cx, cy, cz, nc]
calibPara.Tik = 5e-7;
calibPara.calibSize = calibSize;
if isattr(gMDL,'sTraj'), calibPara.doSift = true; end

[calibPara, extra] = attrParser(calibPara, extra);

f = {'Tik','calibSize','doSift'};
calibPara = chkattrs(calibPara, f); % completing missing fields with []

%%% k2imgPara parsing
[k2imgPara.kPS2PS_fn, k2imgPara.PS2P_fn] = deal([]);
[k2imgPara, extra] = attrParser(k2imgPara, extra);
if isempty(k2imgPara.kPS2PS_fn)
  k2imgPara.kPS2PS_fn = @(kPS)fPSrec(kPS, gMDL.sysEqPara);
end
if isempty(k2imgPara.PS2P_fn)
  k2imgPara.PS2P_fn = @(PS)sosCombine(PS);
end

f = {'kPS2PS_fn','PS2P_fn'};
k2imgPara = chkattrs(k2imgPara, f);

%%% add into gMDL
gMDL.calibPara = calibPara;
gMDL.k2imgPara = k2imgPara;
if ~isempty(extra), warning('un-parsable optionals inputs were found'); end

%% Calibration
if isFullySampled && (isattr(gMDL, 'sMask') || ~gMDL.doGrid)
  [gMDL.patch_c, gMDL.coeff_c] = deal({});
  if opt.doSave, gMDL = saveMDL(gMDL, 'gMDL'); end
else
  % identify patches/constellations
  if isattr(gMDL,'sMask')
    patch_cr = getPatch(gMDL.sMask, pSize, gMDL.pDist, 'grid'); % _r: raw
  else % sTraj
    patch_cr = getPatch(gMDL.sTraj, pSize, gMDL.pDist, 'kdtree');
  end
  
  patch_c = grpPatch(patch_cr, 'doSift',calibPara.doSift);
  
  shifts_c = patch_c(:,1); % relative shifts from ctr to ngbs of all patches
  gMDL.patch_c = patch_c;
  % doSave, do not move this line to before computing patch_c
  if opt.doSave, gMDL = saveMDL(gMDL, 'gMDL'); end
  
  %%% this block calibrates gMDL.coeff
  if ~isequal(size(kcalib), calibPara.calibSize)
    kcalib = mfft3(resize(mifft3(kcalib), calibPara.calibSize));
  end
  nc = size(kcalib, 4);
  
  if isattr(gMDL, 'sMask') % don't change order
    % sMask calib
    gMDL.coeff_c = pinvCalib_sMask(shifts_c, kcalib, calibPara.Tik, gMDL.pMask);
  elseif isattr(gMDL, 'sTraj')
    % sTraj calib
    calib = mifft3(kcalib); % kcalib -> calib, for phs-shift convenience
    
    % FFTTrick utilizes the circulant boundary described in the paper
    nComb = nc*(nc+1)/2;
    calib = single(calib);
    [miscSt.F_c, miscSt.dgrid] = FFTTrick(calib, nComb);
    
    if ~opt.doSave
      gMDL.coeff_c = pinvCalib_sTraj(shifts_c, calib,calibPara.Tik,miscSt);
    else % gMDL is a matfile handle, assume the coeffs set is too large to load
      [nP, ndnc] = deal(size(patch_c,1), 10); % 10 was chosen arbitrarily
      ie = min(nP, floor(nP/ndnc*(1:ndnc)));
      ib = [1, ie(1:end-1)+1];
      
      for ii = 1:numel(ib)
        disp(ii);
        gMDL.coeff_c(ib(ii):ie(ii),1) = ...
          pinvCalib_sTraj(shifts_c(ib(ii):ie(ii)), calib,calibPara.Tik,miscSt);
      end
    end
    parfevalOnAll(@clear, 0, 'mex');
  end
end
%% complete fields and return
f = {'sysEqPara','pDist','doGrid','imSize','pSize','pMask'...
  ,'calibPara','k2imgPara','patch_c','coeff_c'};
% 'Properties' field may appear when doSave
f_ignore = {'sMask', 'sTraj', 'Properties'}; % not to be completed
gMDL = chkattrs(gMDL, f, f_ignore);
end

%% sub-function that does the sMask calib
function [coeff_c] = pinvCalib_sMask(shifts_c, kcalib, Tik, pMask)
pSize = size(pMask);
pMask_float = double(pMask);
pMask_float(ctrInd(pSize)) = inf;
pMask_float = pMask_float(pMask);
pcInd = find(isinf(pMask_float));

nc = size(kcalib, 4);

bigA = reshape(patch2row(kcalib, pMask, 'pad0'), [], nc,nnz(pMask));
bigA = reshape(permute(bigA, [1,3,2]), [], nc*nnz(pMask));
AhA  = bigA'*bigA;
icAllc = pcInd + (0:nc-1)*nnz(pMask);
AhA(icAllc, :) = [];

NhC = AhA(:,icAllc); % (nc*pSize, nc), nc: # of coils
NhN = AhA;
NhN(:, icAllc) = [];

type = cat(2, shifts_c{:}).';

nType = size(type, 1);
coeff_c = cell(nType, 1);

parfor iType = 1:nType
  % selecting out needed elements in AtA and Atb
  typeM  = repmat(type(iType,:),[1, nc]);
  NhNtmp = NhN(typeM(:), typeM(:)); % (nc*nNgb, nc*nNgb)
  NhCtmp = NhC(typeM(:), :);        % (nc*nNgb, nc)
  if Tik
    nrow = size(NhNtmp,1);
    nlam = Tik/nrow; % normalize by matrix norm
    ind = 1:(nrow+1):nrow^2;
    NhNtmp(ind) = NhNtmp(ind) + nlam;
  end
  coeff_c{iType} = pinv(NhNtmp)*NhCtmp;
end
end

%%
% This is gonna take time
function [coeff_c] = pinvCalib_sTraj(shifts_c, calib, Tik, miscSt)
% sub-function for sTraj calib
%INPUTS:
% - shifts_c (np,) cell, each element corresponds to a distinct patch, and
%   contains relative shifts from ctr to ngbs of that patch
% - calib (cx, cy, cz, nc), image domain ACS, i.e. low-res coil images
% - Tik (1,), Tikhonov regularizor lambda
% - miscSt structure,
%   .F_c (nComb,) cell, Hadamard product between conj(calib) and calib
%   .dgrid (ndim,), unit distance in the padded Cartesian grid

nc = size(calib,4);
Ctr = reshape(calib, [], nc);
ns = numel(shifts_c);
if ~ns, disp('no calib needed?'); coeff_c = {}; return; end

%% calibration
coeff_c = cell(ns,1);
[F_c, dgrid] = deal(miscSt.F_c, miscSt.dgrid);

idgrid = single(1./dgrid); % i for inverse
parfor ii = 1:ns
  sh_i = shifts_c{ii};                                          % t1=cputime;
  coeff_c{ii} = LS_fft_mex(sh_i', nc, Tik, F_c, idgrid);
  % coeff_c{ii} = LS_fft_mat(sh_i, nc, Tik, F_c, idgrid); % t2=cputime;
  % disp(num2str([ii, size(NhN,1), t2-t1]));
  % disp('---');
end
end

%% Least Square calibration with the FFT trick
function [coeff] = LS_fft_mat(sh, nc, Tik, F_c, idgrid)
% This function essencially does the samething as LS_fft_mex().
%INPUTS:
% - sh (nngb, ndim), relative shifts from ctr to ngbs
% - nc (1,) #coils
% - Tik (1,), Tikhonov regularizor lambda
% - F_c (nComb,) cell, Hadamard product between conj(calib) and calib
% - idgrid (ndim,), inverse unit distance in the padded Cartesian grid

ndim = size(sh,2);
m2 = (0:2^(ndim)-1)'; % m2: binary digits mask
m2 = bsxfun(@gt, bsxfun(@mod, m2, 2.^(1:ndim)), 2.^(0:ndim-1)-1);

pd = size(F_c{1}); % size of dense spectrum
iofst = sum(bsxfun(@times, m2, cumprod([1, pd(1:end-1)])), 2);
iofst = iofst(:);

% a MATLAB based FFTtrick implementation, not super fast, use mex
[nNgb, ndim] = size(sh);
coilSub = getCoilSub(nc);

fn = @(x)[reshape(bsxfun(@minus, x.', x), [], 1); -x; x];

sh = bsxfun(@times, sh, idgrid); % normalized shifts by grid units

%% indices & coefficients for linear interpolation
coef = ones(2^ndim, nNgb*(nNgb+2), 'single'); % coeffs for linear interpolation

pd = size(F_c{1});
ofstMult = cumprod([1, pd(1:end-1)]);
ofstBase = 0;
for ii = 1:ndim
  sh_tmp = fn(sh(:,ii));
  ofstBase = ofstBase + ofstMult(ii)*floor(sh_tmp);
  coef = coef .* abs(bsxfun(@minus, 1-m2(:,ii), mod(sh_tmp,1)'));
end

cInd = ctrInd(pd);
inds = bsxfun(@plus, cInd+iofst(:), ofstBase(:)'); % inds for lin-interpolation

%%
nNgbnC = nNgb*nc;
[NhN, NhC] = deal(zeros(nNgbnC, 'single'), zeros(nNgbnC,nc, 'single'));
for ii = 1:size(coilSub,1)
  [i1_i, i2_i] = deal(coilSub(ii,1), coilSub(ii,2));
  [b1,e1,b2,e2] = deal((i1_i-1)*nNgb+1, i1_i*nNgb, (i2_i-1)*nNgb+1, i2_i*nNgb);
  tmp1 = sum(single(reshape(F_c{ii}(inds), 2^ndim, [])) .* coef, 1);
  tmp1 = reshape(tmp1, nNgb, nNgb+2);
  [tmp2, tmp3] = deal(tmp1(:,end-1), tmp1(:,end));
  tmp1(:,end-1:end) = [];
  [NhN(b1:e1, b2:e2), NhC(b1:e1, i2_i)] = deal(tmp1, tmp2);
  if i1_i ~= i2_i
    NhN(b2:e2, b1:e1) = tmp1';
    NhC(b2:e2, i1_i)  = conj(tmp3);
  end
end

if Tik
  nNgbnC = nNgb*nc;
  nlam = Tik*nNgb; % normalize by matrix norm
  dinds = 1:(nNgbnC+1):nNgbnC*nNgbnC;  % indices to diagonal comp
  NhN(dinds) = NhN(dinds) + nlam; % adding Tik reg (scaled identity matrix)
end
% [NhN, NhC] = deal(single(NhN), single(NhC));
coeff = pinv(NhN)*NhC;
end

%%
function coilSub = getCoilSub(nc)
% For indexing coil combinations in AhA.
[ich, ic] = ndgrid(1:nc);
ltMask = ich>=ic; % lower triangle mask
coilSub = [ich(ltMask), ic(ltMask)];
end

%%
function [F_c, dgrid] = FFTTrick(calib, nComb)
% F_c are cells of Hadamard product between conj(calib) and calib
[cSize, calibh] = deal(size(calib), conj(calib));
if cSize(3)==1,[nc,cSize,pd]=deal(cSize(4),cSize(1:2),[128,128]);
else, [nc,cSize,pd]=deal(cSize(4),cSize(1:3),[128,128,cSize(3)]);
end % pd set empirically

coilSub = getCoilSub(nc);

dgrid = (cSize-1)./(pd-1); % unit distance of the padded Cartesian grid

F_c = cell(1, nComb);
calib0 = zeros(pd, class(calib));

ncSub = ctrSub(cSize);
fn1 = @(ii)mod((1:cSize(ii))-ncSub(ii), pd(ii))+1; % inds after proper 0-padding
s_c = cellfun(fn1, num2cell(1:numel(cSize)), 'Uni', false); % cell of subscripts

kshift = ctrSub(pd)-1;
for ii = 1:nComb % huge FFT is too slow in parfor
  chc_i = calibh(:,:,:,coilSub(ii,1)).*calib(:,:,:,coilSub(ii,2));
  % trick: replacing 0-field is faster than proper 0-padding
  calib0(s_c{:}) = chc_i;
  % One may be remove the circshift below, then query w/ modded locs, while it
  % may shorten the time of this for-loop, the moding operation may aggravate
  % the cost of the enormous query operations; Also, the modded query only works
  % for interpolation methods that strictly do not involve out-of-bndry values.
  
  F_c{ii} = circshift(fft3(calib0), kshift);
end

end
