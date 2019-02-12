function gMDL = grappaPrep(sMask_sTraj, kcalib, pSize, varargin)
% gMDL = grappaPrep(sMask, kcalib, pSize, varargin)
%or
% gMDL = grappaPrep(sTraj, kcalib, pSize, imSize, varargin)
%Input
% - sMask_sTraj
%    sMask (nx, ny), integer, 0 for not sampled
%    sTraj (nk, ndim+1), cycle/fov, straj(:,1) is a mask for underSampling:
%      0 for undersampled; >0 for sampled and will be used in recon; <0 for
%      sampled but will NOT be used in recon.
% - kcalib (nkx, nky, nkz, nc), ACS on a cycle/FOV Cartesian grid
% - pSize  (2,), size of grappa recon patch
% - imSize (2,), size of the image, required input only for sTraj
%OPTIONAL
% - pDist (1,), used as p-norm, describing patch shape as norm ball
% - sysEqPara struct, describe system eq parameters, see also sysEqPara_()
% - doGrid (t/F), for sTraj, if true, directly reconstruct grid k-space
% (.calibPara wraps below into output)
% - doSift (T/f), for sTraj, if true, sifting too crowded patch neighbors
% - Tik (1,) [0.0], Tikhonov regularizor lambda
% - calibSize (nd, ): resize kcalib for ram safety if input size is too large
% (.k2imgPara wraps below into output)
% - kPS2PS_fn fn_handle, reconstruct coil images from full coil data
% - PS2P_fn   fn_handle, combine coil images into a single one
%Output
% - gMDL struct, besides described in Input/OPTIONAL section, adding fields:
%   .imSize (nx, ny, nz); size of recon'ed image
%   .pSize  (px, py, pz); size of recon patch
%   .pMask  (px, py, pz); patch ball inclusion mask expressed by Cartesian grid
%   .patch_c {numel(coeff),2}, storing the ngb indices used in recon
%   .coeff_c {#distinct rec pat, 1}, coeff of GRAPPA kSpace interpolation

%% prep
%%% gMDL pairing

kcalib = single(kcalib);
kcalib = kcalib/mean(sqrt(sum(abs(reshape(kcalib,[],size(kcalib,4))).^2,1)),2);

gMDL.pDist = inf;
gMDL.sysEqPara = [];

if isinteger(sMask_sTraj) || islogical(sMask_sTraj)
  imSize = size(sMask_sTraj);
  
  [gMDL, extra] = attrParser(gMDL, varargin);
  gMDL.sMask = sMask_sTraj;
  
  isFullySampled = ~nnz(~gMDL.sMask);
elseif isfloat(sMask_sTraj)
  if max(max(abs(sMask_sTraj(:,2:end)))) < 0.6, error('sTraj cycle/FOV'); end
  imSize = varargin{1};
  
  gMDL.doGrid = false;
  [gMDL, extra] = attrParser(gMDL, varargin(2:end));
  gMDL.sTraj = sMask_sTraj;
  
  isFullySampled = ~nnz(~gMDL.sTraj(:,1));
else, error('not supported');
end

[opt.doSave, opt.doFFTTrick] = deal(false, true); % true fftTrick for circ bndry
[opt, extra] = attrParser(opt, extra);

if isempty(gMDL.sysEqPara)
  kfull = sMask_sTraj;
  if isattr(gMDL, 'sMask'), kfull = true(size(kfull));
  else, kfull(:,1) = 1;
  end
  gMDL.sysEqPara = sysEqPara_(kfull,true(imSize));
end
[gMDL.imSize, gMDL.pSize] = deal(imSize, pSize);

gMDL.pMask = normBallMask(pSize, gMDL.pDist);

if isattr(gMDL, 'sTraj') && gMDL.doGrid
  % NOTICE: this assumes fed-in sTraj are fully sampled
  crtMask = imdilate(mloc2mask(gMDL.sTraj, imSize), gMDL.pMask);
  sTraj0 = mask2mloc(crtMask); % add in grid locations into sTraj
  gMDL.sTraj = [gMDL.sTraj; [zeros(nnz(crtMask),1), sTraj0(crtMask, 2:end)]];
end

%%% calibPara pairing
calibSize = size(kcalib); % [sx, sy, sz, nc]
calibPara.Tik = 5e-7;
calibPara.calibSize = calibSize;

if isattr(gMDL,'sTraj')
  [calibPara.doSift, calibPara.doExtend] = deal(true);
end
[calibPara, extra] = attrParser(calibPara, extra);

f = {'Tik','calibSize','doSift','doExtend'};
calibPara = chkfield(calibPara, f); % completing missing fields with []

%%% k2imgPara pairing
[k2imgPara.kPS2PS_fn, k2imgPara.PS2P_fn] = deal([]);
[k2imgPara, extra] = attrParser(k2imgPara, extra);
if isempty(k2imgPara.kPS2PS_fn)
  k2imgPara.kPS2PS_fn = @(kPS)fPSrec(kPS, gMDL.sysEqPara);
end
if isempty(k2imgPara.PS2P_fn)
  k2imgPara.PS2P_fn = @(PS)sosCombine(PS);
end

f = {'kPS2PS_fn','PS2P_fn'};
k2imgPara = chkfield(k2imgPara, f);

%%% add into gMDL
gMDL.calibPara = calibPara;
gMDL.k2imgPara = k2imgPara;
if ~isempty(extra), warning('unpairable optional inputs were found'); end

%% add up .patch_c and .coeff_c
if isFullySampled && (isattr(gMDL, 'sMask') || ~gMDL.doGrid)
  [gMDL.patch_c, gMDL.coeff_c] = deal({});
  if opt.doSave, gMDL = saveMDL(gMDL, 'gMDL'); end
else
  %% identify patches
  if isattr(gMDL,'sMask')
    patch_cr = getPatch(gMDL.sMask, pSize, gMDL.pDist, 'grid'); % _r: raw
  else % sTraj
    patch_cr = getPatch(gMDL.sTraj, pSize, gMDL.pDist, 'kdtree');
  end
  
  patch_c = grpPatch(patch_cr, 'doSift',calibPara.doSift);
  
  shifts_c = patch_c(:,1);
  gMDL.patch_c = patch_c;
  % doSave, do not move this line before calculating patch_c
  if opt.doSave, gMDL = saveMDL(gMDL, 'gMDL'); end
  %% do calibration
  %%% this block generates gMDL.coeff
  if ~isequal(size(kcalib), calibPara.calibSize)
    kcalib = mfft3(resize(mifft3(kcalib), calibPara.calibSize));
  end
  % kcalib = resize(kcalib, calibPara.calibSize); % should resize in img-dom
  nc = size(kcalib, 4);
  
  if isattr(gMDL, 'sMask') % don't change order
    % sMask calib
    gMDL.coeff_c = pinvCalib_sMask(shifts_c, kcalib, calibPara.Tik, gMDL.pMask);
  elseif isattr(gMDL, 'sTraj')
    % sTraj calib
    calib = mifft3(kcalib); % kcalib -> calib, for phs-shift convenience
    
    miscSt.doFFTTrick = opt.doFFTTrick;
    
    % Croping valid region for calib does not make large difference in recon
    if miscSt.doFFTTrick
      prep_fn = @(x)x;
      nComb = nc*(nc+1)/2;
      calib = single(calib);
      [miscSt.F_c, miscSt.dgrid] = FFTTrick(calib, nComb);
    else
      crop_fn = getCrop(calibSize,pSize);
      prep_fn = @(x)crop_fn(mfft3(x));
    end
    miscSt.prep_fn = prep_fn;
    
    if ~opt.doSave
      gMDL.coeff_c = pinvCalib_sTraj(shifts_c, calib, calibPara.Tik, miscSt);
    else % gMDL is a matfile handle, assume the coeffs set is too large to load
      [nP, ndnc] = deal(size(patch_c,1), 10); % 10 was chosen arbitrarily
      ie = min(nP, floor(nP/ndnc*(1:ndnc)));
      ib = [1, ie(1:end-1)+1];
      
      for ii = 1:numel(ib)
        disp(ii);
        gMDL.coeff_c(ib(ii):ie(ii),1) = pinvCalib_sTraj(...
          shifts_c(ib(ii):ie(ii)), calib, calibPara.Tik, miscSt);
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
gMDL = chkfield(gMDL, f, f_ignore);
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

NhC = AhA(:,icAllc); % [ nc*pSize, nc ], nc: # of coils
NhN = AhA;
NhN(:, icAllc) = [];

type = cat(2, shifts_c{:}).';

nType = size(type, 1);
coeff_c = cell(nType, 1);

parfor iType = 1:nType
  % selecting out needed elements in AtA and Atb
  typeM  = repmat(type(iType,:),[1, nc]);
  NhNtmp = NhN(typeM(:), typeM(:)); % [nc*nNgb, nc*nNgb]
  NhCtmp = NhC(typeM(:), :);        % [nc*nNgb, nc]
  if Tik
    nrow = size(NhNtmp,1);
    nlam = Tik/nrow; % normalize by matrix norm
    ind = 1:(nrow+1):nrow^2;
    NhNtmp(ind) = NhNtmp(ind) + nlam;
  end
  coeff_c{iType} = pinv(NhNtmp)*NhCtmp;
end
end

%% sub-function that does the sTraj calib
% This is gonna take time
function [coeff_c] = pinvCalib_sTraj(shifts_c, calib, Tik, miscSt)
% patch_c just to unify interface with cgCalib
% for GRAPPA, this method nolonger transforms img-domain phase-shifted ACS data
% back to the frequency domain and trim, cuz the outputs are empirically similar
nc = size(calib,4);
prep_fn = miscSt.prep_fn;
Ctr = reshape(prep_fn(calib), [], nc);
ns = numel(shifts_c);
if ~ns, disp('no calib needed?'); coeff_c = {}; return; end

%% calibration
coeff_c = cell(ns,1);
if ~miscSt.doFFTTrick
  calib = permute(calib, [1,2,3,5,4]); % -> (nx, ny, nz, 1, nc)
  parfor ii = 1:ns
    sh_i = shifts_c{ii};                                 % t1 = cputime;
    coeff_c{ii} = formLS2(sh_i,Tik,nc,calib,Ctr,prep_fn); % t2 = cputime;
    % disp(num2str([ii, size(NhN,1), t2-t1]));
    % disp('---');
  end
else % FFT trick
  [F_c, dgrid] = deal(miscSt.F_c, miscSt.dgrid);
  
  idgrid = single(1./dgrid);
  parfor ii = 1:ns
    sh_i = shifts_c{ii};                                          % t1=cputime;
    coeff_c{ii} = LS_fft_mex(sh_i', nc, Tik, F_c, idgrid);
    % coeff_c{ii}=formLS1(sh_i,Tik,nc,F_c,dgrid,iofst,m2,coilSub); % t2=cputime;
    % disp(num2str([ii, size(NhN,1), t2-t1]));
    % disp('---');
  end
end
end

%%
function coeff = formLS2(sh, Tik, nc, calib, Ctr, prep_fn)
% ugly implemented function interface, but may save time
% t = cputime;
calibSize = size(calib);
[~, phs] = phs4Shift(sh, 'dims',calibSize(1:3), 'doFFTShift',true);
calib_iter = bsxfun(@times, calib, phs); % (nx,ny,nz,ns,nc);

% disp(cputime - t);
nNgb = size(sh,1);
% Ngb = reshape(prep_fn(calib_iter), [], nc, nNgb);
% Ngb = reshape(permute(Ngb, [1,3,2]), [], nc*nNgb); % -> (nk, nc*nNgb)
Ngb = reshape(prep_fn(calib_iter), [], nc*nNgb); % -> (nk, nc*nNgb);


%% accelerating by resizing
[nrow, ncol] = size(Ngb);
Ngb = imresize(Ngb, [min(floor(2*ncol), nrow), ncol], 'box');
Ctr = imresize(Ctr, [min(floor(2*ncol), nrow), size(Ctr,2)], 'box');
%%
% t = cputime;
NhN = Ngb'*Ngb; % url:goo.gl/Chdjh6, matlab already optimized AtA calc
NhC = Ngb'*Ctr;
% disp(cputime-t);
% t = cputime;
if Tik
  nNgbnC = nNgb*nc;
  nNgbTik = nNgb*Tik; % normalize by matrix norm
  inds = 1:(nNgbnC+1):nNgbnC*nNgbnC;
  NhN(inds) = NhN(inds) + nNgbTik; % adding Tik reg (scaled identity matrix)
end

% pinv is more robust than '\'; pinv(NhN)*NhC is faster than pinv(Ngb)*Ctr
% While w/ small matrix size Ngb\Ctr is faster, profiling (2016a) shows that
% when matrix size gets large (e.g. 3d recon 5x5x5), pinv wins.
% [NhN, NhC] = deal(single(NhN), single(NhC));
coeff = pinv(NhN)*NhC;
% disp(cputime - t);
end

%%
function [coeff] = formLS1(sh, Tik, nc, F_c, dgrid, iofst, m2, coilSub)
% a MATLAB based FFTtrick implementation, not super fast, use mex
[nNgb, ndim] = size(sh);

fn = @(x)[reshape(bsxfun(@minus, x.', x), [], 1); -x; x];

sh = bsxfun(@rdivide, sh, dgrid); % normalized shifts by grid units

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
function [F_c, dgrid] = FFTTrick(calib, nComb)
% F_c are cells of Hadamard product between conj(calib) and calib
[nd, calibh] = deal(size(calib), conj(calib));
if nd(3)==1,[nc,nd,pd]=deal(nd(4),nd(1:2),[128,128]);
else, [nc,nd,pd]=deal(nd(4),nd(1:3),[128,128,nd(3)]);
end % pd set empirically

[ich, ic] = ndgrid(1:nc);
ltMask = ich>=ic; % lower triangle mask
coilSub = [ich(ltMask), ic(ltMask)];

dgrid = (nd-1)./(pd-1); % unit distance of the padded Cartesian grid

F_c = cell(1, nComb);
calib0 = zeros(pd, class(calib));

ncSub = ctrSub(nd);
fn1 = @(ii)mod((1:nd(ii))-ncSub(ii), pd(ii))+1; % inds after proper 0-padding
s_c = cellfun(fn1, num2cell(1:numel(nd)), 'Uni', false); % cell of subscripts

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
