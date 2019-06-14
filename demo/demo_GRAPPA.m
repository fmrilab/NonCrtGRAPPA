function demo_GRAPPA()
% The .mat testing files can be found at 
% https://drive.google.com/drive/folders/1FxF5jcMhL8Z2IzB-i__mb4Q3wDOmi20R
%
% (or one can generate ones own files based on the following comments)
% Here's some general comments about the .mat files:
%   phmA (struct):
%     .PS   (nx, ny, nz, nc), coil images
%   readout_info (struct):
%     .k (nk, nl) cycle/FOV/imSize, complex, normalized spiral trajectory
%     .imSize (2,), matrix size
%   sTraj (nk, 1+ndim) cycle/FOV, sampling trajectory, its 1st column is boolean
%     (0-unsampled, 1-sampled), following columns are coordinates in
%     [x, (y, (z))] direction(s).
%
%% load data
phmA = matfile('phantom_dat.mat');

%% uncomment lines to run specific demos
% do_Cartesian(phmA)
do_General(phmA)

end

function do_Cartesian(phmA)
%% choose a demo setting
pSize = [7,7]; % GRAPPA kernel size in cycle/FOV, do play with this, :)
islice = 5; % pick a slice for demonstration

%% prepare image and ACS
PS = phmA.PS;
PS = PS(:,:,islice,:);
imSize = size(PS);
imSize = imSize(1:3); % exclude #coil

kPS = mfft3(PS); % mfft3 also handles 2D fft properly
kcalib = kPS(86:115, 86:115, :,:); % ACS, center of k-space

%% create sampling mask
[fsMask, usMask] = deal(true(imSize));
usMask(:,1:2:end) = false; % acceleration factor R=2, skip every other readout
usMask(:, 86:115) = true;  % densely sampled center of k-space

%% synthesize a toy dataset
sysEqPara = sysEqPara_(fsMask, true(imSize));
Asyn = sysEqPara.F; % FFT system matrix
% for-loop vectorizeable, kept for consistency w/ do_general()
fkPS = Asyn*reshape(PS,[],size(PS,4));
ukPS = bsxfun(@times, usMask(:), fkPS);

PS_tilde = Asyn'*reshape(fkPS, [], size(PS,4));
aPS_tilde = Asyn'*reshape(ukPS, [], size(PS,4));

%% main, non-Cartesian GRAPPA recon
fP_direct = sosCombine(reshape(PS_tilde, size(PS)));
aP_direct = sosCombine(reshape(aPS_tilde, size(PS)));

% ukPS = reshape(ukPS, [imSize, size(PS,4)]);
uP_direct = test_direct(ukPS, usMask, kcalib, pSize, [], sysEqPara);

%% plots
figure,
subplot(221), imagesc(fP_direct); axis equal, title('fully sampled');
subplot(222), imagesc(aP_direct); axis equal, title('aliased');
subplot(223), imagesc(uP_direct); axis equal, title('directly recon''ed');

%%
% keyboard
end

function do_General(phmA)
%% choose a demo setting
spiral1_radial0 = 1;
do2D1_do3D0 = 1;
pSize = [7,7]; % GRAPPA kernel size in cycle/FOV, do play with this, :)
islice = 5; % pick a slice for demonstration

%% prepare image and ACS
PS = phmA.PS;
if do2D1_do3D0, PS = PS(:,:,islice,:); end
imSize = size(PS);
imSize = imSize(1:3); % exclude #coil

kPS = mfft3(PS); % mfft3 also handles 2D fft properly
kcalib = kPS(91:110, 91:110, :,:); % ACS, center of k-space

%% load a common sampling trajectory
if spiral1_radial0, readout_info = matfile('kSpiral2D.mat');
else,               readout_info = matfile('kRadial2D.mat');
end
[kxy, dcf] = getattrs(readout_info, {'k','w'});

%% create sampling trajectory
kmask = zeros(size(kxy)); % mask for under-sampling
kmask(:,1:2:end,:) = 1;   % acceleration factor R=2, skip every other readout

if do2D1_do3D0  % 2D example
  kTraj = bsxfun(@times, imSize(1:2), [real(kxy(:)), imag(kxy(:))]);
else            % 3D example
  islice = ':';
  str = input('3D demo, confirm your PC has >=32G RAM to proceed [N/y]', 's');
  if ~strcmpi(str, 'y'), error('Not enough RAM for 3D demo'); end

  pSize = [pSize, 3]; % GRAPPA kernel size in cycle/FOV, do play with this, :)
  nz = imSize(3);
  fn_popz = @(x)repmat(x, [1,1,nz]); % populate platters
  [kxy, dcf, kmask] = deal(fn_popz(kxy), fn_popz(dcf), fn_popz(kmask));

  kTraj = bsxfun(@times, imSize(1:2), [real(kxy(:)), imag(kxy(:))]);
  %% simulate 3D spiral/radial
  % The commented line simulates rotated stack-of-spiral/radials.
  % It is RAM demanding, run this test on a server with >32G RAM is recommended
  for iz = 1:nz, kmask(:,:,iz) = circshift(kmask(:,:,iz),[0,0,iz-1]); end

  kz = repmat(reshape((1:nz)-ctrInd(nz), 1,1,[]), [size(kxy,1),size(kxy,2),1]);
  kTraj = [kTraj, kz(:)];
end

[fsTraj, usTraj] = deal([ones(numel(kmask),1), kTraj], [kmask(:), kTraj]);

%% synthesize a toy dataset
sysEqPara = sysEqPara_(fsTraj, true(imSize), dcf);
Asyn = sysEqPara.F; % NUFFT system matrix
% mirt Gmri does not support well 3D coil image input, hence this loop
for ic = 1:size(PS,4), fkPS(:,ic) = Asyn*reshape(PS(:,:,:,ic),[],1); end
ukPS = bsxfun(@times, kmask(:)~=0, fkPS);
for ic = 1:size(PS,4)
  PS_tilde(:,ic) = Asyn'*bsxfun(@times, dcf(:), fkPS(:,ic));
  aPS_tilde(:,ic) = Asyn'*bsxfun(@times, dcf(:), ukPS(:,ic));
end

%% main, non-Cartesian GRAPPA recon
fP_direct = sosCombine(reshape(PS_tilde, size(PS)));
aP_direct = sosCombine(reshape(aPS_tilde, size(PS)));

uP_direct = test_direct(ukPS, usTraj, kcalib, pSize, imSize, sysEqPara);
uP_dogrid = test_dogrid(ukPS, usTraj, kcalib, pSize, imSize);

%% plots
if imSize(3) > 1  % display the picked slice
  [fP_direct, uP_direct, uP_dogrid] = ...
    deal(fP_direct(:,:,islice), uP_direct(:,:,islice), uP_dogrid(:,:,islice));
end

figure,
subplot(221), imagesc(fP_direct); axis equal, title('fully sampled');
subplot(222), imagesc(aP_direct); axis equal, title('aliased');
subplot(223), imagesc(uP_direct); axis equal, title('directly recon''ed');
subplot(224), imagesc(uP_dogrid); axis equal, title('recon''ed onto grid');

%%
% keyboard
end

%%
function uP = test_direct(ukPS, usMask_usTraj, kcalib, pSize, imSize ...
                          , sysEqPara, Tik)
% Directly reconstruct under-sampled k-space
%OUTPUTs:
% - uP (nx, ny, nz), combined final image
% - PS (nx, ny, nz, nc), reconstructed coil images
% - kPSG, reconstructed k-space

if ~exist('Tik', 'var'), Tik = 5e-7; end

if ~isempty(imSize)
  usTraj = usMask_usTraj;
  gMDL = grappaPrep(usTraj, kcalib, pSize, imSize, 'sysEqPara',sysEqPara ...
                    , 'Tik',Tik);
else
  usMask = usMask_usTraj;
  gMDL = grappaPrep(usMask, kcalib, pSize, 'sysEqPara',sysEqPara ...
                    , 'Tik',Tik);
end

[uP, PS, kPSG, ~] = grappa(ukPS, gMDL);

end

%%
function [uP] = test_dogrid(ukPS, usTraj, kcalib, pSize, imSize, Tik)
% Reconstruct a k-space grid
%OUTPUTs:
% - uP (nx, ny, nz), combined final image
% - PS (nx, ny, nz, nc), reconstructed coil images
% - kPSG, reconstructed k-space

if ~exist('Tik', 'var'), Tik = 5e-7; end

sampled = usTraj(:,1) ~= 0;
usTraj = usTraj(sampled, :);
ukPS   = ukPS(sampled,:);

gMDL = grappaPrep(usTraj, kcalib, pSize, imSize, 'doGrid',true, 'Tik',Tik);
[uP, PS, kPSG, ~] = grappa(ukPS, gMDL);

end

