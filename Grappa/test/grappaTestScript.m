function grappaTestScript()
% The .mat testing files is not in the repo due to their size.
% Instead, they can be found at www-personal.umich.edu/~tianrluo/Phantoms/
% (or one can generate ones own files based on the following comments)
% Here's some general comments about the .mat files:
%   phmA (struct):
%     .kPS  (nx, ny, nz, nc), kspace Phantom w/ SMaps (indexed by nc(oil))
%     .PS   (nx, ny, nz, nc), Phantom w/ SMaps
%     .sMap (nx, ny, nz, nc), sensitivity maps
%   readout_info (struct):
%     .k (nk, nl) cycle/FOV/imSize, complex, normalized spiral trajectory
%     .imSize (2,), matrix size
%   mA (struct):
%     .imask_unif (nx, ny), boolean, a uniformly under-sampling mask
%   sTraj (nk, 1+ndim) cycle/FOV, sampling trajectory, its 1st column is boolean
%     (0-unsampled, 1-sampled), following columns are coordinates in
%     [x, (y, (z))] direction(s).
%

%% tests
% uncomment lines to run specific tests

phmA = matfile('phantom_dat.mat');
[PS, kPS] = getattrs(phmA, {'PS', 'kPS'});
kcalib = kPS(91:110, 91:110, :,:); % ACS, center of k-space

spiral1_radial0 = 0;
if spiral1_radial0, readout_info = matfile('kSpiral2D.mat');
else,               readout_info = matfile('kRadial2D.mat');
end

[fov, imSize, k, dcf] = getattrs(readout_info, {'fov','imSize','k','w'});
kTraj = bsxfun(@times, imSize, [real(k(:)), imag(k(:))]);

kmask = zeros(size(k));
kmask(:,1:2:end) = 1;  % config under-sampling

[fsTraj, usTraj] = deal([ones(numel(kmask),1), kTraj], [kmask(:), kTraj]);
imMask = true(imSize);
sysEqPara = sysEqPara_(fsTraj, imMask, dcf);

%% synthesize a toy dataset
Asyn = Gmri(kTraj, imMask);
fkPS = Asyn*reshape(PS, [], size(PS,4));
ukPS = bsxfun(@times, kmask(:)~=0, fkPS);

pSize = [7, 7]; % do play with this, :)

fP_direct = sosCombine(reshape(Asyn'*bsxfun(@times, dcf(:), fkPS), size(PS)));
uP_direct = test_direct(ukPS, usTraj, kcalib, pSize, imSize, sysEqPara);

fP_dogrid = sosCombine(PS);
uP_dogrid = test_dogrid(ukPS, usTraj, kcalib, pSize, imSize);

%% plots
figure,
subplot(221), imagesc(fP_direct); title('fP\_direct');
subplot(222), imagesc(uP_direct); title('uP\_direct');

subplot(223), imagesc(fP_dogrid); title('fP\_dogrid');
subplot(224), imagesc(uP_dogrid); title('uP\_dogrid');

%%
% keyboard
end

%%
function [uP] = test_direct(ukPS, usTraj, kcalib, pSize, imSize, sysEqPara)
% Directly reconstruct under-sampled k-space
%OUTPUTs:
% - uP (nx, ny, nz), combined final image
% - PS (nx, ny, nz, nc), reconstructed coil images
% - kPSG, reconstructed k-space

gMDL = grappaPrep(usTraj, kcalib, pSize, imSize, 'sysEqPara',sysEqPara);
[uP, PS, kPSG, ~] = GrappaNaive(ukPS, gMDL);

end

%%
function [uP] = test_dogrid(ukPS, usTraj, kcalib, pSize, imSize)
% Directly reconstruct under-sampled k-space
%OUTPUTs:
% - uP (nx, ny, nz), combined final image
% - PS (nx, ny, nz, nc), reconstructed coil images
% - kPSG, reconstructed k-space

sampled = usTraj(:,1) ~= 0;
usTraj = usTraj(sampled, :);
ukPS   = ukPS(sampled,:);

gMDL = grappaPrep(usTraj, kcalib, pSize, imSize, 'doGrid',true);
[uP, PS, kPSG, ~] = GrappaNaive(ukPS, gMDL);

end

