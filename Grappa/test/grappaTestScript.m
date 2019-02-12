function grappaTestScript()
% The .mat testing files is not in the repo due to their size.
% Instead, they can be found at www-personal.umich.edu/~tianrluo/Phantoms/
% (or one can generate ones own files based on the following comments)
% Here's some general comments about the .mat files:
%   phmA (struct):
%     .kPS  (nx, ny, nz, nc), kspace Phantom w/ SMaps (indexed by nc(oil))
%     .PS   (nx, ny, nz, nc), Phantom w/ SMaps
%     .sMap (nx, ny, nz, nc), sensitivity maps
%   sprA (struct):
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

testsos();
testsos2crt();

end

%%
function testsos()

phmA = load('phantom_dat.mat');
kcalib = phmA.kPS(86:115, 86:115, :,:);

% kSpiral2D is spiral-out readout
sprA = load('kSpiral2D.mat');
pSize  = [7,7];

fkspace = 200*[real(sprA.k(:)), imag(sprA.k(:))];
kmask = zeros(size(sprA.k));


kmask(1:75,:) = 1; % densely sampled the center
kmask(:,1:2:end) = 1;

sTraj = [kmask(:), fkspace];

imSize = sprA.imSize;
Tik = 5e-7;

% synthesize simu data
imMask = true(imSize);
Asyn = Gmri(fkspace, imMask);

PS = phmA.PS;
kPS = Asyn*reshape(PS, [], size(PS,4));

ukPS = bsxfun(@times, kmask(:)~=0, kPS);

% sysEqPara_() will calc the dcf
sysEqPara = sysEqPara_([ones(size(kmask(:))), fkspace], true(imSize));

kPS2PS_fn = [];

PS2P_fn = [];
% PS2P_fn = @(PS)sMapCombine(PS, phmA.sMap);

doSift = true;
gMDL = grappaPrep(sTraj,kcalib,pSize,imSize, 'doSift',doSift...
  , 'sysEqPara',sysEqPara, 'kPS2PS_fn',kPS2PS_fn, 'PS2P_fn',PS2P_fn, 'Tik', Tik);

[uP, PS, kPSG, arg] = GrappaNaive(ukPS, gMDL);

figure, im(uP);

% keyboard
end
%%
function testsos2crt()
phmA = load('phantom_dat.mat');
kcalib = phmA.kPS(91:110, 91:110, :,:);

% kSpiral2D is spiral-out readout
sprA = load('kSpiral2D.mat');
pSize  = [7,7];

[xx, yy] = ndgrid(-99:100, -99:100);

maskCtr = sqrt(xx.^2 + yy.^2) <= 10; % densely sampled the center
kTrajCtr = [xx(maskCtr), yy(maskCtr)];

R = 2;
kTrajtmp = 200*sprA.k(55:end, 1:R:end); % remove spiral Ctr, and downsample
kTrajXtr = [real(kTrajtmp(:)), imag(kTrajtmp(:))];

m_fn = @(x)ones(size(x,1),1);

kTraj = [      kTrajCtr;          kTrajXtr];
kmask = [m_fn(kTrajCtr); -1*m_fn(kTrajXtr)];

sTraj = [kmask, kTraj];

imSize = sprA.imSize;
Tik = 5e-7;

% synthesize simu data
kTrajSyn = [kTrajCtr; kTrajXtr];
imMask = true(imSize);
Asyn = Gmri(kTrajSyn, imMask);

fPS = phmA.PS;
ukPSi = Asyn*reshape(fPS, [], size(fPS,4));

gMDL = grappaPrep(sTraj, kcalib, pSize, imSize, 'Tik', Tik, 'doGrid',true);

fP = sosCombine(phmA.PS); % raw (not synth) as ground truth

[uP, ~, ukPSo, arg] = GrappaNaive(ukPSi, gMDL); % bingo

figure, im(uP);
end

