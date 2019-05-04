function sysEqPara = sysEqPara_(sMask_sTraj, imMask, dcf, ti,zMap)
% returns a struct('F', ,'dcf', ,'imMask')
% where .F is a fatrix2 object that process NUFFT
% Notice that .F' is not equivalent to iFT functions is Matlab. It does not
% scale down the results. One can divide .F by sqrt(nk) to get a "unitary" NUFFT
% operator.

if nargin == 0, test(); return; end

if ~nonEmptyVarChk('imMask')
  if isfloat(sMask_sTraj), imMask = true(size(mloc2mask(sMask_sTraj)));
  else, imMask = true(size(sMask_sTraj));
  end
end

if ~exist('dcf','var')
  if islogical(sMask_sTraj), dcf = ones(nnz(sMask_sTraj),1);
  elseif isfloat(sMask_sTraj)
    sTraj = sMask_sTraj;
    if numel(sTraj(1,2:end)) == 3 % 3D k-space
      warning('dcf input not found, calculation is expensive, F5 to continue');
      keyboard
    end
    dcf = mdcf_voronoi(sTraj(:,2:end));
  else, error('not supported sMask_sTraj format');
  end
end
if isfloat(sMask_sTraj) && numel(dcf) == numel(sMask_sTraj(:,1))
  % application may input dcf for full-readout, here pick out sampled locs
  disp('picking dcf for sampled locs');
  dcf = dcf(~~sMask_sTraj(:,1));
end

if ~(exist('ti','var') && exist('zMap','var')), [ti,zMap] = deal([]); end

sysEqPara.F = form_FT(sMask_sTraj, imMask, ti, zMap);
sysEqPara.dcf = dcf;
sysEqPara.imMask = imMask;

f = {'F', 'dcf', 'imMask'}; % output should contain this fields
sysEqPara = chkattrs(sysEqPara, f);
end

function test()

kMask = true(3,3);
[kx, ky] = ndgrid(-1:1, -1:1);
sTraj = [kMask(:), kx(:), ky(:)]; % cycle/FOV
dcf = ones(numel(kMask),1);
imMask = true(size(kMask));

sysEqPara1 = sysEqPara_(sTraj, imMask, dcf);
sysEqPara2 = sysEqPara_(kMask, imMask, dcf);

P = ones(size(imMask));

kP_fatrix2ed1 = sysEqPara1.F * P;
kP_fatrix2ed2 = sysEqPara2.F * P;
kP_matrixed  = fft2(P);

P_fatrix2ed1 = sysEqPara1.F'*kP_fatrix2ed1;
P_fatrix2ed2 = sysEqPara2.F'*kP_fatrix2ed2;
P_matrixed  = ifft2(kP_matrixed);

disp('kP_fatrix2ed:')
disp([kP_fatrix2ed1, kP_fatrix2ed2]);
disp('kP_matrixed:');
disp(kP_matrixed(:));

disp('P_fatrix2ed:')
disp([P_fatrix2ed1, P_fatrix2ed2]);
disp('P_matrixed:');
disp(P_matrixed(:));

disp('sysEqPara_():');
disp('- Notice the scaling difference.');
disp('- Check out the documentation if this scaling surprises you');
end
