function FT_fatrix2 = form_FT(sMask_sTraj, m, zMap, dt_ti)
% form fatrix2 Fourier Transform object
% This fn does NOT assume the returned fatrix2 object to be of any certain type,
% e.g., Gmri, etc. Therefore, only access the fatrix2 obj properties at your own
% risk.
%
%INPUTS:
% - sMask_sTraj, sTraj in cycle/FOV
% - m (nx, ny, nz), boolean, imMask
%OPTIONAL:
% - zMap (nx, ny, nz), b0Map or relax_map + 2i*pi*b0Map; Hz, if b0Map
% - dt_ti (1,) or (nt,), Secs, dwell time (dt) or time-order (ti) of sTraj input
%OUTPUTS:
% - FT_fatrix2 (nk, nnz(m)) fatrix2 obj, nk: #sampled k point of sMask_sTraj.

if nargin == 0, test(); return; end

if ~nonEmptyVarChk('zMap'), [zMap,ti,dt_ti] = deal([]);
elseif isreal(zMap) && any(zMap(:) >= 0), zMap = 2i*pi*zMap; % b0Map only
% else,  % input zMap is already the Gmri assumed relax_map + 2i*pi*b0Map
end
if ~isempty(zMap) && ~nonEmptyVarChk('dt_ti'), dt_ti = envMR('get','dt'); end

if islogical(sMask_sTraj)   % sMask input
  if ~isempty(zMap)||~isempty(dt_ti), warning('zMap etc ignored for sMask'); end
  sMask = sMask_sTraj;
  FT_fatrix2 = form_FT_crt(sMask, m);
elseif isfloat(sMask_sTraj) % sTraj input, time ordered
  sTraj = sMask_sTraj;
  k = sTraj(:, 2:end);
  k = k(~~sTraj(:,1), :); % masking
  
  L = 6; % default used by Gmri
  if isrow(m) || iscolumn(m), Nd = numel(m);
  else,                       Nd = size(m);
  end
  Jd = 6*ones(1,numel(Nd));
  n_shift = Nd/2; % not sure if n_shift fails for odd Nd
  nufftArgs = {Nd, Jd, 2*Nd, n_shift, 'table', 2^10, 'minmax:kb'};
  
  if isempty(zMap)
    omega = 2*pi*bsxfun(@times, k, 1./Nd);
    FT_fatrix2 = Gnufft(m, [{omega}, nufftArgs]);
  else
    if ~isempty(dt_ti) % true only when zMap is given
      if isscalar(dt_ti), ti = (0:size(sTraj,1)-1)*dt_ti; % input dt_ti is dt
      else, ti = dt_ti; % input dt_ti is ti
      end
    end
    
    FT_fatrix2 = Gmri(k,m, 'ti',ti,'zmap',2i*pi*zMap, 'L',L, 'nufft',nufftArgs);
  end
else, error('not supported sMask_sTraj format');
end

end

function FT_fatrix2 = form_FT_crt(sMask, imMask)
% input output format mimicing Gmri
assert(isequal(size(sMask), size(imMask))); % dimension must match

[arg.dims, arg.imMask, arg.sMask] = deal(size(imMask), imMask, sMask);

FT_fatrix2 = fatrix2('mask',imMask, 'odim',nnz(sMask), 'does_many',true ...
             , 'forw',@forw_crt, 'back',@back_crt, 'arg',arg);
end

%%
function y = forw_crt(arg, x)
%INPUTS
% - x (nnz(arg.imMask), ...) or (size(imMask), ...)
%OUTPUTS
% - y (nnz(arg.sMask), ...)

if size(x, 1) == nnz(arg.imMask), x = embed(x, arg.imMask); end

if ismatrix(arg.imMask), x = reshape(x, arg.dims(1), arg.dims(2), 1, []); end

y = reshape(mfft3(x), prod(arg.dims), []);
y = y(arg.sMask(:), :);

end

function x = back_crt(arg, y)
%INPUTS
% - y (nnz(arg.sMask), ...)
%OUTPUTS
% - x (nnz(arg.imMask), ...)

y = embed(y, arg.sMask);
if ismatrix(arg.sMask), y = reshape(y, arg.dims(1), arg.dims(2), 1, []); end
% De-normalize, consistent w/ Gmri & Gnufft. This fn does adjoint, not inverse.
x = nnz(arg.sMask)*mifft3(y);
x = embed(reshape(x, nnz(arg.imMask), []), arg.imMask);

end

%%
function test()
imSize = [128,128];
P = phantom(imSize(1));
imMask = true(imSize);
% sMask test
sMask = true(imSize);
F = form_FT(sMask, imMask); % should just be a DFT matrix

figure;
subplot(131);
imagesc(abs(P - mifft3(reshape(F*P, imSize)))), colorbar;

% sTraj tests
sTraj = mask2mloc(sMask); % cycle/FOV
F = form_FT(sTraj, imMask);

subplot(132); % larger error (1e-6) due to single precision of Gmri
imagesc(abs(P - mifft3(reshape(F*P, imSize)))), colorbar;

zMap = zeros(imSize);
F = form_FT(sTraj, imMask, zMap);

subplot(133);
imagesc(abs(P - mifft3(reshape(F*P, imSize)))), colorbar;

disp('form_FT.test() done');

end
