function FT_fatrix2 = form_FT(sMask_sTraj, m, zMap, dt_ti)
% form fatrix2 Fourier Transform object
%INPUTS:
% - sMask_sTraj
% - m (nx, ny, nz), boolean, imMask
%OPTIONAL:
% - zMap (nx, ny, nz), b0Map or relax_map + 2i*pi*b0Map; Hz, if b0Map
% - dt_ti (1,) or (nt,), Secs, dwell time (dt) or time-order (ti) of sTraj input

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
  
  if ~isempty(dt_ti) % true only when zMap is given
    if isscalar(dt_ti), ti = (0:size(sTraj,1)-1)*dt_ti; % input dt_ti is dt
    else, ti = dt_ti; % input dt_ti is ti
    end
  end
  
  L = 6; % default used by Gmri
  imSize = size(m);
  if numel(imSize) == 2, imSize = [imSize, 1]; end

  if imSize(3)~=1, Nd = imSize; else, Nd = imSize(1:2); end
  Jd = 6*ones(size(Nd));
  n_shift = Nd/2; % not sure if n_shift fails for odd Nd
  nufftArgs = {Nd, Jd, 2*Nd, n_shift, 'table', 2^10, 'minmax:kb'};
  
  FT_fatrix2 = Gmri(k,m, 'ti',ti, 'zmap',2i*pi*zMap, 'L',L, 'nufft',nufftArgs);
else, error('not supported sMask_sTraj format');
end

end

function FT_fatrix2 = form_FT_crt(sMask, imMask)
% input output format mimicing Gmri
sDim = size(sMask);
assert(isequal(sDim, size(imMask))); % dimension must match

if all(imMask), fn_imMask = @(x)x;
else, fn_imMask = @(x)bsxfun(@times, x, imMask);
end
fn_col = @(yy)reshape(yy, numel(sMask), []);
mo_sMask = @(yy)yy(sMask(:),:);
% ismatrix checks if input is 2D
if ismatrix(sMask), mi_sMask = @(y)reshape(embed(y,sMask),sDim(1),sDim(2),1,[]);
else, mi_sMask = @(y)embed(y, sMask);
end

forw = @(ob, x)mo_sMask(fn_col(mfft3(fn_imMask(x))));
back = @(ob, y)nnz(sMask)*fn_imMask(mifft3(mi_sMask(y))); % denormalized ifft

FT_fatrix2 = fatrix2('imask',imMask, 'odim',nnz(sMask), 'does_many',true ...
             , 'forw',forw, 'back',back);

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
