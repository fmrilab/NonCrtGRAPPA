function res = patch2row(x, pMask, bMethod, pcSub)
% res = patch2row(x, pMask, doCrop, cSub)
%Inputs
% - x (nx,ny,nz, ..., na,nb,nc, ...), data to be put into patches
% - pMask (px,py,pz, ...), Mask of patch
% - bMethod str: boundary method
% - doCrop ([T]/F), if doCrop, cut out boundaries, deal w/ circshift
% - cSub (ndim,) subscript define the center of the pMask
%Outputs
% - res (npxyz, na,nb,nc, ..., nnz(pMask)), npxyz: #patches in (nx,ny,nz,...)

%% prep and find shifts
pSize = size(pMask);
if ~nonEmptyVarChk('pcSub'),   pcSub = ctrSub(pSize); end
if ~nonEmptyVarChk('bMethod'), bMethod = 'crop'; end

if strcmpi(bMethod, 'pad0')
  x = padarray(padarray(x,pcSub-1,0,'pre'), pSize-pcSub,0,'post');
end

xSize = size(x);

pMaskSub1c = cell(1,numel(pSize));
[pMaskSub1c{:}] = ind2sub(pSize, find(pMask));
pMaskSub1 = cell2mat(pMaskSub1c); % (nnz(pMask), numel(pSize))

shifts = int8(bsxfun(@minus, pcSub, pMaskSub1));

%% do the shifts
% (nx,ny,nz,...,na,nb,nc,...) -> (nx,ny,nz,...,na,nb,nc,...,nnz(pMask))
restmp = doShift(x, shifts);

%% reshape into rows
npdim = (~iscolumn(pMask))* numel(pSize) + iscolumn(pMask);
nxdim = (~iscolumn(x))    * numel(xSize) + iscolumn(x);
ndimX = nxdim - npdim;

if ~strcmpi(bMethod, 'circ') % both 'crop' and 'pad0' will get in here
  cropSubsc = cell(1, npdim);
  for ii = 1:npdim % xSize must be larger than pSize along all dim
    cropSubsc{ii} = pcSub(ii)+(0:(xSize(ii)-pSize(ii)));
  end
  colons_c = cell(1,ndimX +1); % +1 for shifts
  colons_c(:) = {':'};
  cropSubsc = [cropSubsc, colons_c];
  res = restmp(cropSubsc{:}); % do the crop
else
  res = restmp;
end

if ndimX > 0
  shapeTmp = [{[]}, num2cell(xSize((end-ndimX+1):end)), {nnz(pMask)}];
  res = squeeze(reshape(res, shapeTmp{:}));
else
  shapeTmp = [{[]}, {nnz(pMask)}];
  res = reshape(res, shapeTmp{:});
end

end

