function fPS = fPSrec(fkPS, sysEqPara, fpcgPara)
%*Basically a wrapper of qpwls_pcg1, may be modified for other iter meth tho.
% reconstruct fully sampled coil images
% fkPS: data, (nx, ny, nz, nc) for sMask format; (nk, nc) for sTraj format
% eqPara struct: used for sTraj format data recon
%  .imMask (nx, ny, nz): OPT, logical mask of the object in the imaging system
%  .imSize (nd,): OPT, image size, necessary when not providing .imMask
%  .F (nc*nk, nc*nxyz), OPT, image system matrix
%  .sTraj (nk, ndim): OPT, sampling trajectory, necessary when not providing .F
%  .dcf (nc*nk, nc*nk), OPT, density compensation matrix, useful in Ajoint rec
% fpcgPara struct: used for recon with cg algorithm
%  .W (nc*nk, nc*nk), OPT, preconditioner in cg alg
%  .lam, .Reg, regulization pair
%  .niter (1,), number of iteration, w/ .W be applied
%  .niterXtra (1,), number of iteration, w/o .W be applied
% 

[imSize, datSize] = deal(size(sysEqPara.imMask), size(fkPS));
if numel(datSize)>=numel(imSize) && isequal(datSize(1:numel(imSize)), imSize)
  fPS = mifft3(fkPS);
  return;
end

imMask = sysEqPara.imMask;
fkPS = squeeze(fkPS);
[~, nc, nf] = size(fkPS);
% so after embedding, nc with on 4th-dim, nf on 5th-dim
fPS = zeros([nnz(imMask), ones(1, 3-ndims(imMask)), nc, nf]);

if ~nonEmptyVarChk('fpcgPara') % ajoint recon, w/ dcf
  fkPS = bsxfun(@times, sysEqPara.dcf(:), fkPS)/size(sysEqPara.F, 1);
  Fh = sysEqPara.F';
  parfor ii = 1:nc*nf
    fPS(:,ii) = Fh*fkPS(:,ii);
  end
else
  if isempty(fpcgPara.x0)
    x0 = ones(0, nc*nf); % THIS IS AN EMPTY ARRAY
  else
    fpcgPara.x0 = reshape(fpcgPara.x0, [], nc);
    x0 = fpcgPara.x0(imMask, :);
  end
  
  [lam, Reg] = deal(fpcgPara.lam, fpcgPara.Reg);
  if lam && ~isemtpy(Reg), lReg = lam*Reg; else, lReg = 0; end
  
  [niter, niterX] = deal(fpcgPara.niter, fpcgPara.niterXtra);
  [F, W] = deal(sysEqPara.F, fpcgPara.W);
  
  % this implementation uses cg as iter method, w/ dcf as preconditioner
  fn = @(xi,y,W,ni)qpwls_pcg1(xi(:), F,W,y,lReg, 'niter',ni, 'key',0);
  
  parfor ii = 1:nc*nf
    fPS(:,ii) = fn(x0(:, ii), fkPS(:,ii), W, niter);
  end
  
  if niterX > 0
    parfor ii = 1:nc*nf
      fPS(:,ii) = fn(fPS(:,ii), fkPS(:,ii), 1, niterX);
    end
  end
end

% for 2d data, ensure coil-dim is the 4th dim
if ismatrix(imMask), fPS = reshape(fPS, nnz(imMask), 1, nc, nf);
else               , fPS = reshape(fPS, nnz(imMask), nc, nf);
end

fPS = embed(fPS, imMask);

end
