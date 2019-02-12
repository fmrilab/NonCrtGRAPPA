function fpcgPara = fpcgPara_(imMask, varargin)

fpcgPara.x0 = [];
fpcgPara.tol = 1e-6;
fpcgPara.niter = 30;
fpcgPara.lam = 0;
fpcgPara.niterXtra = 5;
fpcgPara.Reg = Cdiff(imMask, 'order', 1);
fpcgPara.W = 1;

fpcgPara = attrParser(fpcgPara, varargin);

f = {'x0','tol','niter','niterXtra','lam','Reg','W'};
fpcgPara = chkfield(fpcgPara, f);
end
