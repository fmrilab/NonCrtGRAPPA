function cInd = ctrInd(Nd, varargin)
Nd = [Nd, varargin{:}];

cInd = sum( (ctrSub(Nd)-1).*[1, cumprod(Nd(1:end-1))] ) + 1;

