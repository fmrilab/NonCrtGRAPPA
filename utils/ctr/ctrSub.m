function cSub = ctrSub(Nd, varargin)
% As a separate fn, ensure consistent behaviour of getting subscripts to the
% center of a Nd Matrix
% 0 1 2 3 4 5 6
% 0 1 2 2 3 3 4
% for consistency w/ fftshift and ifftshift, where the location cSub[1] shifted
% to is assumed as the center

Nd = [Nd, varargin{:}];

cSub = ceil((Nd+1)/2);
cSub(Nd == 0) = 0;

end

