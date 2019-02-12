function subs = ctrSub(Nd, varargin)
% As a separate fn, ensure consistent behaviour of getting subscripts to the
% center of a Nd Matrix
% 0 1 2 3 4 5 6
% 0 1 2 2 3 3 4
% in consistent w/ fftshift and ifftshift, where the location subs [1]
% shifted to is assumed as the center

Nd = [Nd, varargin{:}];

subs = ceil((Nd+1)/2);
subs(Nd == 0) = 0;

end

