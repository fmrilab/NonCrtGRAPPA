function res = mifft(x, n, dim)
% x: image center centered
% res: 0-freq centered, (size(x,dim)+1)/2 (odd), (size(x,dim)/2+1) (even)
if ~exist('n', 'var'), n = []; end
if ~exist('dim', 'var')
  if isrow(x), dim = 2;
  else         dim = 1;
  end
end
res = fftshift(ifft(ifftshift(x,dim),n,dim),dim);
end