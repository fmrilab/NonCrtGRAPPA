function x = ifft3(x, nx, ny, nz)

if numel(size(x)) == 3
  if nargin == 1, x = ifftn(x); else, x = ifftn(x, nx,ny,nz); end; return;
end

if nargin == 1
  x = ifft(ifft2(x), [], 3);
else
  x = ifft(ifft2(x, nx, ny), nz, 3); % self padded by zero, lovely
end

end
