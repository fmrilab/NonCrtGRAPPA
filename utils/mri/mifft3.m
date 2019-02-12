function x = mifft3(x, nx, ny, nz)
% x: 0-freq centered
% res: image center centered

if nargin == 1
  x = fftshift3(ifft3(ifftshift3(x)));
else
  x = fftshift3(ifft3(ifftshift3(x), nx, ny, nz));
end

end
