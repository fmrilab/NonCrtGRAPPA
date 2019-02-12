function x = mifft2(x, nx, ny)
% x: 0-freq centered
% res: image center centered

if nargin == 1
  x = fftshift2(ifft2(ifftshift2(x)));
else
  x = fftshift2(ifft2(ifftshift2(x), nx, ny));
end

end
