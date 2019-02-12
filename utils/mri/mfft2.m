function x = mfft2(x, nx, ny)
% x: image center centered
% res: 0-freq centered 

if nargin == 1
  x = fftshift2(fft2(ifftshift2(x)));
else
  x = fftshift2(fft2(ifftshift2(x), nx, ny));
end

end
