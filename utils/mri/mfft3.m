function x = mfft3(x, nx, ny, nz)
% x: image center centered
% res: 0-freq centered

if nargin == 1
  x = fftshift3(fft3(ifftshift3(x)));
else
  x = fftshift3(fft3(ifftshift3(x), nx, ny, nz));
end

end
