function x = fft3(x, nx, ny, nz)
%% fast

if numel(size(x)) == 3
  if nargin == 1, x = fftn(x); else, x = fftn(x, nx,ny,nz); end; return;
end

if nargin == 1
  x = fft(fft2(x), [], 3);
else
  x = fft(fft2(x, nx, ny), nz, 3); % self padded by zero, lovely
end

%% slow
% ffn converts logical as double, assign double to logical x losses double
% gt3 = false;
% if numel(size(x)) > 3, [gt3, nd] = deal(true, size(x)); end
% % trick: this is faster than fft by indexing along dim>3
% if gt3, x = reshape(x, nd(1), nd(2), nd(3), []); end
% x = fftn(x);
% if gt3
%   x = ifft(x, [], 4);
%   x = reshape(x, nd);
% end

%% slower
% [~,~,~,n4] = size(x);
% if n4 == 1
%   x = fftn(x);
% else
%   for ii = 1:n4
%     x(:,:,:,ii) = fftn(x(:,:,:,ii));
%   end
%   
% end

end
