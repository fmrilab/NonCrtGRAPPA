function res = fftshift2(x)
% gt2 = false;
% if numel(size(x)) > 2, [gt2, nd] = deal(true, size(x)); end
% % trick: this is faster than fftshift by indexing along dim>3
% if gt2, x = reshape(x, nd(1), nd(2), []); end
% x = fftshift(x);
% if gt2
%   ifftshift(x, 3);
%   x = reshape(x, nd);
% end

ns = floor(size(x)/2);
res = circshift(x, ns(1:2)); % not-in-place circshift is faster than in-place

end
