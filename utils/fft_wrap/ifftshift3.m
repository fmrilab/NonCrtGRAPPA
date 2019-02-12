function res = ifftshift3(x)
% gt3 = false;
% if numel(size(x)) > 3, [gt3, nd] = deal(true, size(x)); end
% % trick: this is faster than ifftshift by indexing along dim>3
% if gt3, x = reshape(x, nd(1), nd(2), nd(3), []); end
% x = ifftshift(x);
% if gt3
%   fftshift(x, 4);
%   x = reshape(x, nd);
% end

% [~,~,~,nd4] = size(x);
% if nd4 == 1, x = ifftshift(x);
% else, for ii = 1:nd4, x(:,:,:,ii) = ifftshift(x(:,:,:,ii)); end
% end

ns = ceil(size(x)/2);
% out-place circshift is faster than in-place
res = circshift(x, ns(1:min(numel(ns),3)));

end
