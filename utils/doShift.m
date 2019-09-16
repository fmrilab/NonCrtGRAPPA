function res = doShift(x, shifts, varargin)
% function res = doShift(x, shifts)
%   call by doShift(x, [sx1, sy1, ... ; sx2, sy2, ...; ...])
%   Shifts matrix by overlay linear phase in frequency domain
%Input:
% - x (nx, ny, ...), Nd array to be shifted
% - shifts (ns, Nd): shift along each dimension
%Output:
% - res (nx, ny, ..., ns), Nd-by-ns array after shift
%e.g.:
%>> tmp = zeros(5); tmp(13) = 1;
%>> doShift(tmp, [-1,-1; 1,1]);

dims = size(x);
if numel(varargin) ~= 0 && ~isscalar(shifts)
  error('Unexpected input type');
end
if numel(varargin) > 0, shifts = [shifts, cell2mat(varargin)]; end
if size(shifts,2) > numel(dims)
  error('Prescripted shift exceeds x-dimension');
end

if iscolumn(x), Nd = 1;
else            Nd = ndims(x);
end

if isinteger(shifts) || isequal(shifts(:), int8(shifts(:)))
  ns = size(shifts,1);
  if ns == 1, res = circshift(x, shifts); return; end
  res_c = cell([ones(1,Nd),ns]);
  for ishift = 1:ns
    res_c{ishift} = circshift(x, shifts(ishift, :));
  end
  % res = cell2mat(res_c);
  res = cat(Nd+1, res_c{:});
else
  fx = fftn(x);
  
  [~, phs] = phs4Shift(shifts, 'dims',size(x)); % exponential phs by dflt
  
  fx = bsxfun(@times, fx, phs);
  res = zeros(size(fx), class(x));
  
  colon_c = repmat({':'}, [1,Nd]);
  for ii = 1:size(fx, Nd+1)
    res(colon_c{:},ii) = ifftn(fx(colon_c{:}, ii));
  end
end

end

