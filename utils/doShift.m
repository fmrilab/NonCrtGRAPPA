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
  error('Unsupported non-integer shifts, please open a GitHub issue if needed');
end

end

