function varargout = getattrs(arg1, attrs_c, isPairOut)
% function value_c = getattribute(arg1, attrName_c)
% get values of specified attributes of the input.
%INPUTS:
% - arg1 (1,), struct or object
% - attr_c (1,) or (#attr,), char or cell of chars of attribute names
% - isPairOut [t/F], paired output, e.g. output as {'a', a, 'b', b}
%OUTPUTS:
% - value_c (1,), (#attrName,) or (2*#attrName,), value or cell of values
if nargin == 0, test(); return; end
if ~exist('isPairOut', 'var'); isPairOut = false; end
if isempty(attrs_c), varargout = {{}}; return; end % shortcut

if ischar(attrs_c)
  varargout = {arg1.(attrs_c)}; % varargout requries wrapping into 1 cell
  attrs_c = {attrs_c};
elseif iscell(attrs_c)
  varargout = cell(size(attrs_c));
  % likely too expensive to invoke parfor
  for ia = 1:numel(varargout), varargout{ia} = arg1.(attrs_c{ia}); end
  if nargout <= 1, varargout = {varargout};
  else, assert(nargout == numel(varargout), 'getattrs: unmatched #OUT');
  end
else, error('getattrs: unexpected attrs_c input type');
end

if isPairOut, varargout = {reshape([attrs_c(:)'; [varargout{:}]], 1, [])}; end

end

function test()
f_c = {'x', 'y', 'z'};
v_c = {3,   4,   5};
a = cell2struct(v_c, f_c, 2);

assert(isequal(getattrs(a,f_c), v_c), 'default output test failed');

res_c = getattrs(a, f_c);
assert(isequal(res_c, v_c), 'single output test failed');

res_c = cell(1,numel(f_c));
[res_c{:}] = getattrs(a, f_c);
assert(isequal(res_c, v_c), 'separate output test failed');

res_c = cell(1,numel(f_c));
try [res_c{1:end-1}] = getattrs(a, f_c); %#ok<NASGU>
catch ME
  assert(strcmp(ME.message,'getattrs: unmatched #OUT'),'#OUT test failed');
end

res_c = getattrs(a, f_c, true);
assert(isequal(res_c, reshape({f_c{:}; v_c{:}}, 1, [])), 'PairOut test failed');

disp('getattrs.test() passed');
end
