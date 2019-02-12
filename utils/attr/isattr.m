function tf = isattr(x, str)
% niao
  tf = isfield(x, str) || isprop(x, str);
end
