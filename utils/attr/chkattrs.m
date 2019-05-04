function output = chkattrs(input, attrs_oc, attrs_ignore_c)
%INPUTS:
% - input, st
% - attrs_oc, cell array, list of expected attrs in the output
% - attrs_ignore_c, cell array, list of attrs ignored when checking existence
%OUTPUTS:
% - output, st, input with missing attrs in attrs_oc added as empty, '[]'.

if ~exist('attrs_ignore_c', 'var'), attrs_ignore_c = {}; end
attrs_c = attrs(input);
diff = setdiff(attrs_c, [attrs_oc, attrs_ignore_c]);
if ~isempty(diff)
  disp(diff);
  error('unexpected attributes');
end
output = input;
for ii = 1:numel(attrs_oc)
  % this handles struct, class and 'matlab.io.MatFile'
  % for the extremely expensive matfile query, don't change the logic below
  if isattr(output, attrs_oc{ii}), continue; end
  output.(attrs_oc{ii}) = [];
end
end

