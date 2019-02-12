function structo = chkfield(structi, fieldso, fields_ignore)
if ~exist('fields_ignore', 'var'), fields_ignore = {}; end
fieldsi = fieldnames(structi);
diff = setdiff(fieldsi, [fieldso, fields_ignore]);
if ~isempty(diff)
  disp(diff);
  error('unexpected fields');
end
structo = structi;
for ii = 1:numel(fieldso)
  % this handles struct, class and 'matlab.io.MatFile'
  % for the extremely expensive matfile query, don't change the logic below
  if isattr(structo, fieldso{ii}), continue; end
  structo.(fieldso{ii}) = [];
end
end

