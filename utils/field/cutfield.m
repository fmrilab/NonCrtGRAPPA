function [clip, st] = cutfield(st, fieldsReq, fieldsOpt)

if exist('fieldsOpt', 'var')
  valid = false(numel(fieldsOpt),1);
  for ii = 1:numel(fieldsOpt)
    if isattr(st, fieldsOpt{ii}), valid(ii) = true; end
  end
  fieldsReq = [fieldsReq, fieldsOpt(valid)];
end

for ii = 1:numel(fieldsReq)
  clip.(fieldsReq{ii}) = st.(fieldsReq{ii});
end
st = rmfield(st, fieldsReq);

end
