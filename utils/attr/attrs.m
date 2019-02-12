function attrs_c = attrs(arg1, type)
%OUTPUTS:
% - attrs_c (#attr,), names of fields or properties

if isstruct(arg1)
  attrs_c = fieldnames(arg1);
elseif isobject(arg1)
  attrs_c = properties(arg1);
  if nonEmptyVarChk('type') % getting attrs of specified 'type'
    metaObj  = metaclass(arg1);
    objPropMetaProp = metaObj.PropertyList;
    attrs_c = attrs_c([objPropMetaProp.(type)]);
  end
else,  error('input must be either struct or object');
end

end
