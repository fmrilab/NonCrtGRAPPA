function [clip, clipped] = cutattrs(input, aName_rc, aName_oc)
% cut/clip attrs requrested in aName_rc and valid ones in aName_oc
% This program is not supposed to operate on any object as input
%INPUTS:
% - input st, structure whos fields/attrs to be cutted/clipped
% - aName_rc (nReq,) cell, required attrs to be clipped
% - aName_oc (nOpt,) cell, optional attrs, if exist in input, to be clipped
%OUTPUTS:
% - clip st, structure containing clipped/cutted attrs
% - clipped st, input without attrs in clip removed

if nargin == 0, test(); return; end

if exist('aName_oc', 'var')
  valid = false(numel(aName_oc),1);
  for ii = 1:numel(aName_oc)
    if isattr(input, aName_oc{ii}), valid(ii) = true; end
  end
  aName_rc = [aName_rc, aName_oc(valid)];
end

clip = struct;
for iName = 1:numel(aName_rc)
  aName_i = aName_rc{iName};
  clip.(aName_i) = input.(aName_i);
end
clipped = rmattrs(input, aName_rc);

end

%%
function test()

isSetEqual = @(a_c, b_c)isempty([setdiff(a_c, b_c), setdiff(b_c, a_c)]);

st = struct('a',1, 'b',2, 'c',3, 'd',4, 'e',5);
aName_rc = {'a', 'b'};                   % _r: required; _c: cell
aName_oc = {'d', 'e', 'f'};              % _o: optional
aName_ac = unique([aName_rc, aName_oc]); % _a: all

[st_clip,      st_clipped]      = cutattrs(st, aName_rc, aName_oc);
[aName_clip_c, aName_clipped_c] = deal(attrs(st_clip), attrs(st_clipped));

tf_clip = isSetEqual(attrs(st_clip), intersect(aName_ac, attrs(st)));
for iName = 1:numel(aName_clip_c)
  aName_i = aName_clip_c{iName};
  tf_clip = tf_clip && (st_clip.(aName_i) == st.(aName_i));
end
assert(tf_clip, 'clip stucture failed');

tf_clipped = isSetEqual(attrs(st_clipped), setdiff(attrs(st), aName_ac));
for iName = 1:numel(aName_clipped_c)
  aName_i = aName_clipped_c{iName};
  tf_clipped = tf_clipped && (st_clipped.(aName_i) == st.(aName_i));
end
assert(tf_clipped, 'clipped stucture failed');

disp([mfilename('fullpath'), '.test() passed']);
end
