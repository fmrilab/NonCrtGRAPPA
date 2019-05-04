function [dst, dst_new, dst_unmatched] = mrgattrs(src, dst)
% mrg attrs in src into dst
% This program is not supposed to operate on any object as input
%INPUTS:
% - src st, source structure
% - dst st, destination structure
%OUTPUTS:
% - dst st, structure w/ all fields/attrs in src merged into dst
% - dst_new st, structure w/ fields/attrs in dst updated by those in scr
% - dst_unmatched st, structure w/ fields/attrs in src unmatched to dst input

if nargin == 0, test(); return; end

[aName_sc, aName_dc] = deal(attrs(src), attrs(dst));
for iName = 1:numel(aName_sc)
  aName_i = aName_sc{iName};
  dst.(aName_i) = src.(aName_i);
end

[dst_new, dst_unmatched] = cutattrs(dst, aName_dc);

end

%%
function test()

isSetEqual = @(a_c, b_c)isempty([setdiff(a_c, b_c), setdiff(b_c, a_c)]);

src = struct('a',1, 'b',2);
dst = struct('a',3, 'd',4);

[dst_all, dst_new, dst_unmatched] = mrgattrs(src, dst);

[aName_sc, aName_dc] = deal(attrs(src), attrs(dst));
aName_ac = unique([aName_sc, aName_dc]); % _a: all
[aName_dc_all, aName_dc_new, aName_dc_unmatched] = ...
  deal(attrs(dst_all), attrs(dst_new), attrs(dst_unmatched));

tf_dst_all = isSetEqual(aName_dc_all, aName_ac);
for iName = 1:numel(aName_dc_all)
  aName_i = aName_dc_all{iName};
  if ismember(aName_i, aName_sc)
    tf_dst_all = tf_dst_all && dst_all.(aName_i) == src.(aName_i);
  else
    tf_dst_all = tf_dst_all && dst_all.(aName_i) == dst.(aName_i);
  end
end
assert(tf_dst_all, 'dst_all stucture failed');

tf_dst_new = isSetEqual(aName_dc_new, aName_dc);
for iName = 1:numel(aName_dc_new)
  aName_i = aName_dc_new{iName};
  tf_dst_new = tf_dst_new && dst_new.(aName_i) == dst_new.(aName_i);
end
assert(tf_dst_new, 'dst_new stucture failed');

tf_dst_unmatched = isSetEqual(aName_dc_unmatched, setdiff(aName_sc, aName_dc));
for iName = 1:numel(aName_dc_unmatched)
  aName_i = aName_dc_unmatched{iName};
  tf_dst_unmatched = tf_dst_unmatched && dst_unmatched.(aName_i)==src.(aName_i);
end
assert(tf_dst_unmatched, 'dst_unmatched stucture failed');

disp([mfilename('fullpath'), '.test() passed']);
end

