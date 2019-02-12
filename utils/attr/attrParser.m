function [res, extra_c, unmatched] = attrParser(dflt, argPairs_c, KeepUnmatched)
% A naive way to parse argument
% It may not be a good idea enable this function to parse an obj. For objects,
% setting an attribute is allowed to affect its another attribute. Therefore the
% order of setting, which is unknown to this function, matters.
if nargin == 0, test(); return; end

p = inputParser; % inputParser has no constructor.
if isobject(dflt), KeepUnmatched = false; end
if ~exist('KeepUnmatched','var'), KeepUnmatched = true; end

p.KeepUnmatched = KeepUnmatched;

attr_c = attrs(dflt);
values_c = getattrs(dflt, attr_c);

% may adapt fancier init in the future
for ia = 1:numel(attr_c), addOptional(p, attr_c{ia}, values_c{ia}); end

parse(p, argPairs_c{:});
[res, unmatched] = deal(p.Results, p.Unmatched);
extra_c = reshape(attrs(unmatched), 1, []); % ensuring row cell vector
if ~isempty(extra_c)
  extra_c = reshape([extra_c; getattrs(unmatched, extra_c)], 1, []);
end

end

function test()
disp('##  attrParser.test() ##');
[x,y,z] = deal(3,4,5);
[a.x, a.y] = deal(1,2);
arg_pairs_c = {'x', x, 'y', y};

[~, extra_c, unmatched] = attrParser(a, arg_pairs_c);
disp('--  test: w/o extra:  --');
disp('extra_c: (expect no output for this)');
disp(extra_c);
disp('unmatched:');
disp(unmatched);

KeepUnmatched = false;
try [res, extra_c] = attrParser(a, [arg_pairs_c, {'z', z}], KeepUnmatched);
catch ME % MATLAB ERROR
  if strcmp(ME.identifier,'MATLAB:InputParser:UnmatchedParameter')
    [res, extra_c, unmatched] = attrParser(a, [arg_pairs_c, {'z', z}], true);
  else, disp('attrParser.test() failed'); rethrow(ME);
  end
end
disp('--  test: w/ extra:  --');
disp('matched:');
disp(res);
assert(isequal(res, struct('x',x,'y',y)), 'attrParser: res test failed');
disp('extra_c:');
disp(extra_c);
assert(isequal(extra_c, {'z',z}), 'attrParser: extra_c test failed');
disp('unmatched:');
disp(unmatched);
assert(isequal(unmatched, struct('z',z)), 'attrParser: unmatched test failed');

disp('-- attrParser.test() passed  --');

end
