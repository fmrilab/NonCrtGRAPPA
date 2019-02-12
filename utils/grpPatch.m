function patch_c = grpPatch(patch_c, varargin)
%INPUTS:
% - patch_c (#_patch_types, 2) cell: take iith row as an example
%     (ii,1):
%       sMask case: (#ngb, 1), mask identifying sampled ngbs;
%       sTraj case: (#ngb,nd), ngb-ctr shifts (single float), rowsort
%     (ii,2):
%       array composed as [ctrInds, ngbInds];
%OPTIONAL:
% - doSift (t/F) sift the neighbors by a Cartesian grid
%OUTPUTS:
% - patch_c, same format as its INPUT form, (sifted), unique'd & sorted
% see also: getPatch
arg.doSift = false;
arg = attrParser(arg, varargin);

[shifts_rc, inds_c] = deal(patch_c(:,1), patch_c(:,2));
if arg.doSift
  parfor ii = 1:size(patch_c,1) % parfor
    [shifts_i, inds_i] = deal(shifts_rc{ii}, inds_c{ii});
    [~, si, ~] = unique(round(shifts_i), 'rows');
    [shifts_rc{ii}, inds_c{ii}] = deal(shifts_i(si,:), inds_i([1;si+1]));
  end
end

% group by number of shifts in patch, then check duplicate ones by 'uniquetol'
ngbCnt = cellfun(@numel, inds_c)-1;
[ngrp, shifts_rcc, inds_cc] = grpIntoCell(ngbCnt, shifts_rc, inds_c);
% check 1st non-empty patch to see whether grouping sMask or sTraj
if isfloat(shifts_rc{find(ngbCnt,1)})
  [isMask, ndim] = deal(false, size(shifts_rc{find(ngbCnt,1)},2));
else
  [isMask, ndim] = deal(true, []); % parfor initialization
end

% sub-group by distinct patterns
parfor ii = 1:ngrp % parfor
  [shifts_i, inds_i] = deal(cat(2,shifts_rcc{ii}{:}), cat(1,inds_cc{ii}{:}));
  if ~isMask, shifts_i=reshape(single(shifts_i),ndim*(size(inds_i,2)-1),[]); end
  shifts_i = shifts_i.';
  
  % Identify distinct shifts patterns, remove replicates w/in tol.
  % Shifts is mostly of 1e0 units, tol==5e-3 units should be enough for float.
  % uniquetol only work for float, but faster than unique when tol is dflt value
  itol = 5e-3/1e-6;
  if isMask, [shifts_i, ~, typelbl_i] = unique(shifts_i, 'rows');
  else
    [shifts_i, ~, typelbl_i] = uniquetol(shifts_i/itol, 'ByRows', true);
    shifts_i = shifts_i * itol;
  end
  [~,inds_cc{ii}] = grpIntoCell(typelbl_i,inds_i);
  
  shifts_rc_i = num2cell(shifts_i.', 1).';
  if isMask, shifts_rcc{ii} = shifts_rc_i;
  else shifts_rcc{ii} = cellfun(@(x)reshape(x,[],ndim), shifts_rc_i,'Uni',0);
  end
end

patch_c = [cat(1, shifts_rcc{:}), cat(1, inds_cc{:})];
end

function [ngrp, varargout] = grpIntoCell(grpBy,varargin)
%INPUTS:
% target (n,): column vector that will be sorted
% varargin: each cell is a matrix of (n,...) size
%OUTPUTS:
% varargout: each cell is the corresponding vin cell reordered by sorted target
[grpBy, srtInd_tmp] = sort(grpBy);
grpCnt = histc(grpBy(:), unique(grpBy, 'stable'));
ngrp = numel(grpCnt);
nv = numel(varargin);
varargout = cell(nv,1);
for ii = 1:nv, varargout{ii} = mat2cell(varargin{ii}(srtInd_tmp,:),grpCnt); end
end

