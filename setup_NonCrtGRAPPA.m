
function genericSetup(doSavePath) %#ok<FNDEF> % generic, named intentionally
% This function always assume itself locates in the root dir of the package
if ~exist('doSavePath','var'), doSavePath = false; end
theDir = fileparts(mfilename('fullpath'));
prevDir = cd(theDir);
% mex setup
cd('./GRAPPA/private');
if isempty(which('LS_fft_mex'))
  mex -largeArrayDims -lmwlapack LS_fft_mex.c
end
cd(theDir);

%% generate path
% don't use pwd in genpath, as names in dName_rm_c can appear in it.
dName_s = genpath('.'); % _s: string, not string type though
% appending path w/ filesep makes fn_mtchd no need to worry partial matching
dName_s = strrep(dName_s, pathsep, [filesep,pathsep]);
dName_c = strsplit(dName_s, pathsep)'; % _c: cell
dName_c(end) = []; % dName_s ends w/ pathsep, trailing dName_c an empty {}, rm.

% customize dName_rm_c for different usages
dName_rm_c = {'demo', 'private', '.git', 'test', 'tests', 'arch', 'back'};

fn_mtchd = @(x,y)~isempty(strfind(x,[y,filesep])); %#ok<STREMP>
for iName = 1:numel(dName_rm_c)
  mask = cellfun(@(x)fn_mtchd(x, dName_rm_c{iName}), dName_c);
  dName_c(mask) = [];
end

%% add path
warning('off', 'MATLAB:mpath:packageDirectoriesNotAllowedOnPath');
% 1st char, '.', and last char, filesep, removed before addpath
for iName = 1:numel(dName_c), addpath([theDir, dName_c{iName}(2:end-1)]); end
warning('on', 'MATLAB:mpath:packageDirectoriesNotAllowedOnPath');

%% save?
if doSavePath, savepath; end

%% head back
cd(prevDir);
end
