function mdl = saveMDL(mdl, f)
% save mdl onto disk, then return a matfile() handle of it
if numel(f)<5 || ~strcmp(f(end-3:end), '.mat'), f = [f,'.mat']; end
if exist(f, 'file')
  warning('Existing mdl matfile detected, renaming it to *_old.mat.');
  disp('press any key to continue');
  pause; % ensure the user is aware of renaming
  movefile(f, [f(1:end-4), '_old.mat']);
  disp('renaming done');
end

save(f, '-struct', 'mdl', '-v7.3'); % mdl is used here
mdl = matfile(f, 'Writable', true);

end

