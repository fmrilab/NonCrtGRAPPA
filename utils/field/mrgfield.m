function dst = mrgfield(src, dst)
% Caution: this function also replaces dst fields w/ those contained in src

for fNamesi = fieldnames(src)', dst.(fNamesi{1}) = src.(fNamesi{1}); end

end