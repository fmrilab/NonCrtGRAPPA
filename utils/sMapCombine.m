function P = sMapCombine(PS, sMap)
% function P = sMapCombine(PS, sMap)
% Linearly combining coils images into a single image, using coils sMap.
%INPUTS:
% - PS (nx, ny, nz, nc, ...), coil images. 
% - sMap (nx, ny, nz, nc), coils sensitivity maps
%OUTPUTS:
% - P  (nx, ny, nz), the combined image.

denom = sum(conj(sMap).*sMap,4); % sum up across coil
nom = sum(bsxfun(@times, conj(sMap), PS), 4);
P = bsxfun(@rdivide, nom, denom);
P(isnan(P)) = 0;
end
