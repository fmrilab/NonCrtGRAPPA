function P = sosCombine(PS)
% function P = sosCombine(PS)
% sqrt sum-of-squares combining coils images into a single image.
%INPUTS:
% - PS (nx, ny, nz, nc), coil images;
%OUTPUTS:
% - P  (nx, ny, nz), the combined image.
if size(PS,4) == 1, P = abs(PS); return; end
P = sqrt(sum(abs(PS).^2,4));
end

