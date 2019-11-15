function [ temp_corr_all, temp_corr_masked] = calculate_temp_corr( GT, Recon, mask)
%CALCULATE_TEMP_CORR calculates the temporal correlation coefficient 
%between a ground truth image (GT) and a reconstructed image (Recon).
%
%   Sophie Schauman July 2018 - sophie.schauman@dtc.ox.ac.uk
%
%   INPUTS: 
%   GT            - A 5D ground truth image (X x Y x Z x t x v), spatial 
%                   dimensions, time, and vessel components (VE-ASL)
%   Recon         - A 5D VE-ASL image to be compared with GT
%   mask          - a 4D mask (X x Y x Z x v) for which voxels to include 
%                   in analysis for each vessel component.
%
%   OUTPUTS:
%   temp_corr_all  - Correlation coefficients for each voxel in each vessel
%                   component.
%   temp_corr_masked - list of correlation coefficients within
%                      vessel-specific masks.  Matrix size = nVoxels x v.
%                      Padded with NaN's where the number of voxels within
%                      the different masks differ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(Recon,5) ~= size(mask,4)
    if size(Recon,5) == 2 % if nonVE
        Recon1 = repmat(Recon(:,:,:,:,1), 1, 1, 1, 1, size(mask,4)-1);
        GT1 = repmat(GT(:,:,:,:,1), 1, 1, 1, 1, size(mask,4)-1);
        Recon = cat(5, Recon1, Recon(:,:,:,:,2));
        GT = cat(5, GT1, Recon(:,:,:,:,2));
    else
    error('mask dimensions does not match reconstruction!')
    end
end
    
temp_corr_all = zeros(size(mask));

for v = 1:size(Recon, 5)
    for x = 1:size(Recon,1)
        for y = 1:size(Recon,2)
            for z = 1:size(Recon,3)
                C = corrcoef(squeeze(abs(GT(x,y,z,:,v))),squeeze(abs(Recon(x,y,z,:,v))));
                temp_corr_all(x,y,z,v) = C(1,2);
            end
        end
    end
end
temp_corr_masked =[];
for v = 1:size(Recon, 5)
    masked = temp_corr_all(:,:,:,v);
    masked = masked(logical(mask(:,:,:,v)));
    masked(isnan(masked)) = 0;
    temp_corr_masked = catpad(2, temp_corr_masked, masked);
end

end

