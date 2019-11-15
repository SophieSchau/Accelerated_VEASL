function [ SSIM ] = calculate_SSIM( GT, Recon, temporal_mean, mask )
%CALCULATE_SSIM calculates structural similarity index between a ground
%truth image (GT) and a reconstructed image (Recon).
%
%   Sophie Schauman July 2018 - sophie.schauman@dtc.ox.ac.uk
%
%   INPUTS: 
%   GT            - A 5D ground truth image (X x Y x Z x t x v), spatial 
%                   dimensions, time, and vessel components (VE-ASL)
%   Recon         - A 5D VE-ASL image to be compared with GT
%   temporal_mean - boolean (1 - average over time, 0 - calculate SSIM for
%                   each frame separately)
%
%   OUTPUTS:
%   SSIM          - List of structural similarity values for each vessel 
%                   component (and time frame). (Wang et al 2003, 
%                   Multiscale structural similarity for image quality 
%                   assessment). Both the GT and Recon images are
%                   normalised by their maximum value, so that the dynamic
%                   range is [0, 1] for each SSIM comparison.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4
    mask = ones(size(GT,1), size(GT,2), size(GT,3), size(GT,5));
end

if size(mask,4) ~= size(GT,5)
    temp(:,:,:,1) = logical(sum(mask(:,:,:,1:end-1),4));
    temp(:,:,:,2) = mask(:,:,:,end);
    mask = temp;
end


    
scaling_factor_Recon = squeeze(abs(mean(Recon(:,:,:,1,1:end-1),4))).*...
    squeeze(mask(:,:,:,1:end-1));
scaling_factor_Recon(scaling_factor_Recon==0) = [];
scaling_factor_Recon = prctile(scaling_factor_Recon,95);

scaling_factor_GT = squeeze(abs(mean(GT(:,:,:,1,1:end-1),4))).*...
    squeeze(mask(:,:,:,1:end-1));
scaling_factor_GT(scaling_factor_GT==0) = [];
scaling_factor_GT = prctile(scaling_factor_GT,95);

% SE = strel('disk',3  );
% for v = 1:size(mask,4)
%         mask(:,:,:,v) = imdilate(mask(:,:,:,v),SE);
% end


if temporal_mean % average across temporal dimension
    SSIM = zeros(size(Recon, 5),1);
    for v = 1:size(Recon, 5) %vessel
                        [~, ssim_map] = ssim(abs(mean(Recon(:,:,:,:,v),4))./...
                            scaling_factor_Recon,...
                            abs(mean(GT(:,:,:,:,v),4))./...
                            scaling_factor_GT);
                        ssim_map(~mask(:,:,:,v)) = [];
                        SSIM(v) = mean(ssim_map(:));
    end
else % calculate for each time frame separately
    SSIM = zeros(size(Recon, 5),size(Recon, 4));
    
    for v = 1:size(Recon, 5) %vessel
        for t = 1:size(Recon, 4) %time
            [~, ssim_map] = ssim(abs(Recon(:,:,1,t,v))./...
            scaling_factor_Recon,...
            abs(GT(:,:,1,t,v))./scaling_factor_GT);
        ssim_map(~mask(:,:,:,v)) = [];
        SSIM(v,t) = mean(ssim_map(:));
        end
    end
end

SSIM(isnan(SSIM)) = 0;

end
