function [ SIFT] = calculate_SIFT( GT, Recon,temporal_mean)
%CALCULATE_SIFT calculates Euclidean distance between scale invariant 
%feature transforms of between a ground truth image (GT) and a 
%reconstructed image (Recon). The frames are determined from the ground
%truth.
%
%   Sophie Schauman January 2019 - sophie.schauman@dtc.ox.ac.uk
%
%   INPUTS: 
%   GT            - A 5D ground truth image (X x Y x Z x t x v), spatial 
%                   dimensions, time, and vessel components (VE-ASL)
%   Recon         - A 5D VE-ASL image to be compared with GT
%   temporal_mean - boolean (1 - average over time, 0 - calculate SSIM for
%                   each frame separately)
%
%   OUTPUTS:
%   SIFT          - List of structural similarity values for each vessel 
%                   component (and time frame). (Wang et al 2003, 
%                   Multiscale structural similarity for image quality 
%                   assessment). Both the GT and Recon images are
%                   normalised by their maximum value, so that the dynamic
%                   range is [0, 1] for each SSIM comparison.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if temporal_mean % average across temporal dimension

    SIFT = zeros(size(Recon, 5),1);

    for v = 1:size(Recon, 5) %vessel
        Recon_temp = single(abs(mean(Recon(:,:,:,:,v),4))./...
            max(max(abs(mean(Recon(:,:,:,:,v),4)))));
        GT_temp = single(abs(mean(GT(:,:,:,:,v),4))./...
            max(max(abs(mean(GT(:,:,:,:,v),4)))));
        
        [F_GT, D_GT] = vl_sift(GT_temp);
        [~, D_Recon] = vl_sift(Recon_temp, 'Frames', F_GT);
        
        
        for n = 1:size(D_GT,2)
            SIFT(v) = SIFT(v) + norm(single(D_GT(:,n)-D_Recon(:,n)));
        end
        SIFT(v) = SIFT(v)/size(D_GT,2);
    end
        
else % calculate for each time frame separately
    SIFT = zeros(size(Recon, 5),size(Recon, 4));
    
    for v = 1:size(Recon, 5) %vessel
        for t = 1:size(Recon, 4) %time
            Recon_temp = single(abs(Recon(:,:,:,t,v))./...
                max(max(abs(Recon(:,:,:,t,v)))));
            GT_temp = single(abs(GT(:,:,:,t,v))./...
                max(max(abs(GT(:,:,:,t,v)))));
        
            [F_GT, D_GT] = vl_sift(GT_temp);
            [~, D_Recon] = vl_sift(Recon_temp, 'Frames', F_GT);
        
        
            for n = 1:size(D_GT,2)
                SIFT(v,t) = SIFT(v,t) + norm(single(D_GT(:,n)-D_Recon(:,n)));
            end
            SIFT(v,t) = SIFT(v,t)/size(D_GT,2);
        end  
    end
    


end

