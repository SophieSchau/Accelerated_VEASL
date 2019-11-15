function [ NRMSE ] = calculate_NRMSE( GT, Recon, temporal_mean )
%CALCULATE_NRMSE calculates normalised root mean square error between a 
%ground truth image (GT) and a reconstructed image (Recon).
%
%   Sophie Schauman July 2018 - sophie.schauman@dtc.ox.ac.uk
%
%   INPUTS: 
%   GT            - A 5D ground truth image (X x Y x Z x t x v), spatial 
%                   dimensions, time, and vessel components (VE-ASL)
%   Recon         - A 5D VE-ASL image to be compared with GT
%   temporal_mean - boolean (1 - average over time, 0 - calculate NRMSE for
%                   each frame separately)
%
%   OUTPUTS:
%   NRMSE          - List of normalised root mean square errors (normalised
%                    by the maximum in the vessel component (+ time frame)
%                    of the ground truth image).  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if temporal_mean % average across temporal dimension

    NRMSE = zeros(size(Recon, 5),1);

    for v = 1:size(Recon, 5) %vessel
        error = (abs(Recon(:,:,:,:,v))-abs(GT(:,:,:,:,v)))./...
                max(max(abs(GT(:,:,:,:,v))));
        NRMSE(v) = sqrt(mean(error(:).^2));
    end
else % calculate for each time frame separately
    NRMSE = zeros(size(Recon, 5),size(Recon, 4));
    
    for v = 1:size(Recon, 5) %vessel
        for t = 1:size(Recon, 4) %time
            error = (abs(Recon(:,:,:,t,v))-abs(GT(:,:,:,t,v)))./...
                max(max(abs(GT(:,:,:,t,v))));
            NRMSE(v,t) = sqrt(mean(error(:).^2));
        end
    end
end

end

