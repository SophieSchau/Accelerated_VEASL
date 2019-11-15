function [ HFEN ] = calculate_HFEN( GT, Recon )
%CALCULATE_HFEN calculates high frequency error norm between a 
%ground truth image (GT) and a reconstructed image (Recon). HFEN is 
%computed as the l2-norm of the result obtained by LoG filtering the 
%difference between the recon- structed and reference images.
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
%   HFEN          -  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HFEN = zeros(size(Recon, 5),1);
imsize = size(Recon);
filt = fspecial('log');
if imsize(3) >1
    error('Cannot calculate HFEN for 3D volume (yet). Only 2D implemented.')
end
for v = 1:size(Recon, 5) %vessel
    % normalise intensity range,take absolute values
    Recon_temp = abs(Recon(:,:,:,:,v))/max(max(max(max(abs(Recon(:,:,:,:,v))))));
    GT_temp = abs(GT(:,:,:,:,v))/max(max(max(max(abs(GT(:,:,:,:,v))))));
    image_error = GT_temp - Recon_temp;
    LoG_error = imfilter(reshape(image_error, imsize(1), imsize(2),[]), filt);
    HFEN(v) = norm(LoG_error(:));
end

end

