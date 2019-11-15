function [ kdata, E ] = simulate_acq(GT, kspace, sens, H, noise_sd, seed, sens_from_data, ncoils  ) 
%SIMULATE_ACQ Simulate VE-ASL acquisition
%
%   Sophie Schauman July 2018 - sophie.schauman@dtc.ox.ac.uk
%
%   INPUTS: 
%   GT            - Ground truth complex images (5D - (X x Y x Z x t x v))
%   kspace        - 4D matrix 
%                  (nKPoints x nTPoints x  nEncodings x nSpatialDimensions)
%   sens          - coil sensitivities (4D - (X x Y x Z x nCoils)
%   H             - Vessel-encoding matrix (nEncodings x nVesselComponents)
%   noise_sd      - Standard deviation of noise added to each point in
%                   k-space
%   seed          - seed to use in the random number generator
%   sens_from_data- true if sensitivity maps should be estimated from noisy
%                   data
%   ncoils        - how many coils to compress to 
%
%   OUTPUTS:
%   kdata         - 4D matrix (nKPoints x nTPoints x  nEncodings x nCoils)
%   E             - Encoding operator for modelling acquisition/recon
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 8
    ncoils = size(sens,4); % default = no compression
    sens_from_data = 1;

end
if nargin < 7
    sens_from_data = 1;
end
if nargin < 6
    seed = 1;
end

rng(seed)
% Radial Golden Angle acquisition operator
E   =   xfm_NUFFT_VEASL(size(GT),sens,[],kspace,H);
weight = repmat(E.w,[1 1 1 size(sens, 4)] );

kNoiseless = ((E*GT)./weight)/E.norm;

weight = weight(:,:,:,1:ncoils);

% SNR^2 in k-space is signal power: rms((E*XGT)./E.w).^2, divided by noise
% power: rms(n).^2

n = noise_sd*(randn(E.dsize)+1i*randn(E.dsize))./sqrt(2);
SNR = rms(kNoiseless(:))/noise_sd

kdata = (kNoiseless + n);


%% estimate sensitivity based on noisy 'acquired' data
% because this function was originally built for acquired data the data
% needs to be unweighted going in and reweighted afterwards
if sens_from_data
    disp(['estimating sensitivities from data (compressed to ' num2str(ncoils) ' coils)'])
    [kdata_cc] = compress_coils(kdata, ncoils);
    E = get_acq_operator_from_data( kdata_cc, kspace, size(GT), H) ;
else
    disp(['using provided coil sensitivities'])
    kdata_cc = kdata;
    E = get_acq_operator_from_data( kdata, kspace, size(GT), H, sens) ;
end


kdata = kdata_cc.*weight*E.norm;
end




