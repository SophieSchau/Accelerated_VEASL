function [  ] = simulation_pipeline_f( GTFile, SimID, sensFile, reconSize, H, R, seed, lambda1, lambda2, iter, step, noise_sd, S, maskFile, nsens, savefolderpath, sens_from_data )
%SIMULATION_PIPELINE_F Function for simulating and processing data from radial GA VEASL 
%
%   Sophie Schauman August 2018 - sophie.schauman@dtc.ox.ac.uk
%
%   INPUTS: GTFile - .mat file that contains the ground truth VE-ASL image
%                    (Nx x Ny x Nz x Nt x Nvc)
%           SimID  - identifier/subject number
%           sensFile - .mat file containing sensitivity maps
%                      (Nx x Ny x Nz x Nc)
%           reconSize - length 5 vector containing the size of the
%                       reconstruction (x, y, z, t, v)
%           H      - Hadamard or other matrix for vessel encoding
%           R      - Acceleration factor
%           seed   - seed for random number generator for noise
%           lambda1 - weighting for l1 regularisation term
%           lambda2 - weighting for l2 temporal smoothness term
%           iter   - number of iterations in FISTA
%           step   - stepsize in FISTA
%           noise_sd - standard deviation of noise added to each point in
%                    k-space
%           S      - sparsifying transform for l1-regularisation. Must be
%                  object that can be applied with the *-operator (integer,
%                  2D-matrix, or specifically designed matlab object)
%           maskFile - mask for evaluation of reconstruction quality 
%                      (Nx x Ny x Nz x Nt x Nvc)
%           nsens - number of channels to compress data to
%           savefolderpath - where to save reconstruction and analysis
%                           results
%           sens_from_data - true if sensitivites should be derived from
%                            noisy data rather than using the generative 
%                            sensitivity profiles
%   
%   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0. Check and get inputs
if nargin < 17
    sens_from_data = true; % default = get sensitivity profiles from data 
end
if nargin < 16
    savefolderpath = 'data/reconstructions/sim/';
end

if nargin < 15
    nsens = 32; % default = no coil compression
    
end
if nargin < 14
    maskFile = [];
    mask = repmat(ones(reconSize), 1, 1, 1, size(H,2)); %default = no masking
end
if nargin < 13
    S = 1; % default = no sparsifying transform
end
if nargin < 12
    noise_sd = 0; % default = no noise
end
if nargin < 11
    step = 0.1;
end
if nargin < 10
    iter = 100;
end
if nargin < 9
    lambda2 = 0;
end
if nargin < 8
    lambda1 = 1e-5;
end
if nargin < 7
    seed = 1;
end
if nargin < 6
    R = 1;
end
if nargin < 5
    error(' not enough input arguments provided, [GTFile, SimID, sensFile, reconSize, H] all required');
end

F = load(GTFile);
GT = single(F.Recon);               % ground truth


%%%%%%
F = load(sensFile);
sens = F.sens;              % sensitivities
F = load(maskFile);
mask = F.mask;              % spatial mask for temporal analysis

sensSize = [reconSize(1:3), size(sens,4)];


%% 1. Preporcess raw simulation data and simulate acquisition
sens = sim_pre_process( sens, sensSize); % resize coil sensitivity profiles
sens = sens./max(sens(:));
GT = sim_pre_process( GT, reconSize); % resize GT
GT = GT./max(reshape(GT(:,:,:,:,size(H,2)),[],1))*0.001; % scaling to make similar to in-vivo data

% set trajectory
NGASpokes = ceil(max(reconSize(:))*pi/2/R);
Nsamps = max(reconSize(:))*2;
CGR = -(1-sqrt(5))/2;
Phi = mod( (0:NGASpokes*reconSize(4)-1)' * CGR * pi, 2*pi );
kspace = gen_radial_traj(Phi, Nsamps, []);
kspace = repmat(kspace,size(H,1),1);
kspace = reshape(kspace, NGASpokes*Nsamps, reconSize(4), size(H,1), []);

% simulate acquisition
[ kdata, E ] = simulate_acq(GT, kspace, sens, H, noise_sd, seed, sens_from_data, nsens  ) ;

%% 2. Reconstruct

% FISTA Reconstruction
Recon = fista_L2s(kdata, E, S, lambda1, lambda2, [E.Nd, E.Nt, E.Nvc], iter, step); 


%% 3. Analyse

SSIM_all_frames = calculate_SSIM(GT, Recon, 0);
SSIM_temporal_mean = calculate_SSIM(GT, Recon, 1);

NRMSE_all_frames = calculate_NRMSE(GT, Recon, 0);
NRMSE_temporal_mean = calculate_NRMSE(GT, Recon, 1);

HFEN = calculate_HFEN(GT,Recon);
SIFT_all_frames = calculate_SIFT(GT,Recon,0);
SIFT_temporal_mean = calculate_SIFT(GT, Recon, 1);


if size(GT,4) > 1
    [temp_corr_all, temp_corr_masked] = calculate_temp_corr(GT, Recon, mask);
else
    temp_corr_all = [];
    temp_corr_masked = [];
end

spat_corr = calculate_spat_corr(GT,Recon,1,mask,0);
%% 4. Save

mkdir(savefolderpath);
savefilename = sprintf('Subj%1.0f_Recon%1.0fVC_FOV%3.0f_R%0.1f_%i_Lambda1_%0.9f_Lambda2_%0.9f.mat', SimID, size(H,2), reconSize(1), R, seed ,lambda1, lambda2);
savefilename = [savefolderpath avoidOverwrite(savefilename,savefolderpath,3,1)];

save( savefilename,'Recon', 'SSIM_all_frames', 'SSIM_temporal_mean',...
    'NRMSE_all_frames','NRMSE_temporal_mean','temp_corr_all', 'temp_corr_masked','noise_sd', 'HFEN', 'SIFT_all_frames', 'SIFT_temporal_mean','spat_corr', '-v7.3');

end
