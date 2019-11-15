function [  ] = acquisition_pipeline_f(MeasFileID, SubjID, imsize, H, R, PrepNos, frames, CoilComp, coils, lambda1, lambda2, iter, step, S, GTfile, maskfile,lowmem, savefolderpath) 
%ACQUISITION_PIPELINE_F Function for processing data from radial GA VEASL 
%
%   Sophie Schauman August 2018 - sophie.schauman@dtc.ox.ac.uk
%
%   INPUTS: MeasFileID - .dat TWIX file ID (Siemens) that contains the
%                        acquired data (if left empty, user can choose file
%                        manually
%           SubjID  - identifier/subject number
%           imsize  - length 5 vector containing the size of the
%                       reconstruction (x, y, z, t, v)
%           H       - Hadamard or other matrix for vessel encoding
%           R       - Acceleration factor
%           PrepNos - Which preparation repeats to use (alternative to R)
%           frames  - Which frames to import and reconstruct
%           CoilComp- Number of virtual coils to compress to
%           coils   - pre-calculated coil sensitivity maps
%                     (Nx x Ny x Nz x Nc)
%           lambda1 - weighting for l1 regularisation term
%           lambda2 - weighting for l2 temporal smoothness term
%           iter   - number of iterations in FISTA
%           step   - stepsize in FISTA
%           S      - sparsifying transform for l1-regularisation. Must be
%                    object that can be applied with the *-operator (integer,
%                    2D-matrix, or specifically designed matlab object)
%           GTfile - .mat file with a pre-reconstructed ground truth image
%                    (Nx x Ny x Nz x Nt x Nvc)
%           maskfile - mask for evaluation of reconstruction quality 
%                      (Nx x Ny x Nz x Nt x Nvc)
%           lowmem   - true if the reconstruction should use the table
%                      based NUFFT 
%           savefolderpath - where to save the reconstruction and analysis
%                           result
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 0. Checking inputs
if nargin < 18
    savefolderpath = 'data/reconstructions/in-vivo/';
end
if nargin < 17
    lowmem = false;
end
if nargin < 16
    maskfile = [];
    mask = repmat(ones(imsize), 1, 1, 1, size(H,2));
end
if nargin < 15
    GTfile = [];
    maskfile = [];
    clear mask;
end
if nargin < 14
    S = 1;
end
if nargin < 13
    step = 0.1;
end
if nargin < 12
    iter = 100;
end
if nargin < 11
    lambda2 = 0;
end
if nargin < 10
    lambda1 = 1e-5;
end
if nargin < 9
    coils = [];
end
if nargin < 8
    Coil_Comp = [];
end

if nargin < 7
    frames = 1:imsize(4);
end

if nargin < 6
    PrepNos = [];
end
if nargin < 5
    R = 1;
end
if nargin < 4
    error(' not enough input arguments provided, [MeasFileID, SubjID, imsize, H] all required');
end

if lowmem
    disp('Table lookup used for NUFFT')
end

%% 1. Read in data and create acquisition model
[kdata, E] = acq_pre_process(MeasFileID, imsize, H, frames, 'R', R, 'PrepNos', PrepNos, 'CoilComp', CoilComp, 'Coils', coils, 'lowmem', lowmem);

% for saving file in the end PrepNos can't be empty
if isempty(PrepNos)
    PrepNos = [1 0];
end

%% 2. Reconstruct
% FISTA Reconstruction
Recon = fista_L2s(kdata, E, S, lambda1, lambda2, [E.Nd, E.Nt, E.Nvc], iter, step); 


% %% 3. Analyse
if ~isempty(GTfile) %GT provided
    
    GTStruct = load(GTfile);
    GT = GTStruct.Recon;
    
    if ~isempty(maskfile)
        maskStruct = load(maskfile);
        mask = maskStruct.mask;
    end

    SSIM_all_frames = calculate_SSIM(GT, Recon, 0);
    SSIM_temporal_mean = calculate_SSIM(GT, Recon, 1);

    NRMSE_all_frames = calculate_NRMSE(GT, Recon, 0);
    NRMSE_temporal_mean = calculate_NRMSE(GT, Recon, 1);

    [temp_corr_all, temp_corr_masked] = calculate_temp_corr(GT, Recon, mask);
    
    HFEN = calculate_HFEN(GT,Recon);
    SIFT_all_frames = calculate_SIFT(GT,Recon,0);
    SIFT_temporal_mean = calculate_SIFT(GT, Recon, 1);
    
    %temporary fix!! Use same mask as for temp_corr but dilate to get
    %spat_corr
    SE = strel('disk',1);
    for v = 1:3
        mask(:,:,:,v) = imdilate(mask(:,:,:,v),SE);
    end
    spat_corr = calculate_spat_corr(GT,Recon,1,mask,0);
    
                            
end
%% 4. Save
mkdir(savefolderpath);
savefilename = sprintf('Subj%1.0f_Recon%1.0fVC_FOV%3.0f_R%0.1f_%ito%i_Lambda1_%0.6f_Lambda2_%0.6f.mat', SubjID, size(H,2), imsize(1), R, PrepNos(1), PrepNos(end) ,lambda1, lambda2);
savefilename = [savefolderpath avoidOverwrite(savefilename,savefolderpath',3,1)];
if ~isempty(GTfile)
    save( savefilename,'Recon', 'SSIM_all_frames', 'SSIM_temporal_mean',...
        'NRMSE_all_frames','NRMSE_temporal_mean','temp_corr_all', ...
        'temp_corr_masked','HFEN', 'SIFT_temporal_mean', 'SIFT_all_frames', 'spat_corr','-v7.3');
else
    save( savefilename,'Recon', '-v7.3')
end