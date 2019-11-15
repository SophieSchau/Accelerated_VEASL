function [ kdata, E ] = acq_pre_process(MeasFileID, imsize, H, frames, varargin  ) 
%ACQ_PRE_PROCESS Pre process 2d Golden Angle radial VE-ASL data from 
%Siemens scanner
%
%   Sophie Schauman July 2018 - sophie.schauman@dtc.ox.ac.uk
%
%   INPUTS: 
%   MeasFileID - .dat file identifier. Can be filename with path or just
%                   the MID from the .dat filename if it is on the defined path.
%   imsize        - spatial dimensions of image (vector of length 3)
%   H             - Vessel-encoding matrix (nEncodings x nVesselComponents)
%   frames        - indecies of frames to read in and reconstruct
%   (R            - acceleration factor)
%   (PrepNos      - list of specified prep numbers to use)
%   (CoilComp     - number of virtual coils to use)
%   (Coils        - pre-calculated coils sensitivity matrices)
%
%   OUTPUTS:
%   kdata         - 4D matrix (nKPoints x nTPoints x  nEncodings x nCoils)
%   E             - Encoding operator for modelling acquisition
%
%   DEPENDENCIES:
%   get_K_data, gen_GA_radial_from_data, get_acq_operator_from_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read in extra inputs for reading data 
p   =   inputParser;
p.addParameter('R',       []);
p.addParameter('PrepNos', []);
p.addParameter('CoilComp', []);
p.addParameter('Coils', []);
p.addParameter('lowmem', false);
p.parse(varargin{:});
p   =   p.Results;

R = p.R;
PrepNos = p.PrepNos;
c = p.CoilComp;
coils = p.Coils;
    
%% Read raw data
[kdata, kinfo, R] = get_K_data(MeasFileID,'R', R, 'PrepNos', PrepNos, 'Frames', frames);

%% Coil compression
if exist('c', 'var')
    kdata = compress_coils(kdata, c);
end

%% Estimate acquisition operator
kspace = gen_GA_radial_from_data(kinfo, 'R', R, 'PrepNos', PrepNos, 'Frames', frames);

%%% THIS WILL BE A BUG IF A 5 ENCODING SCHEME EVER IS USED!! ONLY HERE TO
%%% DEAL WITH SUBJECT 1 THAT USED 5 ENCODINGS BUT WE ONLY WANT TO USE 1:4
%%% FOR VESSEL ENCODING or 1 AND 5 FOR NON-VE
if size(kdata,3) == 5
    if size(H, 2) == 4
        warning('5 cycle data treated as 4 cycles')
        kdata = kdata(:,:,1:4,:);
        kspace = kspace(:,:,1:4,:);
    elseif size(H, 2) == 2
        warning('5 cycle data treated as 2 cycles')
        kdata = kdata(:,:,[1,5],:);
        kspace = kspace(:,:,[1,5],:);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E = get_acq_operator_from_data( kdata, kspace, imsize, H, coils, 'lowmem', p.lowmem );

%% Pre-process
% Phase correction (line-by-line)
kdata = reshape(kdata, kinfo.hdr.Config.RawCol, [], size(kdata,2), size(kdata,3), size(kdata,4));

for ii = 1:size(kdata,2) % each line
    for jj = 1:size(kdata,3) % each frame
        for kk = 1:size(kdata,4) % each encoding
            p=angle(mean(reshape(kdata(:,ii,jj,kk,:).*conj(kdata(:,ii,jj,1,:)),[],1)));
            kdata(:,ii,jj,kk,:)=kdata(:,ii,jj,kk,:).*exp(-1j*p);
        end
    end
end

kdata = reshape(kdata, [], size(kdata,3), size(kdata,4), size(kdata,5));


% pre-weighting (to work the same way as simulations)
weight = repmat(E.w,1, 1, 1, E.Nc);
kdata = weight.*kdata;
end




