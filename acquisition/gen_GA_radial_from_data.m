function [ kspace ] = gen_GA_radial_from_data(kinfo, varargin)
%GEN_GA_RADIAL_FROM_DATA Generate 2D radial golden angle sampling positions
%from parsed twix file for vessel encoded ASL images.
%
%
%   Sophie Schauman July 2018 - sophie.schauman@dtc.ox.ac.uk
%   Based on code by Thomas W Okell.
%
%   INPUTS: 
%   kinfo         - twix struct containing parsed information from .dat file
%                   generated with mapVBVD by Philipp Ehses.
%   (R            - acceleration factor)
%   (PrepNos      - list of specified prep numbers to use)
%
%   OUTPUTS:
%   kspace        - 4D matrix 
%                  (nKPoints x nTPoints x  nEncodings x nSpatialDimensions)
%   DEPENDENCIES:
%   GoldenMeans3D.m
%   GoldenRatios3D.m
%%%
% Read in extra inputs for choosing the right spokes 
p   =   inputParser;
p.addParameter('R',       []);
p.addParameter('PrepNos', []);
p.addParameter('Frames',[]);
p.parse(varargin{:});
p   =   p.Results;

R = p.R;
PrepNos = p.PrepNos;
FrameNos = p.Frames;
%%%
% Determine data type
SeqType = kinfo.hdr.MeasYaps.sWiPMemBlock.alFree{12+1};
switch SeqType
    case 2
        Seq2D = true;

    case 4
        Seq2D = false;

    otherwise
        error(['SeqType = ' num2str(kinfo.hdr.MeasYaps.sWiPMemBlock.alFree{12+1}) ': not implemented!']);
end
    
    
    
AcqSpokesPerPhase = kinfo.hdr.Config.NSeg; % 
TotSpokesPerFrame = kinfo.hdr.Config.RawLin;
Nyquist_limit_spokes = ceil(kinfo.hdr.Config.BaseResolution*pi/2);
if ~Seq2D
    Nyquist_limit_spokes = ceil(kinfo.hdr.Config.BaseResolution^2*pi/2);
end
NFrames = kinfo.hdr.Config.NPhs;
AcqPreps = TotSpokesPerFrame / AcqSpokesPerPhase;
NyquistPreps = ceil(Nyquist_limit_spokes / AcqSpokesPerPhase);
S = kinfo.image.dataSize; % Grab the input data size



% Preps to use
if ~isempty(PrepNos) % PrepNos provided
    ReconPreps = length(PrepNos);
elseif ~isempty(R)% PrepNos not provided, need R
    ReconPreps = ceil(NyquistPreps/R);
    PrepNos = 1:ReconPreps;
else
    R = 1;
    ReconPreps = NyquistPreps;
    PrepNos = 1:NyquistPreps;
end

if ~isempty(R) % R provided
    if R == NyquistPreps/ReconPreps %do nothing
    else % update R and throw warning
        R = NyquistPreps/ReconPreps;
        warning('The provided acceleration factor (R) does not use an integer number of acquisitions or does not match the user provided list of preps to use (PrepNos) . Re-calculating R.')
    end
else % R not provided (needs PrepNos)
    R = NyquistPreps/ReconPreps;  
end
    
% Calculate useful parameters
Nsegs = kinfo.hdr.Config.NSeg;
Nspokes = NyquistPreps*AcqSpokesPerPhase/R;
if Nspokes > TotSpokesPerFrame
    error('Not enough data for that acceleration factor!')
end
Nsampl = kinfo.hdr.Config.RawCol;
if isempty(FrameNos)
    FrameNos = 1:kinfo.hdr.Config.NPhs;
end

NVE = kinfo.hdr.Config.NAve;


if Seq2D
    kspace = single(zeros(Nspokes*Nsampl, 2, length(FrameNos)));
else
    kspace = single(zeros(Nspokes*Nsampl, 3, length(FrameNos)));
end

CGR = -(1-sqrt(5))/2; % Chan, MRM 2009
for ii = 1:length(FrameNos)

    NFrame = FrameNos(ii);

    % Identify the actual spoke numbers relative to the start of the
    % readout that we want to grab here
    LinNos = ((NFrame-1)*Nsegs+1):(NFrame*Nsegs);

    % Empty the index of relevant line and phase numbers
    LinIdx = []; PhsIdx = [];

    % Loop through the preps, adding appropriate indices
    for jj = PrepNos

        % Line number at the start of each cardiac phase
        StartLinNo = (jj-1)*Nsegs + 1;

        LinIdx = [LinIdx (StartLinNo+mod(LinNos-1,Nsegs))];
        PhsIdx = [PhsIdx (floor((LinNos-1)/Nsegs)+1)];

    end

    % Calculate the appropriate azimuthal and polar angles
    GRCounter = (LinIdx-1) + (PhsIdx-1) * Nsegs;


    if Seq2D
        Phi = mod( GRCounter * CGR * pi, 2*pi );
        Theta = [];
    else
        [Phi, Theta] = GoldenMeans3D(GRCounter,true);
    end


    % Construct the trajectory
    kspaceframe = gen_radial_traj(Phi,Nsampl, Theta);

    % Change signs to get final image in the right orientation
    if Seq2D
        kspaceframe = [-kspaceframe(:,1) kspaceframe(:,2)];
    else
        kspaceframe = [-kspaceframe(:,1) kspaceframe(:,2) kspaceframe(:,3)];
    end

    kspace(:,:,ii) = kspaceframe;

end

    kspace = permute(kspace,[1 3 4 2]);
    kspace = repmat(kspace,1, 1, NVE, 1);


