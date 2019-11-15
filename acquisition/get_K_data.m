function [ kdata, kinfo, R ] = get_K_data( MeasFileID, varargin)
%GET_K_DATA Get suitably undersampled kdata, from 2D radial golden angle 
%vessel encoded ASL images from Siemens raw .dat file 
%
%   Sophie Schauman July 2018 - sophie.schauman@dtc.ox.ac.uk
%   Based partially on code by Thomas W Okell.
%
%   INPUTS: 
%   MeasFileID - .dat file identifier. Can be filename with path or just
%                   the MID from the .dat filename if it is on the defined path.
%   (R            - acceleration factor)
%   (PrepNos      - list of specified prep numbers to use)
%   (Frames       - indecies of frames to read in)
%
%   OUTPUTS:
%   kdata         - 4D matrix (nKPoints x nTPoints x  nEncodings x nCoils)
%   kinfo         - twix object with parsed information from .dat file
%   R             - R used (might have changed from input if input doesn't
%                   work with data.
%
%   DEPENDENCIES:
%   mapVBVD by Philipp Ehses (philipp.ehses@tuebingen.mpg.de) for reading 
%   .dat files. Available on:
%   https://github.com/CIC-methods/FID-A/tree/master/inputOutput/mapVBVD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
    % Read in extra inputs for reading data 
    p   =   inputParser;
    p.addParameter('R',       []);
    p.addParameter('PrepNos', []);
    p.addParameter('Frames', []);
    p.parse(varargin{:});
    p   =   p.Results;
    
    R = p.R;
    PrepNos = p.PrepNos;
    Frames = p.Frames;
    
    % Read the raw data headers
    twix_obj = mapVBVD(MeasFileID,'ignoreSeg');
    kinfo = twix_obj;
    
    % Determine data type
    SeqType = twix_obj.hdr.MeasYaps.sWiPMemBlock.alFree{12+1};
    switch SeqType
        case 2
            Seq2D = true;

        case 4
            Seq2D = false;

        otherwise
            error(['SeqType = ' num2str(twix_obj.hdr.MeasYaps.sWiPMemBlock.alFree{12+1}) ': not implemented!']);
    end
    
    % Extract info from the header
    % 
    % Vocabulary:
    % Prep    - One cycle of VE-ASL encoding and acquisition
    % Phase   - Each prep readout is divided into phases of length TempResn
    % Frame   - Data from multiple preps in same phase
    % Spoke   - Line through k-space
    % Segment - Each Prep is divided into segments of length TR (one line).
    
    AcqSpokesPerPhase = twix_obj.hdr.Config.NSeg; % 
    TotSpokesPerFrame = twix_obj.hdr.Config.RawLin;
    Nyquist_limit_spokes = ceil(twix_obj.hdr.Config.BaseResolution*pi/2);
    if ~Seq2D
        Nyquist_limit_spokes = ceil(twix_obj.hdr.Config.BaseResolution^2*pi/2);
    end
    NFrames = twix_obj.hdr.Config.NPhs;
    AcqPreps = TotSpokesPerFrame / AcqSpokesPerPhase;
    NyquistPreps = ceil(Nyquist_limit_spokes / AcqSpokesPerPhase);
    S = twix_obj.image.dataSize; % Grab the input data size

    
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
    SpokesPerFrame = NyquistPreps*AcqSpokesPerPhase/R;
    if SpokesPerFrame > TotSpokesPerFrame
        error('Not enough data for that acceleration factor!')
    end
    
    disp(['Spokes per frame           = ' num2str(SpokesPerFrame)])
    disp(['Spokes needed for Nyquist  = ' num2str(Nyquist_limit_spokes)])
    disp(['Undersampling factor       = ' num2str(R)])
    disp(['Using ' num2str(ReconPreps) ' of ' num2str(AcqPreps) ' Preps'])
    
    % Initialise kdata and kspace if requested
    kdata = single(zeros([S(1:2), SpokesPerFrame, S(4:6) 1 S(8:end), length(Frames)]));   
if isempty(Frames)
    Frames = 1:NFrames;
end
    for ii = Frames

        Frame = ii;

        % Identify the actual spoke numbers relative to the start of the
        % readout that we want to grab here
        LinNos = ((Frame-1)*AcqSpokesPerPhase+1):(Frame*AcqSpokesPerPhase);

        % Empty the index of relevant line and phase numbers
        LinIdx = []; PhsIdx = [];

        % Loop through the preps, adding appropriate indices
        for jj = PrepNos

            % Line number at the start of each cardiac phase
            StartLinNo = (jj-1)*AcqSpokesPerPhase + 1;

            LinIdx = [LinIdx (StartLinNo+mod(LinNos-1,AcqSpokesPerPhase))];
            PhsIdx = [PhsIdx (floor((LinNos-1)/AcqSpokesPerPhase)+1)];

        end

        
        % Grab the raw data with the appropriate indices
        disp(['Reading in raw data from frame ' num2str(ii) ' of ' num2str(length(Frames)) '...'])
        kdataframe = single(zeros([S(1) S(2) length(LinIdx) S(4:6) 1 S(8:end)]));

        % Order of raw data:
        %  1) Columns
        %  2) Channels/Coils
        %  3) Lines
        %  4) Partitions
        %  5) Slices
        %  6) Averages
        %  7) (Cardiac-) Phases
        %  8) Contrasts/Echoes
        %  9) Measurements
        % 10) Sets
        % 11) Segments
        % 12) Ida
        % 13) Idb
        % 14) Idc
        % 15) Idd
        % 16) Ide

        for kk = 1:length(LinIdx)
            kdataframe(:,:,kk,:,:,:,:,:,:,:,:,:,:,:,:,:) = twix_obj.image(:,:,LinIdx(kk),:,:,:,PhsIdx(kk),:,:,:,:,:,:,:,:,:);
        end

        kdata(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,ii) = kdataframe;
    end

    kdata = permute(kdata,[1 3 5 17 6 2 4 7:16]);
    kdata = reshape(kdata, [], length(Frames), S(6), S(2));


end

