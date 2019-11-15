function [ GT ] = sim_pre_process( GT, dimensions,varargin)
%SIM_PRE_PROCESS Pre process simulation ground truth (make right size)
%
%   Sophie Schauman July 2018 - sophie.schauman@dtc.ox.ac.uk
%
%   INPUTS: 
%   filePath      - file to use as ground truth (e.g. Thomas Okell's 2016
%                   NMR paper data. 5D files (X x Y x Z x t x VesselComps))
%   dimensions    - what size the data should be (if any dimension is 1 it
%                   will take a sum across that dimension instead of 
%                   slicing through it.
%   (VesselCombine- list of which vessel channels should be combined)
%
%   OUTPUTS:
%   GT            - a 5D matrix to use as ground truth in simulations
%                   (X x Y x Z x t x VesselComps)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read in extra inputs for reading data 
p   =   inputParser;
p.addParameter('VesselCombine',       []);
p.parse(varargin{:});
p   =   p.Results;

combineV = p.VesselCombine;

%% Resize spatial dimensions
if dimensions(1)>=size(GT,1)
    GT = padarray(GT, [(dimensions(1)-size(GT,1))/2]);
elseif dimensions(1) == 1
    GT = sum(GT,1);
else
    margin = (size(GT,1)-dimensions(1))/2;
    GT = GT(margin+1:end-margin, :,:,:,:);
end
if dimensions(2)>=size(GT,2)
    GT = padarray(GT, [0 (dimensions(2)-size(GT,2))/2]);
elseif dimensions(2) == 1
    GT = sum(GT,2);
else
    margin = (size(GT,2)-dimensions(2))/2;
    GT = GT(:,margin+1:end-margin,:,:,:);
end

if dimensions(3)>=size(GT,3)
    GT = padarray(GT, [0 0 (dimensions(3)-size(GT,3))/2]);
elseif dimensions(3) == 1
    GT = sum(GT,3);
else
    margin = (size(GT,3)-dimensions(3))/2;
    GT = GT(:,:,margin+1:end-margin,:,:);   
end

%% Resize temporal dimension
if size(GT,4) == dimensions(4)
    %do nothing
elseif size(GT,4) > dimensions(4)
    % cut off time
    GT = GT(:,:,:,1:dimensions(4),:);
elseif size(GT,4) < dimensions(4)
    % add empty frames
    frames = zeros([dimensions(1:3),dimensions(4)-size(GT,4),size(GT,5)]);
    GT = cat(4, GT, frames);
end

%% Resize vessel component dimension
if ~isempty(combineV)
    combined = sum(GT(:,:,:,:,combineV),5)./length(combineV);
    GT(:,:,:,:,combineV(1)) = combined;
    GT(:,:,:,:,combineV(2:end)) = [];
end
GT = abs(GT)/max(abs(GT(:)));
