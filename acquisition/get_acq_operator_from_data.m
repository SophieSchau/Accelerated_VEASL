function [ E ] = get_acq_operator_from_data( kdata, kspace, imsize, H, coils, varargin ) 
%GET_ACQ_OPERATOR_FROM_DATA Generates a linear operator E for modelling 
%acquisition
%
%   Sophie Schauman July 2018 - sophie.schauman@dtc.ox.ac.uk
%   Based partially on code by Mark Chiew.
%
%   INPUTS: 
%   kdata         - 4D matrix (nKPoints x nTPoints x  nEncodings x nCoils)
%   kspace        - 4D matrix 
%                  (nKPoints x nTPoints x  nEncodings x nSpatialDimensions)
%   imsize        - spatial dimensions of image (vector of length 3)
%   H             - Vessel-encoding matrix (nEncodings x nVesselComponents)
%   coils         - precalculated coil sensitivity maps (5D matrix
%                   (X by Y by Z by 1 by nCoils)
%   varargin      - 'lowmem' is a logical that implements the 'table'
%                    method in the NUFFT's if true. It reduces memory
%                    burden, but slows down the application of the NUFFT.
%
%   OUTPUTS:
%   E             - Encoding operator for modelling acquisition
%
%   DEPENDENCIES:
%   xfm_NUFFT_VEASL,
%   adaptive_estimate_sens
%   fista
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p   =   inputParser;
p.addParamValue('lowmem',   false,                   @islogical);   
p.parse(varargin{:});
p   =   p.Results;
lowmem = p.lowmem;
clear p


if ~exist('coils','var')
        coils = [];
end
    if ~isempty(coils)
        E   =   xfm_NUFFT_VEASL([imsize(1) imsize(2) imsize(3) size(kdata,2),...
            size(kdata,3)],coils,[],kspace,H, 'lowmem',lowmem);
    else


        % Estimate sensitivity maps
        E0 = xfm_NUFFT_VEASL([imsize(1) imsize(2) imsize(3) size(kdata,2) 1, 1],ones([imsize(1) imsize(2) imsize(3)]),[],...
            kspace(:,:,1,:),1,'lowmem',lowmem);
        
        T = calcToeplitzEmbedding(E0);
        ims = zeros([imsize(1) imsize(2) imsize(3)  size(kdata,4)],'single');

        for c = 1:size(kdata,4)
            disp(['Calculating coil sensitivity map ' num2str(c) ' of ' num2str(size(kdata,4))])
            kdata2 = single(E0'.*((mean(kdata(:,:,:,c),3).*E0.w).*E0.norm));
            ims(:,:,:,c) =   mean(fista(kdata2, E0, 1, 0, ...
                [E0.Nd, size(kdata,2),1], 10, 0.5,T),4);
        end
        disp(['Adaptive combine estimation of coil sensitivities'])
        sens = adaptive_estimate_sens('data', permute(ims,[4,1,2,3]),...
            'kernel', 10, 'thresh', 0.05);


        % Create E
        E   =   xfm_NUFFT_VEASL([imsize(1) imsize(2) imsize(3) size(kdata,2),...
            size(kdata,3)],sens,[],kspace,H,'lowmem',lowmem);
    end
end



