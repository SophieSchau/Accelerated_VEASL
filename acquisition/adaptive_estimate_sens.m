function varargout = adaptive_estimate_sens(varargin);
%
%
% sens = adaptive_estimate_sens(data, kernel)
%
%   Mark Chiew
%   May 2014
%   Last Updated: Oct 2014
%
%   Estimate coil sensitivities using Adaptive Combine (Walsh, MRM 2000)
%
%   data is [Nc x Nx x Ny x Nz] complex data array
%   data can also be a string, in which case complex array is generated 
%   kernel is a radius 
%
%   n_parts is an optional parameter for parallelisation, how many partial 
%   problems to split this into
%   part is the current part to be worked on
%
%   sens is a [Nc, Nx, Ny, Nz] complex sensitivity array
%   Phases are anchored to the sensitivity of the first channel
%   and sensitivities are normalised

p   =   inputParser;

p.addParameter('data',  []);
p.addParameter('dir',   [], @(x) exist(x)==7);
p.addParameter('kernel',4,  @isscalar);
p.addParameter('parts', 1,  @isscalar);
p.addParameter('part',  1,  @isscalar);
p.addParameter('z',     0,  @isscalar);
p.addParameter('thresh',0,  @isscalar);

p.parse(varargin{:});

data    =   double(p.Results.data);
data_dir=   p.Results.dir;
kernel  =   p.Results.kernel;
n_parts =   p.Results.parts;
part    =   p.Results.part;
z       =   p.Results.z;
thresh  =   p.Results.thresh;

clear p;

if isempty(data) && isempty(data_dir)
    error('No valid data specified');
end

if isempty(data) && ~isempty(data_dir)
    files=dir([data_dir '/tmean_c*']);
    load([data_dir '/' files(1).name()])
    data=zeros([length(files) size(est)]);
    data(1,:,:,:)=est;
    for i = 2:size(data,1)
        load([data_dir '/' files(i).name()]);
        data(i,:,:,:)=est;
    end
end

if z
    data    =   data(:,:,:,z);
end

dims    =   size(data);
ncoils  =   dims(1);

N   =   prod(dims(2:end));
n   =   ceil(N/n_parts);
idx =   (part-1)*n+1:min(part*n, N); 

sens    =   zeros(size(data));
mask    =   zeros(size(data(1,:,:,:)));

for i = idx
    [x,y,z] =   ind2sub(dims(2:end),i);
    ii      =   getROIidx(dims(2:end), [x,y,z], kernel);
    R       =   data(:,ii)*data(:,ii)';
    [V,D]   =   eigs(R,1);

    %   Pin coil-phase to first channel
    sens(:,i)   =   V*exp(-1j*angle(V(1)));
    mask(i)     =   sqrt(D);
end

if nargout == 1
    varargout{1}    =   permute(sens.*(abs(mask)>thresh*max(abs(mask(:)))),[2,3,4,1]);
elseif nargout == 2
    varargout{1}    =   sens(:,idx);
    varargout{2}    =   idx;
elseif nargout == 3
    varargout{1}    =   permute(sens,[2,3,4,1]);
    varargout{2}    =   [];
    varargout{3}    =   mask;
end



function idx = getROIidx(dims, centre, r)

if length(dims)==2
    centre  =   centre(1:2);
    [sx, sy]    =   ndgrid(-r:r,-r:r);
    xy = bsxfun(@plus, centre, [sx(:) sy(:)]);
    xy = bsxfun(@min, max(xy,1), dims); 
    xy = unique(xy,'rows');                     

    test    =   sum((xy-repmat(centre,size(xy,1),1)).^2,2).^0.5;
    roi     =   find(test <= r);

    idx =   sub2ind(dims, xy(roi,1), xy(roi,2));
elseif length(dims) == 3
    [sx, sy, sz]    =   ndgrid(-r:r,-r:r,-r:r);
    xyz = bsxfun(@plus, centre, [sx(:) sy(:) sz(:)]);
    xyz = bsxfun(@min, max(xyz,1), dims); 
    xyz = unique(xyz,'rows');                     

    test    =   sum((xyz-repmat(centre,size(xyz,1),1)).^2,2).^0.5;
    roi     =   find(test <= r);

    idx =   sub2ind(dims, xyz(roi,1), xyz(roi,2), xyz(roi,3));
end
