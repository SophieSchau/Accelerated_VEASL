function [ kspace ] = gen_radial_traj( Phi, NSamps, Theta )
%GEN_RADIAL_TRAJ generates a list of kspace positions sampled by a radial
%trajectory
%
%   Sophie Schauman July 2018 - sophie.schauman@dtc.ox.ac.uk
%   Based on code by Thomas W Okell.
%
%   INPUTS: 
%   Phi           - list of angles (to y axis) to sample along (2D)
%   NSamps        - number of samples along each spoke
%   (Theta        - additional list of angles (to z-axis) to sample along (3D))
%
%   OUTPUTS:
%   kspace        - 2D matrix 
%                  (nKPoints x nSpatialDimensions)
%
%%%

if isempty(Theta) % 2D
    NSpokes = length(Phi);
    dr = 2*pi/NSamps;
    r   =   (-pi:dr:(NSamps-1)*dr-pi)';
    
    % Initialise the trajectory
    kspace = zeros(NSpokes*NSamps,2); % radians/pixel
    
    for ii = 1:NSpokes
        Idx = (ii-1)*NSamps+1:ii*NSamps;
        kspace(Idx,:) = r * [sin(Phi(ii)) cos(Phi(ii))];
    end
    
    
else % 3D
    NSpokes = length(Phi);
    dr = 2*pi/NSamps;
    r   =   (-pi:dr:(NSamps-1)*dr-pi)';
    
    % Initialise the trajectory
    kspace = zeros(NSpokes*NSamps,3); % radians/pixel
        
    
    for ii = 1:NSpokes
        Idx = (ii-1)*NSamps+1:ii*NSamps;
        kspace(Idx,:) = r * [sin(Phi(ii))*sin(Theta(ii)) cos(Phi(ii))*sin(Theta(ii)) cos(Theta(ii))];
    end
    
    
end



end

