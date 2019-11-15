% Returns azimuthal and polar angles for N 3D golden means radial spokes as per Chan et
% al. MRM 2009.  If N is an array, this is taken as the golden ratio index
% array.

function [Azi, Polar] = GoldenMeans3D(N,ReverseOddLines)

if nargin < 2; ReverseOddLines = false; end

% Define the increments, from Chan et al
Phis = GoldenRatios3D;

% Calculate Polar and Azimuthal angles
% NB. apparent error in Chan paper figure here - Beta is the angle from the
% kz axis = Polar angle in Siemens terms
if length(N) == 1
    m = (0:(N-1))';
else
    m = N(:);
end

kz = mod(m*Phis(1),1); % Use GR to evenly distribute between 0 and 1
% Can potentially invert the sign and add pi to every other azimuthal angle
% to get samples in both directions if desired, but leave this until
% later...

Polar = acos(kz); 
Azi   = mod(m*Phis(2),1) * 2 * pi;

% Reverse every other line if requested
if ReverseOddLines
    OddIdx = logical(mod(m,2));
    
    % Add pi to the azimuthal angle
    Azi(OddIdx) = mod( Azi(OddIdx) + pi, 2*pi);
    
    % Reverse kz
    Polar(OddIdx) = acos(-kz(OddIdx)); 
end

%% Old version
% % Define the increments, from Chan et al (??may need more precision here??)
% Phis = [0.4656 0.6823];
% 
% % Calculate mPhis modulo 1
% mPhis = mod(((1:N)')*Phis,1);
% 
% % Convert to x,y,z
% Alphas = 2*pi*mPhis(:,2);
% %Betas  = asin(mPhis(:,1)); % NB. apparent error in Chan paper here.
% 
% z = mPhis(:,1);
% xysize = sqrt(1^2 - z.^2);
% x = xysize .* cos(Alphas);
% y = xysize .* sin(Alphas);
% 
% UnitVecs = [x y z];

