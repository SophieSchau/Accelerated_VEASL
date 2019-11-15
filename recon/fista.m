function [u, cost] = fista(d, xfm, W, lambda, im_size, niter, step, T)

%   Mark Chiew  
%   May 2017
%
%   Implementation of the FISTA iteration scheme for L1-Regularised problems
%   Adapted from Beck and Teboulle, SIAM J Im Sci 2009
%
%   This will solve a problem of the form:
%
%   min{u} |xfm*u - d|_2 + lambda*|W*u|_1
%
%   where u is the estimate,  W is an invertible  transform (often Wavelet), 
%   or identity transform
%   xfm is the measurement operator, lambda is the L1-weighting
%   and d is the measured raw k-space data

%   To make things easier, we work with our estimate u in the main loop
%   in the sparse domain, so that our measurement operator is actually 
%   effectively xfm*W'
%   Then we simply transform by W' as a final step
%   
%   Another way of viewing this is that we perform a change of variables
%   such that u_new = W*u
%
%   Constant stepsize variant
%
% Edited by Sophie Schauman Nov 2018

%   Initialise
    u   =   single(W*zeros(im_size));
    v   =   u;
    t   =   1;
    up  =   u;
    cost = zeros(niter,1);

    
    if nargin < 8
        d2 = single(xfm'.*d);
        clear d;
        T = calcToeplitzEmbedding(xfm);
    else
        d2 = d;
        clear d;
    end

%   Main loop
for ii = 1:niter

    %   Data consistency
%     u   =   v + step*(W*(xfm'.*(d - xfm*(W'*v))));
    u   =   v + step* (W  * ( d2 -   mtimes_Toeplitz(xfm, T, W'*v)));
       
    %   Solve proximal sub-problem
    u   =   shrink(u, lambda*step);

    %   Compute momentum parameter
    t2  =   (1+sqrt(1+4*t^2))/2;

    %   Compute momentum update
    v   =   u + ((t-1)/t2)*(u-up);
%    v = u; %ista

    
    %   Update variables
    t   =   t2;
    up  =   u;
    
    
% if mod(ii,10) == 0   
%     diff = d - xfm*(W'*u);
%     c1 = 0.5*norm(diff(:))^2;
%     c2 = lambda * norm(u(:),1);
%     cost(ii) = c1+c2;
%     fprintf(['iter: ', num2str(ii), '   cost: ', num2str(cost(ii)), '\n']);
%  end
    fprintf(['iter: ', num2str(ii), '\n']);
%     
%     imshow(abs(u(:,:,1,1,1)),[]);
%     drawnow
    
end

u   =   W'*u;
end


function y = shrink(x, thresh)
    y = exp(1j*angle(x)).*max(abs(x)-thresh,0);
end
