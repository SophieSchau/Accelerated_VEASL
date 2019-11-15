function [u, cost] = fista_L2s(d, xfm, W, lambda1, lambda2, im_size, niter, step)

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
u   =   W*zeros(im_size, 'single');
v   =   W*(xfm'.*d);
t   =   1;
up  =   u;
cost = zeros(niter,1,'single');


d2 = W'*v;

T = calcToeplitzEmbedding(xfm);

z = mtimes_Toeplitz(xfm, T, W'*v);


fprintf(1, '%-5s %-16s %-16s %-16s %-16s\n', 'Iter','DataCon','L1','L2','Cost');

%   Main loop
for ii = 1:niter

    %   Data consistency
%     u   =   v + step*(W*(xfm'.*(d - xfm*(W'*v))));
    u   =   v + step* (W  * ( d2 -   z - lambda2*R1(W'*v)));
       
    %   Solve proximal sub-problem
    u   =   shrink(u, lambda1*step);

    %   Compute momentum parameter
    t2  =   (1+sqrt(1+4*t^2))/2;

    %   Compute momentum update
    v   =   u + ((t-1)/t2)*(u-up);
%    v = u; %ista

    
    %   Update variables
    t   =   t2;
    up  =   u;
    z   =   mtimes_Toeplitz(xfm, T, W'*v);
    
    
    %   Error terms and cost function
    err1(ii)  =   (W'*u(:))'*(z(:)-2*d2(:))+d(:)'*d(:);
    err2(ii)  =   lambda1*norm(W'*u(:),1);
    err3(ii)  =   lambda2*sum(abs(reshape(R1(W'*u),[],1)).^2);
    cost(ii)  =   0.5*err1(ii) + err2(ii) + 0.5*err3(ii);

    %   Display iteration summary data
    fprintf(1, '%-5d %-16G %-16G %-16G %-16G\n', ii, err1(ii), err2(ii), err3(ii), cost(ii));
%     imagesc(squeeze(max(abs(mean(u(:,:,:,:,1),4)),[],2)))
%     drawnow
    
%     imshow(abs(u(:,:,1,1,1)),[]);
%     drawnow
    
end

u   =   W'*u;
end


function y = shrink(x, thresh)
    y = exp(1j*angle(x)).*max(abs(x)-thresh,0);
end

function x = R1(x)
    %x =  -1*circshift(x,-1,2) + 2*x - 1*circshift(x,1,2);
    if size(x,4) > 1
        temp = -1*circshift(x,-1,4) + 2*x - 1*circshift(x,1,4);
        x = cat(4, x(:,:,:,1,:) - x(:,:,:,2,:), temp(:,:,:,2:end-1,:), x(:,:,:,end,:)-x(:,:,:,end-1,:));
    else
        x = 0;

    end
end