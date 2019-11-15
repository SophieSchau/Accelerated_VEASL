classdef (Abstract) xfm
properties (SetAccess = protected, GetAccess = public)
    Nd      =   [];
    Nt      =   0;
    msize   =   []; 
    adjoint =   0;
    dsize   =   [];
    S       =   [];
    Nc      =   0;
    M;
end

methods
function res = xfm(dims, coils, fieldmap_struct)
    %   Initialise dimensions 
    res.Nd      =   dims(1:3);
    res.Nt      =   dims(4);
    res.msize   =   [prod(res.Nd) res.Nt];
    res.adjoint =   0;

    %   Initialise sensitivity encoding operator
    if isempty(coils)
        coils   =   ones([res.Nd 1]);
    end
        res.S   =   sensEncodingMatrix(coils);
        res.Nc  =   res.S.Nc;

    %   Initialise multi-frequency interpolation off-resonance correction operator
    if isempty(fieldmap_struct)
        res.M   =   mfi();
    else
        for i = 1:length(fieldmap_struct)
            M(i)    =   mfi(fieldmap_struct(i).field, ...
                            fieldmap_struct(i).t, ...
                            fieldmap_struct(i).L, ...
                            fieldmap_struct(i).idx);
        end
        res.M   =   M;
    end
end

function res = ctranspose(a)
    
    a.adjoint   = xor(a.adjoint,1);
    res         = a;

end

function step = max_step(xfm)
    step    =   norm(ones(xfm.msize(1),1))/norm(mtimes(xfm',mtimes(xfm,ones(xfm.msize(1),1),1),1));
end

function est = cg(xfm, d, tol, iters)

    %   Performs iterative SENSE recon using built-in lsqr
    %   Input d should be shaped like the output of mtimes

    if nargin < 3
        tol     =   [];
    end
    if nargin < 4
        iters   =   100;
    end

    [est, flag, relres] =   lsqr(@(x,mode) afun(x,mode,xfm), reshape(d,[],1), tol, iters);
    est =   reshape(est, xfm.msize);
end

function est = tmean(xfm, d, tol, iters)

    %   Performs iterative SENSE recon using built-in lsqr
    %   Input d should be shaped like the output of mtimes

    if nargin < 3
        tol     =   [];
    end
    if nargin < 4
        iters   =   100;
    end

    [est, flag, relres] =   lsqr(@(x,mode) bfun(x,mode,xfm), reshape(d,[],1), tol, iters);
end

function est = iter(xfm, d, w, optfn, tol, iters, T, est0)

    %   Performs iterative recon using built-ins
    %   Input d should be shaped like the output of mtimes
    %   Solves normal equation

    if nargin < 3 || isempty(w)
        w       =   ones(xfm.Nc,1);
    end
    if nargin < 4
        optfn   =   @minres;
    end
    if nargin < 5
        tol     =   [];
    end
    if nargin < 6
        iters   =   100;
    end
    if nargin < 7
        T       =   [];
    end
    if nargin < 8
        est0    =   [];
    end

    d   =   xfm'*(d.*reshape(w,1,1,[]));
    [est, flag, relres, iter] =   optfn(@(x,mode) cfun(x,xfm,T,w), reshape(d,[],1), tol, iters, [], [], est0);
    est =   reshape(est, xfm.msize);
    
    if iter < iters
        disp(['Converged after ' num2str(iter) ' iterations'])
    else
        disp(['Max iterations (' num2str(iters) ') reached'])
    end

end

function res = times(a,b)
    if a.adjoint
        res =   reshape(mtimes(a,b), [a.Nd(1:2) a.Nd(3), a.Nt]);
    else
        res =   mtimes(a,b);
    end
end

end

methods (Static)
function b = fftfn(b,dims)
    d   =   sqrt(size(b));
    for i = intersect(1:ndims(b), dims)
        b   =   fftshift(fft(ifftshift(b, i), [], i), i)/d(i);
    end
end
function b = fftfn_ns(b,dims)
    d   =   sqrt(size(b));
    for i = intersect(1:ndims(b), dims)
        b   =   fft(b, [], i)/d(i);
    end
end

function b = ifftfn(b,dims)
    d   =   sqrt(size(b));
    for i = intersect(1:ndims(b), dims)
        b   =   fftshift(ifft(ifftshift(b, i), [], i), i)*d(i);
    end
end
function b = ifftfn_ns(b,dims)
    d   =   sqrt(size(b));
    for i = intersect(1:ndims(b), dims)
        b   =   ifft(b, [], i)*d(i);
    end
end
function x = size(b)
    [x(1) x(2) x(3) x(4)]   =   size(b);
end
end

end

function y = afun(x, mode, xfm)
    if strcmp(mode, 'transp')
        y   =   reshape(xfm'*reshape(x,xfm.dsize),[],1);
    else
        y   =   reshape(xfm*reshape(x,xfm.msize),[],1);
    end
end

function y = bfun(x, mode, xfm)
    if strcmp(mode, 'transp')
        y   =   mean(xfm'*reshape(x,xfm.dsize),2);
    else
        y   =   reshape(xfm*repmat(x,1,xfm.Nt),[],1);
    end
end

function y = cfun(x, xfm, T, w)
    y   =   reshape(mtimes_Toeplitz(xfm, T, x, w),[],1);
end
