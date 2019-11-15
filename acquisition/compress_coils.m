function [kdata_cc] = compress_coils(kdata, c)
%COMPRESS_COILS Performs coil compression using singular value
%decomposition. (Buehrer, MRM 2007)
%
%   Sophie Schauman Oct 2018 
%
%   INPUTS: 
%   kdata         - raw data from scanner (extracted from twix file)
%                   4D matrix (nKPoints x nTPoints x  nEncodings x nCoils)
%   c             - number of virtual coils to create
%
%   OUTPUTS:
%   kdata_cc      - new compressed k-data. 4D matrix.
%                   (nKPoints x nTPoints x  nEncodings x nVirtualCoils)
%
%   DEPENDENCIES:
%   lsvd.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[u,s,~] = lsvd(reshape(kdata,[],size(kdata,4)),c);
kdata_cc = reshape(u*s, size(kdata,1), size(kdata,2), size(kdata,3), []);

end
