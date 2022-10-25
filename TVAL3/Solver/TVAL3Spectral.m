% TVAL3Spectral.m
% 
% Function to perform reconstruction based on L1 minimization for 
% 3D matrices based on the model:
%   min sum ||D_i u||, s.t. Au = b
%
% Input:
% A             measurement matrix (i.e. structure of measurement)
%               defined either as an explicit matrix (greedy on RAM),
%               or by an implicit function
% b             measurement vector (i.e. data, 1D)
% p, q, r       dimensions of 3D output data
% opts          
%
%
% Output:
% normVal       norm
%
% Copyright (C) 2015-2017  Nicolas Pavillon, Osaka University

function [U,out] = TVAL3Spectral(A,b,p,q,r,opts)

setVersion(1);

if ~isfield(opts,'TVL2')
    opts.TVL2 = false;
end
if opts.TVL2
    warning('L2 model is not implemented for spectral reconstructions');
end

[U, out] = ftvcs_alp3(A,b,p,q,r,opts);

