function [recImage, sparseData] = AxiCS(image, ratio, method, PSF, cutoff)
% Function takes a sub-sampled image measurement in, 
% arranges the data points in a grid of the same size as the final resolution, 
% and performs inversion through compressed sensing to retrieve the full image.

%Input:
%image          original lower-resolution image
%ratio          sub-sampling ratio
%method         reconstruction method (1: spatial, 2: PSF included)
%PSF            point-spread function of the imaging system (required if method==2)
%               should have same sampling as final image
%cutoff         cut-off value for truncated deconvolution (minimal value for division in spectral domain)

%Output
%recImage       reconstructed image
%sparseData     original sample points in grid of final size

%Notes:
% Assumes a square subsampling. The ratio is given as the sub-sampling factor
% Example: a 512x512 image to reach a 1024x1024 is a ratio 2.
%
% Copyright (C) 2015-2017  Nicolas Pavillon, Osaka University

path(path,genpath('./TVAL3'));
if (any(ratio < 1)) % modify to suit 2d ratio
    error('ratio must be > 1.');
end
if (method ~= 1 && method ~= 2)
    error('method must be 1 or 2');
end
if (method == 2)
    if (~exist('PSF', 'var') || isempty(PSF))
        error('PSF array is required with method 2')
    end
    if (~exist('cutoff', 'var') || isempty(cutoff))
        %Default value (rather high to ensure convergence)
        cutoff = 0.2;
    end
end


%% Preparing indices in which to put data
% this part was modified to suit square FOV from DMD
% mode to suit uneven sparisity along x-y direction@20221103
s = size(image); 
[subX, subY] = meshgrid( floor(ratio(2)/2)+1:ratio(2):ratio(2)*s(2), ...
                         floor(ratio(1)/2)+1:ratio(1):ratio(1)*s(1)); 
                     
%Computing indices, sparse raw_data to reconstruction dimension
newIdx = sub2ind(ratio.*s, subY(:), subX(:)); 

%% Preparing sparse data
sparseData = zeros(ratio.*s);
sparseData(newIdx) = image(:);

%Finally, randomize indices
newIdx = newIdx(randperm(length(newIdx)));

%% Preparing parameters of CS algorithm
p = size(sparseData, 1); q = size(sparseData, 2); % p x q is the size of image
N = p*q;

%Measurement matrix
switch(method)
    case 1 %Spatial
        A = @ (z, mode) dfA_sf(z, newIdx, N, mode);
    case 2 %PSF included
        %Prepare the OTF mask
        [OTF, mask] = makeOTF(PSF, [p, q], cutoff);
    	A = @ (z, mode) dfA_fPsfD(z, [p, q], newIdx, OTF, mask, mode);              
end
    
%Measurement, modelled by measurement matrix
f = A(sparseData(:), 1); % original is 1

%% Preparing parameters
clear opts
opts.mu = 2^4;
opts.beta = 2^4;
%To push quality in case of low-noise signal
% opts.mu = 2^12; opts.beta = 2^9;
opts.tol = 1E-3;
opts.maxit = 300;
opts.TVnorm = 1;
opts.nonneg = true;
opts.isreal = true;

%% Reconstruction itself
recImage = TVAL3(A,f,p,q,opts);

%End of ConfocalCS

%% Functions
%Spatial implicit functions
function y = dfA_sf(x,picks,N,mode)
switch mode
    case 1
        y = A_sf(x,picks);
    case 2
        y = At_sf(x, N, picks);
    otherwise
        error('Unknown mode passed to f_handleA!');
end

%PSF-included implicit functions
function y = dfA_fPsfD(x,dim,picks,OTF,mask, mode)
switch mode
    case 1
        y = A_fPsfF(x, dim, picks, OTF);
    case 2
        y = At_fPsfD(x, dim, picks, OTF, mask);
    otherwise
        error('Unknown mode passed to f_handleA!');
end
