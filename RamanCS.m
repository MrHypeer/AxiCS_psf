function [recSpectra, waveNbCrop] = RamanCS(spectra, waveNb, ratio, method, PSF, cutoff, doCrop)
% Function takes a sub-sampled Raman measurement in,
% arranges the data points in a grid of the same size as the final resolution, 
% and performs inversion through compressed sensing to retrieve the full image.

%Input:
%spectra        original lower-resolution image (assumes line sub-sampling)
%               dimension are assumed to be [spectral, x, sub-sampled y]
%waveNb         wave number axis
%ratio          sub-sampling ratio
%method         reconstruction method (1: spatial, 2: PSF included)
%PSF            point-spread function of the imaging system (required if method==2)
%               should have same sampling as final image
%cutoff         cut-off value for truncated deconvolution (minimal value for division in spectral domain)
%doCrop         region of spectrum to use (1: fingerprint, 2: lipids, 3: full without silent region)

%Output
%recSpectra     reconstructed spectrum
%waveNbCrop     cropped wave number axis depending on region used

%Notes:
% Assumes a line subsampling. The ratio is given as the sub-sampling factor
% Example: a 512x1024 image to reach a 1024x1024 is a ratio 2.
%
% Copyright (C) 2015-2017  Nicolas Pavillon, Osaka University

path(path,genpath('./TVAL3'));
if (ratio < 1)
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

%% First cutting silent region
if (doCrop)
    silentLb = 1800;
    silentUb = 2700;
    [~,lb] = min(abs(waveNb-silentLb));
    [~,ub] = min(abs(waveNb-silentUb));

	switch (doCrop)
        case 1 %Fingerprint
        waveNbCrop = waveNb(1:lb);
        dataCrop = spectra(1:lb, :, :);
        case 2 %Lipids
        waveNbCrop = waveNb(ub:end);
        dataCrop = spectra(ub:end, :, :);
        case 3 %Fingerprint + Lipids
        waveNbCrop = waveNb([1:lb, ub:end]);
        dataCrop = spectra([1:lb, ub:end], :, :);
	end
end      

%% Preparing indices in which to put data
s = size(dataCrop);
%Preparing coordinates in the direction of lines
lineVal = 2:ratio:ratio*s(3);
%Preparing coordinates in other dimensions

%First dimension s(1): This one does not repeat values and is is repeated s(2)*s(3) times
sVal = 1:s(1); 
sVal = repmat(sVal, 1, s(2)*s(3));
%Second dimension s(2): This one repeats values s(1) times and is repeated s(3) times
xVal = 1:s(2); 
xVal = repmat(xVal, s(1), 1);
xVal = xVal(:)';
xVal = repmat(xVal, 1, s(3));
%Third dimension s(3): This one repeats values s(1)*s(2) times and is not repeated
colVal = repmat(lineVal, s(1)*s(2), 1);
colVal = colVal(:)';

%Preparing shuffling values
shuffleVal = ceil(rand(length(colVal), 1)*3) - 2;
newIdx = sub2ind([s(1) s(2) s(3)*ratio], sVal, xVal, colVal+shuffleVal');

%% Preparing sparse data
sparseData = zeros(s(1), s(2), ratio*s(3));
%Assigning data
sparseData(newIdx) = dataCrop(:);

%Finally, randomize indices
newIdx = newIdx(randperm(length(newIdx)));

%Clean up 
clear colVal lineVal shuffleVal xVal

%% Preparing CS algorithm
p = size(sparseData, 2); q = size(sparseData, 3); % p x q is the size of image
r = size(sparseData, 1); %Spectral dimension
N = p*q*r;

%Measurement matrix
switch(method)
    case 1
        A = @ (z, mode) dfA_sf(z, newIdx, N, mode);
	case 2
        dim = [r, p, q];
        [OTF, mask] = makeOTF(PSF, dim(2:3), cutoff);
        
    	A = @ (z, mode) dfA_fPsfD(z, dim, newIdx, OTF, mask, mode);  
end

%Measurement, modelled by measurement matrix
f = A(sparseData(:), 1);

%Cleaning as much memory as possible before running algorithm
clear spectra dataCrop sparseData sVal waveNb mask

%Preparing parameters
clear opts
opts.mu = 2^7; %2^8 by default, but adjusting for noisier Raman
opts.beta = 2^5;
opts.tol = 1E-3;
opts.maxit = 50; %Reducing maximum of iterations to speed up reconstruction
opts.TVnorm = 1;
opts.nonneg = false; %Spectra are background-subtracted, so negative is possible
opts.isreal = true;

recSpectra = TVAL3Spectral(A,f,r,p,q,opts);

%End of RamanCS

function y = dfA_sf(x,picks,N,mode)
switch mode
    case 1
        y = A_sf(x,picks);
    case 2
        y = At_sf(x, N, picks);
    otherwise
        error('Unknown mode passed to f_handleA!');
end

function y = dfA_fPsfD(x,dim,picks,OTF,mask, mode)
switch mode
    case 1
        y = A_fPsfF(x, dim, picks, OTF);
    case 2
        y = At_fPsfD(x, dim, picks, OTF, mask);
    otherwise
        error('Unknown mode passed to f_handleA!');
end
