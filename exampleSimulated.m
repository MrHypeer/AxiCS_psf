%exampleSimulated.m
%
% Example script showing the use of spatial sub-sampling with TVAL3 through
% simulated data.
%
% The script is given as a function to allow the definition of internal
% functions.
%
% Copyright (C) 2015-2017  Nicolas Pavillon, Osaka University

function exampleSimulated()

addpath('./Simul');
path(path,genpath('./TVAL3'));

%% Set general parameters
method = 4;     %Choose reconstruction type: 1 (Fourier) 2 (Hadamard), 
                %3 (sub-spatial), 4 (PSF included)
doNoise = 1;    %Add noise
doSmooth = 1;   %Use smoothing (i.e. PSF)
FWHM = 5;       %Full width at half maximum of smoothing Gaussian
doShow = 1;     %Activates display of images
signalStrength = 500; %Amplitude of signal, influences amount of noise
ratio = 4;

%% Set random generator (for reproducible results)
%Comment out for different random realizations
if (~exist('randSeed.mat', 'file'))
    s = rng;
    save('randSeed.mat',  's');
else
    load randSeed.mat
end
rng(s);

%% Prepares output text
switch(method)
    case 1
        methodText = 'Fourier ensemble';
    case 2
        methodText = 'Hadamard matrix';
    case 3
        methodText = 'Spatial sub-sampling';
    case 4
        methodText = 'Spatial with PSF included';
end
if (doNoise)
    noiseText = ['Poisson noise (max. ', int2str(signalStrength), ' photon per pixel)'];
else
    noiseText = 'no noise';
end
if (doSmooth)
    smoothText = ['Gaussian smoothing (FWHM=', int2str(FWHM), ')'];
else
    smoothText = 'no smoothing';
end
disp(['Testing reconstruction with ', methodText, ' and ', num2str(100/ratio^2, '%.2f'), '% of data.']);
disp(['   Using ', smoothText, ' and ', noiseText]);

%% Prepare pattern
pattern = createPattern();
%Normalizing
pattern = pattern-min(pattern(:));
pattern = pattern/max(pattern(:));
%Leaving 10% of amplitude for background
pattern = pattern*0.9*signalStrength + 0.1*signalStrength;

patternHard = pattern;
if(doSmooth)
    [x, y] = meshgrid(-3:3, -3:3);
    sigma = FWHM/(2*sqrt(2*log(2)));
    kernel = 1/(sigma*sqrt(2*pi)) * exp ( - (x.^2 + y.^2)./(2*sigma^2) );
    kernel = kernel./sum(kernel(:));
    pattern = imfilter(pattern, kernel);
else
    %Unit kernel
    kernel = zeros(3, 3);
    kernel(2, 2) = 1;
end
    
%Add noise if requested
noiseless = pattern;
if (doNoise ~= 0)
    %Poisson noise
    noise = PoissonFct(pattern(:)); %Poisson noise is not available in built-in Matlab
    pattern = reshape(noise, size(pattern, 1), size(pattern, 2));
end

%Normalizing 110% of amplitude to account for possible noise additions
pattern = pattern./(1.1*signalStrength);
patternHard = patternHard/(1.1*signalStrength);
noiseless = noiseless/(1.1*signalStrength);

%% Initializations and measurement matrix
p = size(pattern, 1); q = size(pattern, 2); % p x q is the full size of image
N = p*q;      %Full data size
M = round(1/ratio*N);
%Ensuring that M is even
if (M/2 ~= round(M/2))
    M = M+1;
end

% Random permutations and selection generation
P = randperm(N);
Q = randperm(N/2-1)+1;
freq = Q(1:M/2);       %freq is for Fourier. 2x less samples because of frequency redundancy (1 sample is 2 points, real and imag)
picks = P(1:M);

%Measurement matrix
switch(method)
    case 1
        A = @(z, mode) dfA_f(z, freq, N, P, mode);
    case 2
        A = @(z, mode) dfA_fWH(z, picks, P, mode);
    case 3
        A = @ (z, mode) dfA_sf(z, picks, N, mode);
    case 4
        dim = [p, q];
        cutoff = 0.20;
        [OTF, mask] = makeOTF(kernel, dim, cutoff);
        
        A = @ (z, mode) dfA_fPsfD(z, dim, picks, OTF, mask, mode);
end

%% Reconstruction
% Original image
I = pattern;

% Observation, modelled by measurement matrix
f = A(I(:), 1);

% First guess, based on inversion of observations
x0 = A(f, 2);
I0 = reshape(x0, p, q);

%% Run TVAL3
clear opts
opts.mu = 2^8; opts.beta = 2^5;
opts.tol = 1E-3;
opts.maxit = 300;
opts.TVnorm = 1;
opts.nonneg = true;
opts.isreal = true;
opts.disp = false;

tic;
U = TVAL3(A,f,p,q,opts);
t = toc;

%Computing rms error (recovered vs original)
rms = sqrt(1/N * sum(sum((U-I).^2)) );
disp(['Error to original is: ', num2str(rms, '%.2f')]);
snr = 10*log10( sum(sum(U.^2)) ./ sum(sum( (U-I).^2 )) );
disp(['SNR of image (vs original): ', num2str(snr, '%.2f')]);

if (doShow) %Showing results
    subplot(3, 2, 1), imagesc(I), title('Original image'), colorbar();
    subplot(3, 2, 2), imagesc(I0), title('Initial guess'), colorbar();
    subplot(3, 2, 3), imagesc(U), title('Recovered signal'), colorbar();
    subplot(3, 2, 4), imagesc(U-I), title('Difference'), colorbar();
    if (doNoise)
        subplot(325), imagesc(U-noiseless), title('Difference to noiseless'), colorbar();
    else
        subplot(325), plot(0, 0), colorbar();
    end
    if (doSmooth)
        subplot(326), imagesc(U-patternHard), title('Difference to non-smoothed'), colorbar();
    else
        subplot(326), plot(0, 0), colorbar();
    end
    disp(['Execution time is: ', num2str(t, '%.3f'), ' s']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fourier ensemble
function y = dfA_f(x,freq,N,P,mode)
switch mode
    case 1
        y = A_fFT(x, freq, P);
    case 2
        y = At_fFT(x, N, freq, P);
    otherwise
        error('Unknown mode passed to f_handleA!');
end

%Hadamard
function y = dfA_fWH(x,freq,P,mode)
switch mode
    case 1
        y = A_fWH(x, freq, P);
    case 2
        y = At_fWH(x, freq, P);
    otherwise
        error('Unknown mode passed to f_handleA!');
end

%Spatial subsampling
function y = dfA_sf(x,picks,N,mode)
switch mode
    case 1
        y = A_sf(x,picks);
    case 2
        y = At_sf(x, N, picks);
    otherwise
        error('Unknown mode passed to f_handleA!');
end

%Deconvolution
function y = dfA_fPsfD(x,dim,picks,OTF,mask, mode)
switch mode
    case 1
        y = A_fPsfF(x, dim, picks, OTF);
    case 2
        y = At_fPsfD(x, dim, picks, OTF, mask);
    otherwise
        error('Unknown mode passed to f_handleA!');
end
