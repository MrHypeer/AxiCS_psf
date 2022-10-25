%exampleRaman.m
%
% Example script showing the use of the RamanCS function based on demo
% data.
% See README.txt about where to get the demo data (not included with code)
%
% Copyright (C) 2015-2017  Nicolas Pavillon, Osaka University

%% Choose reconstruction type: 1 (spatial) or 2 (PSF included) and region
doCrop = 2; %1 (only fingerprint), 2 (only lipids), 3 (both)
recType = 2; %Deconvolution case

disp('Warning: 3D reconstruction can use a large amount of memory.');

%% Retrieve data
load('Data/Raman_Data/3-3_bkg.mat');

cutoff = 0.15;
ratio = 3;
%Reading PSF
psfIm = imread('Data/Raman_Data/PSF_avg.tif');
    
%Pixel size in PSF domain
dx = 103.5; %nm
%Emission of beads
lambda = 561; %nm, measured as average in 900-1000 cm-1

%Pixel size of actual images (binned images)
dx2 = 207; %nm

%Create the xy map (assuming even-sized square image)
l = floor(size(psfIm)./2);
x = -(l(1)-1)*dx:dx:l(1)*dx; %x sampling of measured PSF
xi = 0:dx2:7*dx2; %sampling of PSF (15x15 matrix)
xi = [-xi(end:-1:2), xi];

%Interpolation in object space
[xx, yy] = meshgrid(x, x);
[xxi, yyi] = meshgrid(xi, xi);
PSF = interp2(xx, yy, psfIm, xxi, yyi);

%Force border ones (outside measured range if any) to be at background level
nanIdx = find(isnan(PSF(:)));
PSF(nanIdx) = median(psfIm(:));

tic
if (doCrop < 3)
    [rec, waveNbRec] = RamanCS(spectra, waveNb, ratio, recType, PSF, cutoff, doCrop);
else
    %Two reconstructions that we concatenate
    [rec1, waveNbRec1] = RamanCS(spectra, waveNb, ratio, recType, PSF, cutoff, 1);
    [rec2, waveNbRec2] = RamanCS(spectra, waveNb, ratio, recType, PSF, cutoff, 2);
    rec = cat(1, rec1, rec2);
    waveNbRec = [waveNbRec1, waveNbRec2];
    clear rec1 rec2
end
t = toc;

switch(doCrop)
    case 1
        lb = 1;
        [~,ub] = min(abs(waveNb-1800));
        w = 1659;
    case 2
        ub = length(waveNb);
        [~,lb] = min(abs(waveNb-2700));
        w = 2930;
    case 3
        [~,lb] = min(abs(waveNb-1800));
        [~,ub] = min(abs(waveNb-2700));
        w = 2930;
end

[~, peak] = min(abs(waveNb-w));
subplot(221), imagesc(squeeze(spectra(peak, :, :))'), title(['Measured low-res (', int2str(w),' cm-1)']);
if (doCrop < 3)
    subplot(222), plot(waveNb(lb:ub), spectra(lb:ub, 146, 70)), title('Single-pixel measured spectrum');
else
    subplot(222), plot(waveNbRec, spectra([1:lb, ub:end], 146, 70)), title('Single-pixel measured spectrum');
end
[~, peak2] = min(abs(waveNbRec-w));
subplot(223), imagesc(squeeze(rec(peak2, :, :))'), title(['Reconstruction (', int2str(w), ' cm-1)']);
subplot(224), plot(waveNbRec, rec(:, 146, 208)), title('Single-pixel reconstructed spectrum');
disp(['Reconstruction time: ', num2str(t), ' s']);
