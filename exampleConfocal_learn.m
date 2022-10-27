%exampleConfocal.m
%
% Example script showing the use of the ConfocalCS function based on demo
% data.
% See README.txt about where to get the demo data (not included with code)
%
% Copyright (C) 2015-2017  Nicolas Pavillon, Osaka University

clear all
%% Choose reconstruction type: 1 (spatial) or 2 (PSF included)
recType = 2;

%% Retrieve data
imageName = 'im2'; %other: im2

data_folder = '../data_psf/Confocal_Data/';
im = imread([data_folder, imageName, '-256.tif']);
imFull = imread([data_folder, imageName, '-1024.tif']);
if (recType == 2)
    psfIm = imread(strcat(data_folder,'5-PSF_avg_centered.tif'));
    cutoff = 0.15;      %Cut-off     
else
    PSF = [];
    cutoff = 0;
end
ratio = 4;          %Sub-sampling ratio

%% Adjusting PSF size to measurement wavelength
if (recType == 2)
    %Pixel size for PSF measured by beads
    dx = 20.716; %nm
    %Emission of beads
    lambda = 567; %nm
    
    %Pixel size of actual images (cells, ER-tracker)
    dx2 = 172.63; %nm, px size for ER-tracker with zoom 1.2x
    lambda2 = 511; %nm, ER-tracker
    
    %Scaling the pixel size with the wavelength
    % why pixel size is related with wavelength???
    dx3 = dx2*lambda/lambda2;
    
    %Create the xy map (assuming even-sized square image)
    l = floor(size(psfIm)./2);
    x = -(l(1)-1)*dx:dx:l(1)*dx; %x sampling of measured PSF
    xi = 0:dx3:7*dx3; %sampling of PSF in object space (15x15 matrix)
    xi = [-xi(end:-1:2), xi];

    %Interpolation in object space
    [xx, yy] = meshgrid(x, x);
    [xxi, yyi] = meshgrid(xi, xi);
    PSF = interp2(xx, yy, psfIm, xxi, yyi);
    %Force border ones (outside measured range if any) to be at background level
    nanIdx = find(isnan(PSF(:)));
    PSF(nanIdx) = median(psfIm(:));
end

tic
[rec, samples] = ConfocalCS(im, ratio, recType, PSF, cutoff);
t = toc;

subplot(221), imagesc(im), title('Measured low-res image')
subplot(222), imagesc(samples), title('Sampling grid');
subplot(223), imagesc(rec), title('Reconstruction')
subplot(224), imagesc(imFull), title('Measured full-res image (comparison)');
disp(['Reconstruction time: ', num2str(t), ' s']);

