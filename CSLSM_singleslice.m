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
imageName = 'wfsample1'; %other: im2

data_folder = '../data_psf/Confocal_Data/';
im = imread([data_folder, imageName, '-512.tif']);
imFull = imread([data_folder, imageName, '-1024.tif']);
if (recType == 2)
    psfIm = double(imread(strcat(data_folder,'5-PsfGen4.tif')));
    cutoff = 0.15;      %Cut-off     
else
    PSF = [];
    cutoff = 0;
end
ratio = [2 2];          %Sub-sampling ratio % change to x-y ratio

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

figure(1)
subplot(221), imagesc(im), title('Measured low-res image')
subplot(222), imagesc(samples), title('Sampling grid');
subplot(223), imagesc(rec), title('Reconstruction')
subplot(224), imagesc(imFull), title('Measured full-res image (comparison)');
disp(['Reconstruction time: ', num2str(t), ' s']);

%% evaluate the reconstruction performance
% normalize two contrasting objects
idx = randi(size(imFull,1));
rec_norm = mapminmax(rec,0,1);
imFull_norm = mapminmax(imFull,0,1);
% plot intensity profile along sample
figure(2)
plot(1:size(rec,2), rec_norm(idx,:), 1:size(imFull,2), circshift(imFull_norm(idx,:),-2))
legend('reconstructed','full resolution')
title(strcat('line across x=',num2str(idx)))
% calculate PSNR
[psnr_value, snr_value] = psnr(rec_norm, imFull_norm)
RMSE = sqrt(immse(rec_norm, imFull_norm))




