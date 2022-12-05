%%%
% generate high_res and low_res using psf connnnnnvolution

clear all;
%% Step1: load original image, psf and convolution
target_folder = 'D:\OneDrive - City University of Hong Kong\CUHK\08 project\14 CS\data_psf\Confocal_Data\append\USAF\';
addpath(target_folder)
bmpName = 'USAF1951G7.tif';
psfname = 'PSF-5pxl-64.tif';
imgHighRes = double(imread(bmpName));
psf = imread(psfname);
noise_flag = 1; % noise switch
mean_gaussian = [0 0.1 0.2];
var_gaussian = [0.01 0.03 0.05 0.07];


%% Step2: generate file name
% change file index in target file
splitname_bmp = split(bmpName,'.');
splitname_psf = split(psfname,'.');

oldFolder = cd(target_folder);
order = 1;
name2save = strcat(cell2mat(splitname_bmp(1)),cell2mat(splitname_psf(1)),...
    '_noisestack',num2str(order));
name2save_tif = strcat(name2save,'.tif');
name2save_txt = strcat(name2save,'.txt');

% file alreaady exist, automatically +1
while ~isempty(dir(name2save_txt)) % means already file has same name
    order = 1 + order;
    name2save = strcat(cell2mat(splitname_bmp(1)),...
        cell2mat(splitname_psf(1)), '_noisestack',num2str(order));
    name2save_tif = strcat(name2save,'.tif');
    name2save_txt = strcat(name2save,'.txt');
end
cd(oldFolder)


%% Step3: scale psf into image space and make OTF
% wavelength information
lambda_psf = 511;
lambda_img = 511;
% pixel size information
[x_psf, y_psf] = size(psf);
dx_psf = 100; % nm, same as dz_psf
[x_imgHigh, y_imgHigh] = size(imgHighRes);
dx_imgHigh = 730; %unit, nm

%Create the xy map (assuming even-sized square image), low res
dx2 = dx_imgHigh*lambda_psf/lambda_img;
l = floor([x_psf, y_psf]./2);
x = -(l(1)-1)*dx_psf:dx_psf:l(1)*dx_psf; %x sampling of measured PSF
xi = 0:dx2:7*dx2; %sampling of PSF in object space (15x15 matrix)
xi = [-xi(end:-1:2), xi];

%Interpolation in object space, scale if dx_psf != dx_imgHigh
[xx, yy] = meshgrid(x, x);
[xxi, yyi] = meshgrid(xi, xi);
PSF = interp2(xx, yy, psf, xxi, yyi); 

% PSF = psf;

%Force border ones (outside measured range if any) to be at background level
nanIdx = isnan(PSF(:));
PSF(nanIdx) = median(psf(:));

% make OTF
[OTF, mask] = makeOTF(PSF, [x_imgHigh, y_imgHigh], 0.2); % high resolution, 0.2-cutoff value in fourier domain

%high resolution
xs = fft2(imgHighRes);
xs = xs.*OTF;
xs = ifft2(xs); 


%% Step4: add noise and save data
if noise_flag % add noise
    img = addGaussianNoise(xs,mean_gaussian, var_gaussian, ...
        strcat(target_folder, name2save_txt)); % txt record starts from two
end
stack_data = cat(3, imgHighRes, img); % add noise-less image at first frame
saveMultipageTiff(stack_data, strcat(target_folder, name2save_tif));


%% Step5: show original and psf-convoluted 
figure(1);
subplot(1,3,1); imagesc(imgHighRes); title('Original image')
subplot(1,3,2); imagesc(xs); title('Psf-convoluted image')
subplot(1,3,3); imagesc(squeeze(stack_data(:,:, randi(size(stack_data,3))))); title('Noise-added image')






