%%%
% generate high_res and low_res using psf connnnnnvolution

clear all;
%% %%%%%%%%%%
% Step1: load original image, psf and convolution

%%% generate bead sample
% img_bead = zeros(1024,1024);
% for ii_bead = 1:100
%     idx = randi(1024*1024);
%     img_bead(idx) = 1;
% end
% imwrite(img_bead,'bead.tif')
%%%
bmpName = 'sharpAndSmooth.tif';
imgHighRes = imread(bmpName);
psf = imread('PSF-3pxl.tif');
suffix = '_3pxl_sub';
suffix2 = '_3pxl.tif';
noise_flag = 1; % noise switch
mean_gaussian = 0;
var_gaussian = 0.07;
% wavelength information
lambda_psf = 511;
lambda_img = 511;

% pixel size information
[x_psf, y_psf] = size(psf);
dx_psf = 100; % nm, same as dz_psf
[x_imgHigh, y_imgHigh] = size(imgHighRes);
dx_imgHigh = 100; %unit, nm



%% %%%%%%%%%
% Step2: scale psf into image space and make OTF

%Create the xy map (assuming even-sized square image), low res
dx2 = dx_imgHigh*lambda_psf/lambda_img;
l = floor([x_psf, y_psf]./2);
x = -(l(1)-1)*dx_psf:dx_psf:l(1)*dx_psf; %x sampling of measured PSF
xi = 0:dx2:7*dx2; %sampling of PSF in object space (15x15 matrix)
xi = [-xi(end:-1:2), xi];

% %Interpolation in object space, scale if dx_psf != dx_imgHigh
% [xx, yy] = meshgrid(x, x);
% [xxi, yyi] = meshgrid(xi, xi);
% PSF = interp2(xx, yy, psf, xxi, yyi); 

PSF = psf;

%Force border ones (outside measured range if any) to be at background level
nanIdx = isnan(PSF(:));
PSF(nanIdx) = median(psf(:));

% make OTF
[OTF, mask] = makeOTF(PSF, [x_imgHigh, y_imgHigh], 0.2); % high resolution, 0.2-cutoff value in fourier domain

%% %%%%%%
% Step3: convolute psf and image
%high resolution
xs = fft2(imgHighRes);
xs = xs.*OTF;
xs = ifft2(xs); 


%% %%%%
% Step4: show original and psf-convoluted 
figure(1);
subplot(1,2,1); imagesc(imgHighRes); title('Original image')
subplot(1,2,2); imagesc(xs); title('Psf-convoluted image')

%% %%%%
% Step5: save data   
name = split(bmpName,'.'); % pattern name

% generate subsampling using interpolation, low-res iamge should have even
% dim
for ii_subsam = 1:9
   if mod(floor(1024/(ii_subsam+1)),2) == 0
       width = floor(1024/(ii_subsam+1));
       dx = 1024/width;
       if width*(ii_subsam+1) < 1024
           width = width + 2;
       end
   else
       width = ceil(1024/(ii_subsam+1));
       dx = 1024/width;
       if width*(ii_subsam+1) < 1024
           width = width + 2;
       end
   end
%     x = -(l(1)-1)*dx:dx:l(1)*dx; %x sampling of measured PSF
%     xi = 0:dx3:7*dx3; %sampling of PSF in object space (15x15 matrix)
%     xi = [-xi(end:-1:2), xi];
%     x_low = -(width/2-1)*dx:dx:width/2*dx;
    x_low = -(1024/2-1):1024/width:1024/2;
    x_high = -(1024/2-1):1024/2;
    [xx_low, yy_low] = meshgrid(x_low, x_low);
    [xx_high, yy_high] = meshgrid(x_high,x_high);
   
   img_LowRes = interp2(xx_high, yy_high, xs, xx_low, yy_low, 'cubic',0);
   
   filename = strcat(name(1),suffix,num2str(ii_subsam+1),'.tif');

   % add noise
   if noise_flag == 1
       img_LowRes = imnoise(rescale(img_LowRes), 'gaussian', mean_gaussian,...
           var_gaussian);
       filename = strcat(name(1),suffix,num2str(ii_subsam+1),'noise',...
           num2str(var_gaussian),'.tif');
   end
   
   temp = zeros(width, width, 2);
   temp(:,:,1) = img_LowRes;
   figure(2)
   imagesc(img_LowRes); title('low-resolution after trim')
   
   saveMultipageTiff(temp(:,:,1),cell2mat(filename))
end

xs2 = zeros(1024,1024,2);
xs2(:,:,1) = xs;
saveMultipageTiff(xs2(:,:,1),strcat(cell2mat(name(1)),suffix2))
% saveMultipageTiff(xs4(:,:,1),[cell2mat(name(1)),'-1024.tif'])









