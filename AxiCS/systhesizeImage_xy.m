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
bmpName = 'USAF-50.bmp';
imgHighRes = imread(bmpName);
psf = imread('psf-limo-2.tif');
suffixLow = '_limo2_sub';
suffixHigh = '_limo2.tif';
noise_flag = 1; % noise switch
mean_gaussian = 0;
var_gaussian = 0.03;
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

PSF = double(psf);

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

[x_pxl, y_pxl] = size(xs);
for ii_subsam = 1:1
   if mod(floor(x_pxl/(ii_subsam+1)),2) == 0
       height = floor(x_pxl/(ii_subsam+1));
       dx = x_pxl/height;
       if height*(ii_subsam+1) < x_pxl
           height = height + 2;
       end
   else
       height = ceil(x_pxl/(ii_subsam+1));
       dx = x_pxl/height;
       if height*(ii_subsam+1) < x_pxl
           height = height + 2;
       end
   end
   if mod(floor(y_pxl/(ii_subsam+1)),2) == 0
       width = floor(y_pxl/(ii_subsam+1));
       dx = y_pxl/width;
       if width*(ii_subsam+1) < y_pxl
           width = width + 2;
       end
   else
       width = ceil(y_pxl/(ii_subsam+1));
       dx = y_pxl/width;
       if width*(ii_subsam+1) < y_pxl
           width = width + 2;
       end
   end

%     x = -(l(1)-1)*dx:dx:l(1)*dx; %x sampling of measured PSF
%     xi = 0:dx3:7*dx3; %sampling of PSF in object space (15x15 matrix)
%     xi = [-xi(end:-1:2), xi];
%     x_low = -(width/2-1)*dx:dx:width/2*dx;
    y_low = -(y_pxl/2-1):y_pxl/width:y_pxl/2;
    y_high = -(y_pxl/2-1):y_pxl/2;
    x_low = -(x_pxl/2-1):x_pxl/height:x_pxl/2;
    x_high = -(x_pxl/2-1):x_pxl/2;
    [xx_low, yy_low] = meshgrid(y_low, x_low);
    [xx_high, yy_high] = meshgrid(y_high, x_high);
   
   img_LowRes = interp2(xx_high, yy_high, xs, xx_low, yy_low, 'cubic',0);
   
   filename = strcat(name(1),suffixLow,num2str(ii_subsam+1),'.tif');

   % add noise
   if noise_flag == 1
       img_LowRes = imnoise(rescale(img_LowRes), 'gaussian', mean_gaussian,...
           var_gaussian);
       img_LowRes = rescale(img_LowRes,0,255);
       filename = strcat(name(1),suffixLow,num2str(ii_subsam+1),'noise',...
           num2str(var_gaussian),'.tif');
   end
   
   temp = zeros(height, width, 2);
   temp(:,:,1) = img_LowRes;
   figure(2)
   imagesc(img_LowRes); title('low-resolution after trim')
   
   saveMultipageTiff(temp(:,:,1),cell2mat(filename))
end

xs2 = zeros(x_pxl,y_pxl,2);
xs2(:,:,1) = xs;
saveMultipageTiff(xs2(:,:,1),strcat(cell2mat(name(1)),suffixHigh))
% saveMultipageTiff(xs4(:,:,1),[cell2mat(name(1)),'-1024.tif'])









