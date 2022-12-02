%%%%%%%%
% this code generates bead sample in 512*512 image
% density can be controlled using
%%%%%%%%
clear all

%% set params
sparsity = 2; 
show_flag = 0;
add_psf = 1;
signalStrength = 255;
gaussianNoise_var = 0.07; % gaussian noise
gaussianNoise_mean = 0.1;
poissionNoise = 0; % add or not possiaon noise
psf_name = 'PSF-5pxl.tif';
target_folder = 'data_psf\Confocal_Data\append\beadGen\';


%% working folder changing
current_path = pwd;
pos_v = strfind(current_path,'\');
p = current_path(1:pos_v(length(pos_v)-1)); % -1: delete the last character '/' or '\' 
%% bead generation
image = makeBead(sparsity, show_flag);
%% soft image using psf
if add_psf
	PSF = imread(psf_name);
    %Force border ones (outside measured range if any) to be at background level
    nanIdx = isnan(PSF(:));
    PSF(nanIdx) = median(PSF(:));
    % make OTF
    [OTF, mask] = makeOTF(PSF, [size(image,1), size(image,2)], 0.2); % high resolution, 0.2-cutoff value in fourier domain
    %convolute psf and image
    imageSoft = fft2(image);
    imageSoft = imageSoft.*OTF;
    imageSoft = ifft2(imageSoft);  
end

%% normalization
imageSoft = imageSoft-min(imageSoft(:));
imageSoft = imageSoft/max(imageSoft(:))/1.2; % divide 1.1 for noise addition
imageClean = imageSoft;
%% add gaussian noise
if gaussianNoise_var ~= 0
    disp(['add Gaussian noise, var: ', num2str(gaussianNoise_var)])
    imageSoft = imnoise(imageSoft, 'gaussian', gaussianNoise_mean, gaussianNoise_var);
end
% % Leaving 10% of amplitude for background
% imageSoft = imageSoft*0.9*signalStrength + 0.1*signalStrength;

%% add possion noise
disp('add Poisson noise')
if poissionNoise == 1
    imageSoft = imnoise(imageSoft, 'poisson');
end

%% plot
figure(1);
imagesc(imageSoft); title(['generated bead with Gaussian noise: ',...
    num2str(gaussianNoise_mean),'/',num2str(gaussianNoise_var)])
xlabel('x/pxl');ylabel('y/pxl')
% figure(2)
% plot(imageSoft(243,:))

%% save data into tiff
psfname_split = split(psf_name,'.');
name2save_noise = strcat('bead_spar',num2str(sparsity), '_', ...
    cell2mat(psfname_split(1)), '_gauM', num2str(gaussianNoise_mean)...
    ,'V',num2str(gaussianNoise_var),'poi',num2str(poissionNoise),'.tif');
name2save_clean = strcat('bead_spar',num2str(sparsity),'_', ...
    cell2mat(psfname_split(1)), '.tif');

saveMultipageTiff(imageClean, strcat(p, target_folder,name2save_clean))
saveMultipageTiff(imageSoft, strcat(p, target_folder, name2save_noise))
