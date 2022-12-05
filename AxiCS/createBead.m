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
gaussianNoise_var = [0.01 0.03 0.05 0.07]; % gaussian noise
gaussianNoise_mean = [0 0.1 0.2];
poissionNoise = 0; % add or not possiaon noise
psf_name = 'PSF-5pxl.tif';
target_folder = 'data_psf\Confocal_Data\append\beadGen\';

%% define name of saved file
% working folder changing
current_path = pwd;
pos_v = strfind(current_path,'\');
p = current_path(1:pos_v(length(pos_v)-1)); % -1: delete the last character '/' or '\' 

oldFolder = cd([p, target_folder]);
psfname_split = split(psf_name,'.');
order = 1;
name2save = strcat('bead_spar',num2str(sparsity),'_', ...
    cell2mat(psfname_split(1)), '_noisestack',num2str(order));
name2save_tif = strcat(name2save,'.tif');
name2save_txt = strcat(name2save,'.txt');

% file alreaady exist, automatically +1
while ~isempty(dir(name2save_txt)) % means already file has same name
    order = 1 + order;
    name2save = strcat('bead_spar',num2str(sparsity),'_', ...
        cell2mat(psfname_split(1)), '_noisestack',num2str(order));
    name2save_tif = strcat(name2save,'.tif');
    name2save_txt = strcat(name2save,'.txt');
end
cd(oldFolder)
%% bead generation
image = makeBead(sparsity, show_flag);
% soft image using psf
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
stack_data = zeros(512,512, length(gaussianNoise_mean)*length(gaussianNoise_var));
count = 1; % start from one, corresponding with tiff stack
if sum(gaussianNoise_var) ~= 0
    for ii_mean = 1: length(gaussianNoise_mean)
        for ii_var = 1:length(gaussianNoise_var)
%             disp(['add Gaussian noise, mean/var: ', num2str(gaussianNoise_mean),...
%                 '/',num2str(gaussianNoise_var)])
            count = count + 1;
            temp = imnoise(imageSoft, 'gaussian', gaussianNoise_mean(ii_mean), ...
                gaussianNoise_var(ii_var));
            %% add possion noise
            if poissionNoise == 1
                disp('add Poisson noise')
                temp = imnoise(temp, 'poisson');
            end
            stack_data(:,:,length(gaussianNoise_var)*(ii_mean-1)+ii_var) = ...
                temp;
            fileID = fopen(strcat(p, target_folder, name2save_txt),'a');

            formatSpec = '%2.0u: gaussian nosie mean/var is %4.2f/%4.2f, poisson noise is %4.2f\n';
            fprintf(fileID, formatSpec,...
                [count, gaussianNoise_mean(ii_mean), gaussianNoise_var(ii_var), poissionNoise]);
            fclose(fileID);
        end
    end
end
stack_data = cat(3, imageClean, stack_data);
% % Leaving 10% of amplitude for background
% imageSoft = imageSoft*0.9*signalStrength + 0.1*signalStrength;

%% plot
figure(1);
imagesc(imageSoft); title(['generated bead with Gaussian noise: ',...
    num2str(gaussianNoise_mean),'/',num2str(gaussianNoise_var)])
xlabel('x/pxl');ylabel('y/pxl')


%% save data into tiff
saveMultipageTiff(stack_data, strcat(p, target_folder,name2save_tif))
