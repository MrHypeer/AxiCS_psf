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
data_folder = '../data_psf/Confocal_Data/append/';
imageName = ''; %other: im2
lowRes_name = 'FOVStack_Dz09Lay36_lowRes.tif';
fullRes_name = 'FOVStack_Dz09Lay36.tif';
xPxl = 512;
im_raw = readMultiTiff([data_folder, imageName, lowRes_name],0); % subsample data
imFull_raw = readMultiTiff([data_folder, imageName, fullRes_name],0); % full resolution data

%% reconstruction volume along x-axis
imFull_rec = zeros(size(imFull_raw));
samples_rec = imFull_rec;
psnr_yz = zeros(1,xPxl);
psnr_xy = zeros(1,size(imFull_raw,3));
tic
for ii_x = 1:xPxl
    im = squeeze(im_raw(:,ii_x,:));
    imFull = squeeze(imFull_raw(:,ii_x,:));
    if (recType == 2)
        psfIm = imread(strcat(data_folder,'5-PsfGen1.tif'));
        cutoff = 0.15;      %Cut-off     
    else
        PSF = [];
        cutoff = 0;
    end
    ratio = [2 2];          %Sub-sampling ratio % change to x-y ratio

    %% Adjusting PSF size to measurement wavelength
    if (recType == 2)
        %Pixel size for PSF measured by beads
        dx = 100; %nm
        %Emission of beads
        lambda = 520; %nm

        %Pixel size of actual images (cells, ER-tracker)
        dx2 = 200; %nm, px size for ER-tracker with zoom 1.2x
        lambda2 = 570; %nm, ER-tracker

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

    
    [rec, samples] = ConfocalCS(im, ratio, recType, PSF, cutoff);
    psnr_yz(ii_x) = psnr(rec,imFull);
    imFull_rec(:,ii_x,:) = rec;
    samples_rec(:,ii_x,:) = samples;
    
end
t = toc;

%% compare PSNR and SNR along x-y plane

for ii_z = 1:size(imFull_raw,3)
   psnr_xy(ii_z) = psnr(imFull_rec(:,:,ii_z),imFull_raw(:,:,ii_z)); 
end
psnr(imFull_rec(:,:,:),imFull_raw(:,:,:));
%% save and show result 
str = split(fullRes_name,'.');
saveMultipageTiff(imFull_rec, [data_folder, str{1},'_PSFCS_psfGen.tif'])


slice_idx = randi([100,400]); % randomly select a slice to see the result
figure(1)
subplot(221), imagesc(squeeze(im_raw(:,slice_idx,:))), title('Measured low-res image')
subplot(222), imagesc(squeeze(samples_rec(:,slice_idx,:))), title('Sampling grid');
subplot(223), imagesc(squeeze(imFull_rec(:,slice_idx,:))), title('Reconstruction')
subplot(224), imagesc(squeeze(imFull_raw(:,slice_idx,:))), title('Measured full-res image (comparison)');
disp(['Reconstruction time: ', num2str(t), ' s']);
sgtitle(['idx:',num2str(slice_idx)])

%% plot intensity profile along sample
% normalize two contrasting objects
line_idx = 20;
rec_normline = mapminmax(rec(line_idx,:),0,1);
imFull_normline = mapminmax(imFull(line_idx,:),0,1);
% figure(2)
% plot(1:size(rec,2), rec_normline, 1:size(imFull,2), circshift(imFull_normline,-4))
% legend('reconstructed','full resolution')



%% %%%%%%%%%%%%%%%%%%%%%%%appendix function
function stacked_img = readMultiTiff(filename, twochannelflag)
    %% this section reads single .tiff file with same name in format (x,y,z)
    tiffname = filename;
    twochannel_flag = twochannelflag;
    info = imfinfo(tiffname);

    if length(info) >= 1
        frame_num = length(info);
        width = info(1).Width;
        height = info(1).Height;
    end

    stacked_img = zeros(height, width, frame_num);

    for ii_frame = 1:frame_num
        stacked_img(:,:,ii_frame) = imread(tiffname,ii_frame);
    end
end

function saveMultipageTiff(data, filename)
% this code compress data to stacked tiff, data should be in (x,y,z)
% format
    t = Tiff(filename, 'w');
    tagstruct.ImageLength = size(data, 1);
    tagstruct.ImageWidth = size(data, 2);
    tagstruct.Compression = Tiff.Compression.None;
    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 64;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    for ii=1:size(data,3)
       setTag(t,tagstruct);
       write(t,data(:,:,ii));
       writeDirectory(t);
    end
%     t.setTag(tagstruct); % used for multichannel writing
%     t.write(data);
    t.close();
end


