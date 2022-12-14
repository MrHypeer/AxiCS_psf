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
ratio = [2 2];          %Sub-sampling ratio % change to x-y ratio

%% Retrieve data
data_folder = '..\data_psf\Confocal_Data\append\USAF\';
stackName = 'USAF1951G7PSF-5pxl-64_noisestack2.tif'; %other: im2
name2save = 'USAF1951G7PSF-5pxl-64_noisestack2_CR2_recparam.tif';
dat2save = 'USAF1951G7PSF-5pxl-64_noisestack2_recparam.mat';
log2save = 'USAF1951G7PSF-5pxl-64_noisestack2_recparam.txt';

im_Full = readMultiTiff([data_folder, stackName],0); % full resolution data
im_ori = squeeze(im_Full(:,:,1)); % original data without noise
im_raw = im_Full(1:ratio(1):end,1:ratio(2):end,:); % subsampling data using CR=2
rec_z = 1; % if 1, rec along z layers, if 0, rec along x-axis
opt_tuning = 1; % if 1, opts will adaptively change whether using grid search or automatical calculation

%% reconstruction volume along x-axis
x_full = size(im_Full,1);
y_full = size(im_Full,2);
z_full = size(im_Full,3);

recPxl = x_full;
if rec_z
    recPxl = z_full;
end

imFull_rec = []; % used for normalization
samples_rec = [];
psnr_rec = []; % compare with noise-less img
psnr_ortrec = [];
count = 0;

tic
fileID = fopen(strcat(data_folder, log2save),'w'); % save reconstruciton log
for ii_rec = 1:recPxl
    disp(['current layer:',num2str(ii_rec),' /',num2str(recPxl)])
    im = squeeze(im_raw(:,:,ii_rec));
%     imFull = squeeze(imFull_raw(:,ii_x,:));
    if (recType == 2)
        psfIm = imread(strcat(data_folder,'PSF-5pxl-64.tif')); % 
        cutoff = 0.15;      %Cut-off     
    else
        PSF = [];
        cutoff = 0;
    end

    %% Adjusting PSF size to measurement wavelength
    if (recType == 2)
        %Pixel size for PSF measured by beads
        dx = 100; %nm
        %Emission of beads
        lambda = 567; %nm

        %Pixel size of actual images (cells, ER-tracker)
        dx2 = 730; %nm, px size for ER-tracker with zoom 1.2x
        lambda2 = 511; %nm, ER-tracker

        %Scaling the pixel size with the wavelength
        % why pixel size is related with wavelength???
        dx3 = dx2*lambda/lambda2;

        %Create the xy map (assuming even-sized square image)
        l = floor(size(psfIm)./2);
        x = -(l(1)-1)*dx:dx:l(1)*dx; %x sampling of measured PSF
        xi = 0:dx3:21*dx3; %sampling of PSF in object space (15x15 matrix)
        xi = [-xi(end:-1:2), xi];

        %Interpolation in object space
        [xx, yy] = meshgrid(x, x);
        [xxi, yyi] = meshgrid(xi, xi);
        PSF = interp2(xx, yy, double(psfIm), xxi, yyi);
        %Force border ones (outside measured range if any) to be at background level
        nanIdx = find(isnan(PSF(:)));
        PSF(nanIdx) = median(psfIm(:));    
    end

    % set reconstruction params
    clear opts
    if opt_tuning % change three: mu, beta, tol
        for mu = 4:13
            for beta =  4:13
                for tol = [1e-1 1e-2 1e-3 1e-4 1e-5 1e-6]
                    count = count + 1;
                    opts.mu = 2^mu;
                    opts.beta = 2^beta;
                    opts.tol = tol;
                    opts.maxit = 300;
                    opts.TVnorm = 1;
                    opts.nonneg = true;
                    opts.isreal = true;
                    
                    % execute reconstruction
                    [rec, samples] = AxiCS(im, ratio, recType, PSF, cutoff,opts);
                    temp_psnr = psnr(rec./max(max(rec)),im_ori./max(max(im_ori)));
                    psnr_rec = cat(1, psnr_rec, temp_psnr);
                    imFull_rec = cat(3,imFull_rec, rec);
                    samples_rec = cat(3,samples_rec, samples);
                    
                    % write reconstruction details
                    formatSpec = '%2.0u \t %2.0u th rec layer,\t mu:%4.2f,\t beta:%4.2f,\t tol:%4.2e,\t psnr:%4.2f\n';
                    fprintf(fileID, formatSpec,...
                            [count ii_rec, opts.mu, opts.beta, opts.tol, temp_psnr]);
                end
            end
        end
    else % default params sutiable for tpm
        count = count + 1;
        opts.mu = 2^4;
        opts.beta = 2^4;
        %To push quality in case of low-noise signal
        % opts.mu = 2^12; opts.beta = 2^9;
        opts.tol = 1E-3;
        opts.maxit = 300;
        opts.TVnorm = 1;
        opts.nonneg = true;
        opts.isreal = true;
        
        % execute reconstruction
        [rec, samples] = AxiCS(im, ratio, recType, PSF, cutoff,opts);
        temp_psnr = psnr(rec./max(max(rec)),im_ori./max(max(im_ori)));
        psnr_rec = cat(1, psnr_rec, temp_psnr);
        imFull_rec = cat(3,imFull_rec, rec);
        samples_rec = cat(3,samples_rec, samples);
        
        % write reconstruction details
        formatSpec = '%2.0u \t %2.0u th rec layer,\t mu:%4.2f,\t beta:%4.2f,\t tol:%4.2e,\t psnr:%4.2f\n';
        fprintf(fileID, formatSpec,...
                            [count ii_rec, opts.mu, opts.beta, opts.tol, temp_psnr]);
    end
end
fclose(fileID); % close file
t = toc;

%% compare PSNR and SNR along orthonogal reconstruction direction
% for ii_ort = 1:x_full
%    psnr_ortrec(ii_ort) = psnr(imFull_rec(:,:,ii_ort),imFull_raw(:,:,ii_ort)); 
% end
% psnr(imFull_rec(:,:,:),imFull_raw(:,:,:));
%% save result and data
saveMultipageTiff(imFull_rec, strcat(data_folder, name2save))
save(strcat(data_folder, dat2save), 'psnr_rec')

%% show result
% slice_idx = randi(x_full); % randomly select a slice to see the result
% figure(1)
% subplot(221), imagesc(squeeze(im_raw(:,:,slice_idx))), title('Measured low-res image')
% subplot(222), imagesc(squeeze(samples_rec(:,:,slice_idx))), title('Sampling grid');
% subplot(223), imagesc(squeeze(imFull_rec(:,:,slice_idx))), title('Reconstruction')
% temp = squeeze(imFull_raw(slice_idx,:,:));
% subplot(224), imagesc(temp'), title('Measured full-res image (comparison)');
% disp(['Reconstruction time: ', num2str(t), ' s']);
% sgtitle(['idx:',num2str(slice_idx)])

%% plot intensity profile along sample
% normalize two contrasting objects
% line_idx = 20;
% rec_normline = mapminmax(rec(line_idx,:),0,1);
% imFull_normline = mapminmax(imFull(line_idx,:),0,1);
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


