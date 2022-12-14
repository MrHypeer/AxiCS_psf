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
imageName = 'bead'; %other: im2

psf_pxl = [3,4,5,6];
sub_samp = [2,3,4,5,6,7,8,9,10];
metrics = zeros(length(sub_samp),3,length(psf_pxl));

data_folder = '../data_psf/Confocal_Data/append/simulationGeneration/';
for ii_psf = 1:length(psf_pxl)
    stack_data = zeros(1050,1050,length(sub_samp)+1);
    name2save = strcat(data_folder, imageName, '_', num2str(psf_pxl(ii_psf)),...
        'pxl_stack.tif');
    for ii_subsamp = 1:length(sub_samp)
        im = imread(strcat(data_folder, imageName, '_', num2str(psf_pxl(ii_psf)),...
        'pxl_sub',num2str(sub_samp(ii_subsamp)),'.tif'));
        imFull = imread(strcat(data_folder, imageName, '_',num2str(psf_pxl(ii_psf)),'pxl.tif'));

        ratio = [sub_samp(ii_subsamp) sub_samp(ii_subsamp)];          %Sub-sampling ratio % change to x-y ratio

        if (recType == 2)
            psfIm = double(imread(strcat(data_folder,'PSF-',num2str(psf_pxl(ii_psf)),'pxl.tif')));
            cutoff = 0.15;      %Cut-off     
        else
            PSF = [];
            cutoff = 0;
        end

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
        %     PSF = interp2(xx, yy, psfIm, xxi, yyi);
            PSF = psfIm; % not sampling when dx_psf == dx_imgHigh

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
        % resize reconstructed iamge into original size
        [x_rec, y_rec] = size(rec);
        [x_ori, y_ori] = size(imFull);
        xx_rec = -(x_ori/2-1):x_ori/(x_rec+1):x_ori/2;
        xx_ori = -(x_ori/2-1):x_ori/2;
        [XX_rec,YY_rec] = meshgrid(xx_rec, xx_rec);
        [XX_ori, YY_ori] = meshgrid(xx_ori, xx_ori);
        imFull_trim = interp2(XX_ori, YY_ori, imFull, XX_rec, YY_rec);
        % normalize two contrasting objects, shown in scaled size
        idx = randi(size(rec,1));
        rec_norm = mapminmax(rec,0,1);
        imFull_norm = mapminmax(imFull_trim,0,1);
        % plot intensity profile along sample
        figure(2)
        plot(1:size(rec,2), rec_norm(idx,:), 1:size(imFull_trim,2), circshift(imFull_norm(idx,:),0))
        legend('reconstructed','full resolution')
        title(strcat('line across x=',num2str(idx)))
        % calculate PSNR
        [psnr_value, snr_value] = psnr(rec_norm, imFull_norm);
        RMSE = sqrt(immse(rec_norm, imFull_norm));
        disp(strcat('psf/sub:',num2str(psf_pxl(ii_psf)),'/',num2str(sub_samp(ii_subsamp))))
        disp(strcat('psnr,snr,rmse: ',num2str(psnr_value),',',num2str(snr_value),',',num2str(RMSE)))
        metrics(ii_subsamp,:,ii_psf) = [psnr_value, snr_value, RMSE];
        %% save image

        stack_data(1:size(rec,1),1:size(rec,2),ii_subsamp+1) = imFull_trim;
    end
    stack_data(1:size(rec,1),1:size(rec,2),1) = imFull_trim;
    saveMultipageTiff(stack_data, name2save)
end

%% %%%%%%%%%%%%%%%
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





