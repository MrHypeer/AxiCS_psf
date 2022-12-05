
close all; clear all
%% list all data to be processed
% sort data in natural order
file_fodler = 'D:\OneDrive - City University of Hong Kong\CUHK\08 project\14 CS\data_psf\Confocal_Data\append\beadGen\';
oldFolder = cd (file_fodler);

file_pattern = 'bead_spar2_PSF-5pxl_noisestack*.tif';
d = dir(file_pattern);
nameCell = cell(length(d)-2,1);

for i = 1:length(d)
    temp = d(i).name;
    nameCell{i} = temp;
end
d2 = sort_nat(nameCell);

cd(oldFolder)
% compose new file with natural order
% for ii_file = 1:length(d2)
%    d2(ii_file) = strcat(splitname(1),'_',cell2mat(d2(ii_file)),'.tif');
% end

%% load data and calculate snr value
snr_method = 2; % 1-whole volume, 2-single slice along x-axis
z_fisrt = 0; % 0-full raw image(x,y,z), 1-rec/subsample image(z,y,x)
naturalz = 1; % 0-snr along z-axis, 1-snr along x-axis

% save snr value into a matrix, file_num * rec layer
snr = [];
data_full = readMultiTiff(strcat(file_fodler, cell2mat(d2(2))),0);

for ii_stack = 1:1
    data = readMultiTiff(strcat(file_fodler, cell2mat(d2(ii_stack))),0);
    snr = [snr;calsnr(data, z_fisrt, snr_method,naturalz)];
    % calculate PSNR with original image(no-noise image)
    psnr_val = zeros(1, size(data,3));
    for ii_psnr = 1:size(data,3)
        psnr_val(ii_psnr) = psnr(squeeze(data(:,:,ii_psnr)), squeeze(data_full(:,:,1)));
    end
end





%% plot data
% figure(1)
% plot(snr');
% title('full-res SNR comp. with diff. integration frames')
% legend(patt_extr(1:ii_stack))
% axis([1 size(snr,2) min(min(snr)) max(max(snr))+1])


%% analysis from saved data
% clear all
% load('snr_pollen_yz')
% bar_data = cat(3,snr_yz_raw, snr_yz_recraw, snr_yz_rec);
% % plot bar graph under same avg number, each slice has three data points
% idx_avg = 7; % [2 5 8 11 20 30 40]
% figure (2) 
% plot(squeeze(bar_data(idx_avg,:,:)))
% legend('raw','raw-sub','rec')
% title('SNR comp., sum: 40')
% xlabel('rec. plane index/-'); ylabel('SNR/dB')
% axis([1 size(bar_data,2) min(min(min(bar_data))) max(max(max(bar_data)))+1])

%%
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


function snr = calsnr(data, z_first, snr_method, naturalz)
%% decide data dimesion depending on rec or raw image
if z_first == 1 %if z is first dim, swap axis to match original full image
    data = permute(data,[3,2,1]);
    x_pxl = size(data,1);
    y_pxl = size(data,2);
    z_pxl = size(data,3);
elseif z_first == 0 %raw_image
    x_pxl = size(data,1);
    y_pxl = size(data,2);
    z_pxl = size(data,3);
else
    error('please input right img types');
end

%% plot to test rightness of data dimension
% figure(1)
% subplot(121)
% imagesc(squeeze(data(150,:,:)))
% subplot(122)
% imagesc(squeeze(data2(150,:,:)))

%% calculate snr using threshold
% set global signal and noise level
noise_level = multithresh(data); % set binary value to separate noise and signal
seg_I = imquantize(data,noise_level); % split and segment image
switch snr_method
    case 1 % calculate snr in volume
        signal_idx = seg_I == 2; % signal index, label == 2
        noise_idx = seg_I == 1; % nosie index, label == 1

        % calculate power of noise and signal and calculte SNR
        P_sig = rms(data(signal_idx))^2;
        P_noise = rms(data(noise_idx))^2;
        snr = 10*log10(P_sig/P_noise); 

        % disp snr value
        formatSpec = strcat(filename," | SNR: %4.2f \n");
        fprintf(formatSpec, snr);
    case 2 % calculate snr in every single slice using global threshold
        if naturalz
            rec_pxl = z_pxl;
        else
            rec_pxl = x_pxl;
        end
        
        snr = zeros(1,rec_pxl);
        for ii_x = 1:rec_pxl
            temp = squeeze(data(ii_x,:,:));
%             figure(2)
%             imagesc(squeeze(temp))
            signal_idx = temp > noise_level;
            noise_idx = temp < noise_level;                
            P_sig = rms(temp(signal_idx))^2;
            P_noise = rms(temp(noise_idx))^2;
            snr(ii_x) = 10*log10(P_sig/P_noise);
            if ~sum(sum(signal_idx))
                snr(ii_x) = 0;
            end
%             disp([num2str(ii_x),': snr=',num2str(P_sig),'/',num2str(P_noise),'=',num2str(snr(ii_x))])
        end
    otherwise 
        error('please select right calculation domain')
end
end


function [cs,index] = sort_nat(c,mode)
%sort_nat: Natural order sort of cell array of strings.
% usage:  [S,INDEX] = sort_nat(C)
%
% where,
%    C is a cell array (vector) of strings to be sorted.
%    S is C, sorted in natural order.
%    INDEX is the sort order such that S = C(INDEX);
%
% Natural order sorting sorts strings containing digits in a way such that
% the numerical value of the digits is taken into account.  It is
% especially useful for sorting file names containing index numbers with
% different numbers of digits.  Often, people will use leading zeros to get
% the right sort order, but with this function you don't have to do that.
% For example, if C = {'file1.txt','file2.txt','file10.txt'}, a normal sort
% will give you
%
%       {'file1.txt'  'file10.txt'  'file2.txt'}
%
% whereas, sort_nat will give you
%
%       {'file1.txt'  'file2.txt'  'file10.txt'}
%
% See also: sort

% Version: 1.4, 22 January 2011
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Real_email = regexprep(Email,{'=','*'},{'@','.'})


% Set default value for mode if necessary.
if nargin < 2
    mode = 'ascend';
end

% Make sure mode is either 'ascend' or 'descend'.
modes = strcmpi(mode,{'ascend','descend'});
is_descend = modes(2);
if ~any(modes)
    error('sort_nat:sortDirection',...
        'sorting direction must be ''ascend'' or ''descend''.')
end

% Replace runs of digits with '0'.
c2 = regexprep(c,'\d+','0');

% Compute char version of c2 and locations of zeros.
s1 = char(c2);
z = s1 == '0';

% Extract the runs of digits and their start and end indices.
[digruns,first,last] = regexp(c,'\d+','match','start','end');

% Create matrix of numerical values of runs of digits and a matrix of the
% number of digits in each run.
num_str = length(c);
max_len = size(s1,2);
num_val = NaN(num_str,max_len);
num_dig = NaN(num_str,max_len);
for i = 1:num_str
    num_val(i,z(i,:)) = sscanf(sprintf('%s ',digruns{i}{:}),'%f');
    num_dig(i,z(i,:)) = last{i} - first{i} + 1;
end

% Find columns that have at least one non-NaN.  Make sure activecols is a
% 1-by-n vector even if n = 0.
activecols = reshape(find(~all(isnan(num_val))),1,[]);
n = length(activecols);

% Compute which columns in the composite matrix get the numbers.
numcols = activecols + (1:2:2*n);

% Compute which columns in the composite matrix get the number of digits.
ndigcols = numcols + 1;

% Compute which columns in the composite matrix get chars.
charcols = true(1,max_len + 2*n);
charcols(numcols) = false;
charcols(ndigcols) = false;

% Create and fill composite matrix, comp.
comp = zeros(num_str,max_len + 2*n);
comp(:,charcols) = double(s1);
comp(:,numcols) = num_val(:,activecols);
comp(:,ndigcols) = num_dig(:,activecols);

% Sort rows of composite matrix and use index to sort c in ascending or
% descending order, depending on mode.
[unused,index] = sortrows(comp);
if is_descend
    index = index(end:-1:1);
end
index = reshape(index,size(c));
cs = c(index);
end