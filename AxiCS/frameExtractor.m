%%%%%%%%%%%%%%%%%%%%%%%
%%% this code extracts different frames from different construction stack
%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all
%% list all data to be processed
% choose data name and slice to grab
file_origin = 'PollenStack.tif';
file_pattern = 'PollenStack_avg*_sub2_CR2.tif'; % load differnet rec images
name2save_xy = 'PollenStack_goodrec_xy38.tif';
name2save_yz = 'PollenStack_goodrec_yz281.tif';
idx_xy = 38; % set slice to generate stack
idx_yz = 281;

% sort file in natural order
d = dir(file_pattern);
nameCell = cell(length(d)-2,1);

for i = 1:length(d)
%     disp(d(i).name);
    temp = d(i).name;
    splitname = split(temp,'_');
    nameCell{i} = char(splitname(2));
end
d2 = sort_nat(nameCell);
patt_extr = d2;

% compose new file with natural order
for ii_file = 1:length(d2)
   d2(ii_file) = strcat(splitname(1),'_',cell2mat(d2(ii_file)),'_',...
       splitname(3),'_',splitname(4));
end

%% load, re-arrange and save data 
% load data
stack_data_xy = [];
stack_data_yz = [];
for ii_stack = 1:length(d2)
    temp = readMultiTiff(cell2mat(d2(ii_stack)),0);
    stack_data_xy = cat(3, stack_data_xy, normMatrix(squeeze(temp(idx_xy,:,:))));
    stack_data_yz = cat(3, stack_data_yz, normMatrix(squeeze(temp(:,:,idx_yz))));
end

data_ori = readMultiTiff(file_origin,0);
stack_data_xy = cat(3, stack_data_xy, normMatrix(transpose(squeeze(data_ori(:,:,idx_xy)))));
stack_data_yz = cat(3, stack_data_yz, normMatrix(transpose(squeeze(data_ori(idx_yz,:,:)))));

% save data
saveMultipageTiff(stack_data_xy, name2save_xy)
saveMultipageTiff(stack_data_yz, name2save_yz)

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


function norm_value = normMatrix(input)
%%% this function morms input into range 0-255(uint8), canbe 3d and 2d input
%%% default input is (x,y,z)
norm_value = uint8(zeros(size(input)));
if length(size(input)) == 3
    % process 3d data
    for ii_stack = 1:size(input,3)
        max_int = max(max(input(:,:,ii_stack)));
        min_int = min(min(input(:,:,ii_stack)));
        norm_value(:,:,ii_stack) = (input(:,:,ii_stack)-min_int)./(max_int ...
            -min_int)*255;
    end
elseif length(size(size(input) == 2))
    % process 2d data
    max_int = max(max(input));
    min_int = min(min(input));
    norm_value = (input-min_int)./(max_int ...
        -min_int)*255;
else 
    % error
    error('Input dim shoule be 2D or 3D')
end
end