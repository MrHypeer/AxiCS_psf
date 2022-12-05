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