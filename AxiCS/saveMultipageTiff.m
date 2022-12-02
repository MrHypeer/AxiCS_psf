function saveMultipageTiff(data, filename)
% this code compress data to stacked tiff, data should be in (x,y,z)
% format
    if length(size(data)) == 2
        temp(:,:,2) = data;
        data = temp(:,:,2);
    end
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