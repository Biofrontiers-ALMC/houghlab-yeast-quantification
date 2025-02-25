%Process images to measure change in pH

clearvars
clc

inputDir = 'Z:\Microscopy\Yeast\Sup35\20250205_JA_SEP-mRuby';
outputDir = 'Z:\Microscopy\Yeast\Sup35\20250224 pH JWT';
files = dir(fullfile(inputDir, '*.nd2'));

storeTime = zeros(1, numel(files));
storeRatio = zeros(1, numel(files));

for iF = 1:numel(files)

    reader = BioformatsImage(fullfile(files(iF).folder, files(iF).name));

    Iegfp = getPlane(reader, 1, 1, 1);
    Itritc = getPlane(reader, 1, 2, 1);

    %Mask from Itritc
    cellMask = imbinarize(Itritc, 'adaptive');

    cellMask = imopen(cellMask, strel('disk', 2));
    cellMask = imfill(cellMask, 'holes');

    cellMask = imclearborder(cellMask);

    cellMask = bwareaopen(cellMask, 100);

    dd = -bwdist(~cellMask);
    dd(~cellMask) = Inf;
    dd = imhmin(dd, 0.5);

    LL = watershed(dd);

    cellMask(LL == 0) = 0;

    cellMask = bwareaopen(cellMask, 400);

    Iout = showoverlay(Itritc, bwperim(cellMask));


    [~, fn] = fileparts(reader.filename);
    imwrite(Iout, fullfile(outputDir, [fn, '_mask.tif']), 'Compression', 'none');

    cellData_red = regionprops(cellMask, Itritc, 'MeanIntensity');
    cellData_green = regionprops(cellMask, Iegfp, 'MeanIntensity');

    for iCell = 1:numel(cellData_red)

        cellData(iCell).TRITCint = cellData_red(iCell).MeanIntensity;
        cellData(iCell).EGFPint = cellData_green(iCell).MeanIntensity;

        cellData(iCell).ratio = cellData_green(iCell).MeanIntensity./cellData_red(iCell).MeanIntensity;

    end

    storeRatio(iF) = mean([cellData.ratio]);
    
    %Convert to POSIX time, which is number of seconds since Jan 1, 1970
    storeTime(iF) = convertTo(datetime(files(iF).date), 'posixTime');

end

%%
minTime = min(storeTime);
storeTime = storeTime - minTime;

[tt, idx] = sort(storeTime);

storeRatio_sorted = storeRatio(idx);

plot(tt/60, storeRatio_sorted)

ylabel('Mean ratio')
xlabel('Minutes')


