clearvars
clc

inputDir = 'Z:\Microscopy\Yeast\Sup35\20250211_JA_Sup35MstartWT_4hrinduction';
outputDir = 'Z:\Microscopy\Yeast\Sup35\20250214 Analysis JWT\MATLAB';

% inputDir = 'D:\Projects\ALMC Tickets\Hough\data';
% outputDir = 'D:\Projects\ALMC Tickets\Hough\processed\20250214_c';

%Size of spot filter
sigma1 = 1/(1 + sqrt(2)) * 6;
sigma2 = 1/(1 + sqrt(2)) * 12;

if ~exist(outputDir, 'dir')
    mkdir(outputDir)
end

%% Process all files in directory

files = dir(fullfile(inputDir, '*.nd2'));

for iFile = 1:numel(files)

    fprintf('[%s] Processing %s (file %d of %d)...\n', ...
        datetime, files(iFile).name, iFile, numel(files))

    reader = BioformatsImage(fullfile(files(iFile).folder, files(iFile).name));

    [~, outputFN] = fileparts(reader.filename);

    mask = zeros(reader.height, reader.width, reader.sizeZ, 'logical');
    spotMask = zeros(reader.height, reader.width, reader.sizeZ, 'logical');
    image = zeros(reader.height, reader.width, reader.sizeZ, 'uint16');

    for iZ = 1:reader.sizeZ

        I = getPlane(reader, iZ, 1, 1);

        image(:, :, iZ) = I;

        %Find cells
        currMask = imbinarize(I, 'adaptive');
        currMask = imopen(currMask, strel('disk', 3));
        currMask = imfill(currMask, 'holes');

        mask(:, :, iZ) = currMask;

        %Find spots
        df1 = imgaussfilt(image(:, :, iZ), sigma1);
        df2 = imgaussfilt(image(:, :, iZ), sigma2);

        DoG = df1 - df2;

        currSpotMask = DoG > 200;

        tmp_spotData = regionprops(currSpotMask, 'Circularity', 'PixelIdxList');

        for iSpot = 1:numel(tmp_spotData)

            if tmp_spotData(iSpot).Circularity < 0.8

                currSpotMask(tmp_spotData(iSpot).PixelIdxList) = false;

            end

        end

        spotMask(:, :, iZ) = currSpotMask;


        %%imshow(mask)
    end

    %%
    mask3 = imopen(mask, strel('sphere', 4));

    %%
    dd = -bwdist(~mask3);
    dd(~mask3) = Inf;

    dd = imhmin(dd, 0.8);

    LL = watershed(dd);
    LL(~mask3) = 0;

    % volshow(LL)
    % %% Spot finding
    % sigma1 = 6;
    % sigma2 = 3;
    % 
    % df1 = imgaussfilt3(image, sigma1);
    % df2 = imgaussfilt3(image, sigma2);
    % 
    % DoG = df2 - df1;
    % 
    % spotMask = DoG > 250;


    %%
    % largeObjs = bwareaopen(spotMask, 500, 26);
    % 
    % spotMask(largeObjs) = false;

    %volshow(spotMask)

    %%
    mask3(LL == 0) = false;
    %volshow(mask3)


    %% Visualize by generating zStack of images
    for iZ = 1:size(image, 3)

        imgOut = showoverlay(image(:, :, iZ), bwperim(mask3(:, :, iZ)), 'Color', [0 1 0]);
        imgOut = showoverlay(imgOut, bwperim(spotMask(:, :, iZ)), 'Color', [1 0 1]);

        if iZ == 1
            imwrite(imgOut, fullfile(outputDir, [outputFN, '.tiff']), 'Compression', 'none')
        else
            imwrite(imgOut, fullfile(outputDir, [outputFN, '.tiff']), 'Compression', 'none', 'WriteMode', 'append')
        end

    end


    celldata = regionprops(mask3, image, 'PixelIdxList', 'PixelList', 'PixelValues');
    spotData = regionprops(spotMask, image, 'Centroid', 'MeanIntensity', 'PixelValues', 'PixelIdxList');

    spotCentroids = cat(1, spotData.Centroid);

    %%
    for iCell = 1:numel(celldata)

        celldata(iCell).NumSpots = 0;
        celldata(iCell).TotalSpotInt = 0;
        celldata(iCell).SpotPixelIdxList = {};
        celldata(iCell).DiffusePixelIdxList = celldata(iCell).PixelIdxList;

        for iSpot = 1:numel(spotData)
            %Check if spot is in cell
            if ismember(round(spotCentroids(iSpot, :)), celldata(iCell).PixelList, 'rows')

                %Store spot information as well

                celldata(iCell).NumSpots = celldata(iCell).NumSpots + 1;
                celldata(iCell).TotalSpotInt = celldata(iCell).TotalSpotInt + ...
                    sum(spotData(iSpot).PixelValues);

                if isempty(celldata(iCell).SpotPixelIdxList)
                    celldata(iCell).SpotPixelIdxList = {spotData(iSpot).PixelIdxList};
                else
                    celldata(iCell).SpotPixelIdxList{end + 1} = spotData(iSpot).PixelIdxList;
                end

                %Remove the spot indices from the indices used to calculate the
                %diffuse signal - this is probably not very efficient right now
                for ii = 1:numel(spotData(iSpot).PixelIdxList)
                    idxMatch = celldata(iCell).DiffusePixelIdxList == spotData(iSpot).PixelIdxList(ii);
                    celldata(iCell).DiffusePixelIdxList(idxMatch) = [];
                end

            end
        end

        %Calculate total intensity in cell vs total intensity in spot - need to
        %exclude regions in spots
        celldata(iCell).TotalIntDiffuse = sum(image(celldata(iCell).DiffusePixelIdxList));
        celldata(iCell).TotalIntSpot = 0;
        for iSpots = 1:numel(celldata(iCell).SpotPixelIdxList)
            celldata(iCell).TotalIntSpot = ...
                celldata(iCell).TotalIntSpot + sum(image(celldata(iCell).SpotPixelIdxList{iSpots}));
        end
    end

    %% Save data


    save(fullfile(outputDir, [outputFN, '.mat']), 'celldata', 'mask3', 'spotMask')

    fprintf('\b DONE\n')
end


%%









