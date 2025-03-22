function analyzeImages(input, outputDir, varargin)
%ANALYZEIMAGES  Analyze images of fluorescent yeast
%
%  ANALYZEIMAGES(INPUT) will process image(s) of fluorescent yeast,
%  attempting to quantify protein concentration in both diffuse and
%  punctate form. INPUT can be a string specifying a single file or
%  directory, or a list of files or directories to process. Output files
%  will be saved in the current working directory.
%
%  ANALYZEIMAGES(INPUT, OUTPUTDIR) will alternatively output files in the
%  directory specified. If the path in OUTPUTDIR does not exist, it will be
%  created.
%
%  ANALYZEIMAGES(..., 'Property', 'Value') property/value pairs can be used
%  to change default parameters. The following properties are currently
%  defined:
%    * 'CellSensitivity' (Default: 0.2) - Higher numbers mean more pixels
%      considered cells
%
%  Example:
%
%  ANALYZEIMAGES('Z:\Microscopy\Yeast\Sup35\20250211_JA_Sup35MstartWT_4hrinduction', ...
%     'Z:\Microscopy\Yeast\Sup35\20250214 Analysis JWT\MATLAB');

%% Parse inputs
if ischar(input) || isStringScalar(input)
    input = {input};
end

ip = inputParser;
ip.addParameter('Smoothing', 3);
ip.addParameter('CellSensitivity', 0.08);
ip.addParameter('SpotSensitivity', 1.75);
parse(ip, varargin{:});

fileList = {};

fileCtr = 0;
for iInputs = 1:numel(input)

    if isfile(input(iInputs))
        fileCtr = fileCtr + 1;
        fileList{fileCtr} = input{iInputs}; %#ok<AGROW>

    elseif isfolder(input(iInputs))

        files = dir(fullfile(input{iInputs}, '*.nd2'));

        for iFile = 1:numel(files)
            fileCtr = fileCtr + 1;
            fileList{fileCtr} = fullfile(files(iFile).folder, files(iFile).name);            %#ok<AGROW>
        end

    else
        error('analyzeImages:UnknownPath', ...
            '''%s'' is not a valid file or folder path.', input{iInputs});
    end
end

if ~exist('outputDir', 'var')
    outputDir = '';
end

if ~exist(outputDir, 'dir')
    mkdir(outputDir)
end

%% Process files

for iFile = 1:numel(fileList)

    %Get name and folder current file
    [fpath, fn, fext] = fileparts(fileList{iFile});
    [~, fdir] = fileparts(fpath);

    %Print working statement
    fprintf('[%s] Processing %s (file %d of %d)...\n', ...
        datetime, fullfile(fdir, [fn, fext]), iFile, numel(fileList));

    %% Process individual file

    %Create a BioformatsImage object to read the file
    reader = BioformatsImage(fileList{iFile});

    %Create a matrix to store image data
    imageData = zeros(reader.height, reader.width, reader.sizeZ, 'uint16');

    %Read in 3D image data
    for iZ = 1:reader.sizeZ

        imageData(:, :, iZ) = getPlane(reader, iZ, 1, 1);

    end

    %Calculate a median filtered version for segmentation only
    imageDataSmoothed = medfilt3(imageData, ...
        [ip.Results.Smoothing, ip.Results.Smoothing, ip.Results.Smoothing]);

    %Calculate a cell threshold intensity
    bgLvl = mode(imageDataSmoothed, 'all');
    thLvl = bgLvl + ip.Results.CellSensitivity * bgLvl;

    %Make the cell mask
    mask = imageDataSmoothed > thLvl;

    mask = imopen(mask, strel('sphere', 4));
    mask = imfill(mask, 4, 'holes');

    dd = -bwdist(~mask);
    dd(~mask) = Inf;
    dd = imhmin(dd, 1);

    LL = watershed(dd);

    mask(LL == 0) = 0;

    mask = imerode(mask, strel('sphere', 2));
    mask = imclearborder(mask, 4);

    %showoverlay(imadjust(imageData(:, :, 16)), bwperim(mask(:, :, 16)));

    %---Find vacuoles and bright regions---

    cellData = regionprops3(mask, imageDataSmoothed, 'MeanIntensity', 'VoxelIdxList', 'VoxelValues', 'Centroid');

    for iCell = 1:height(cellData)

        %Calculate the vacuole threshold
        vacTh = 1.2 * cellData(iCell, :).MeanIntensity;

        %Find voxels which are below this threshold
        isVacuole = cellData(iCell, :).VoxelValues{:} < vacTh;

        % %Update the mask
        % vacIdxList = cellData(iCell, :).VoxelIdxList{:};
        % vacIdxList = vacIdxList(isVacuole);
        % vacMask(vacIdxList) = true;

        %Recalculate the mean cell intensity for bright region detection
        vxVal = cellData(iCell, :).VoxelValues{:};
        netMeanIntensity = mean(vxVal(~isVacuole));  %Note that this is intensity of smoothed image
        spotTh = ip.Results.SpotSensitivity * netMeanIntensity;

        %Find voxels which are below this threshold
        isSpot = cellData(iCell, :).VoxelValues{:} > spotTh;

        % %Update the mask
        % spotIdxList = cellData(iCell, :).VoxelIdxList{:};
        % spotIdxList = spotIdxList(isSpot);
        % spotMask(spotIdxList) = true;

        %---Put together the final cell data---

        %Include metadata about each cell
        finalCellData(iCell).Filename = reader.filename;
        d = System.IO.File.GetLastWriteTime(reader.filename);
        finalCellData(iCell).CreationDate = datetime(d.Year, d.Month, d.Day, d.Hour, d.Minute, d.Second);
        finalCellData(iCell).PosixTime = convertTo(finalCellData(iCell).CreationDate, 'posixTime');

        %Measure data from full cell
        %Note: This is larger than the final cell because some of the
        %surrounding region is classified as "vacuole"
        finalCellData(iCell).Centroid = cellData(iCell, :).Centroid;
        finalCellData(iCell).RawCellIdxList = cellData(iCell, :).VoxelIdxList{:};
        finalCellData(iCell).RawCellVolume = numel(finalCellData(iCell).RawCellIdxList);
        finalCellData(iCell).RawCellMeanIntensity = mean(imageData(finalCellData(iCell).RawCellIdxList), 'all');

        %Measure diffuse intensity (excl. vacuoles and spots)
        isDiffuse = ~isVacuole & ~isSpot;
        finalCellData(iCell).DiffuseIdxList = finalCellData(iCell).RawCellIdxList(isDiffuse);
        finalCellData(iCell).DiffuseVolume = numel(finalCellData(iCell).DiffuseIdxList);
        finalCellData(iCell).DiffuseMeanIntensity = mean(imageData(finalCellData(iCell).DiffuseIdxList), 'all');

        %Measure spot intensities (excl. vacuoles in case of overlap)
        isOnlySpot = ~isVacuole & isSpot;
        finalCellData(iCell).SpotIdxList = finalCellData(iCell).RawCellIdxList(isOnlySpot);
        finalCellData(iCell).SpotVolume = numel(finalCellData(iCell).SpotIdxList);
        finalCellData(iCell).SpotMeanIntensity = mean(imageData(finalCellData(iCell).SpotIdxList), 'all');

        %Measure intensity within diffuse region and spots
        isDiffuseAndSpot = isDiffuse | isOnlySpot;
        finalCellData(iCell).NetCellIdxList = finalCellData(iCell).RawCellIdxList(isDiffuseAndSpot);
        finalCellData(iCell).NetCellVolume = numel(finalCellData(iCell).NetCellIdxList);
        finalCellData(iCell).NetCellMeanIntensity = mean(imageData(finalCellData(iCell).NetCellIdxList), 'all');

    end

    %---Generate outputs---%

    %Renormalize the image data for output image
    imageDataNorm = double(imageData);
    imageDataNorm = (imageDataNorm - min(imageDataNorm, [], 'all')) ./ (max(imageDataNorm, [], 'all') - min(imageDataNorm, [], 'all'));

    %Remake the masks to double-check
    finalFullCellMask = false(size(imageDataNorm));
    finalDiffuseMask = false(size(imageDataNorm));
    finalSpotMask = false(size(imageDataNorm));

    %Save data
    save(fullfile(outputDir, [fn, '.mat']), 'finalCellData', 'finalFullCellMask', 'finalDiffuseMask', 'finalSpotMask');

    for iCell = 1:numel(finalCellData)

        finalFullCellMask(finalCellData(iCell).RawCellIdxList) = true;
        finalDiffuseMask(finalCellData(iCell).DiffuseIdxList) = true;
        finalSpotMask(finalCellData(iCell).SpotIdxList) = true;

    end

    for iZ = 1:size(imageData, 3)

        imgOut = showoverlay(imageDataNorm(:, :, iZ), bwperim(finalFullCellMask(:, :, iZ)), 'Color', [0 1 0]);

        if any(finalDiffuseMask(:, :, iZ), 'all')
            imgOut = showoverlay(imgOut, finalDiffuseMask(:, :, iZ), 'Color', [1 1 0], 'Opacity', 20);
        end
        if any(finalSpotMask(:, :, iZ), 'all')
            imgOut = showoverlay(imgOut, bwperim(finalSpotMask(:, :, iZ)), 'Color', [1 0 1]);
        end

        if iZ == 1
            imwrite(imgOut, fullfile(outputDir, [fn, '.tiff']), 'Compression', 'none')
        else
            imwrite(imgOut, fullfile(outputDir, [fn, '.tiff']), 'Compression', 'none', 'WriteMode', 'append')
        end

    end

    clearvars finalCellData

end




% 
% %% Process files
% 
% %Size of spot filter
% sigma1 = 1/(1 + sqrt(2)) * 6;
% sigma2 = 1/(1 + sqrt(2)) * 12;
% 
% if ~exist(outputDir, 'dir')
%     mkdir(outputDir)
% end
% 
% for iFile = 1%:numel(fileList)
% 
%     %Print statement showing which file we're working on
%     fprintf('[%s] Processing %s (file %d of %d)...\n', ...
%         datetime, fileList{iFile})
% 
%     %Create a BioformatsImage object to read ND2 files
%     reader = BioformatsImage(fileList{iFile});
% 
%     %Generate an output filename
%     [~, outputFN] = fileparts(reader.filename);
% 
%     %Initialize variables to store 3D masks and images
%     %cellMask = zeros(reader.height, reader.width, reader.sizeZ, 'logical');
%     spotMask = zeros(reader.height, reader.width, reader.sizeZ, 'logical');
%     image = zeros(reader.height, reader.width, reader.sizeZ, 'uint16');
% 
%     %Read in the image
%     for iZ = 1:reader.sizeZ
% 
%         image(:, :, iZ) = getPlane(reader, iZ, 1, 1);
% 
%     end
% 
%     imageDataSmoothed = medfilt3(imageData, [3 3 3]);
% 
%     bgLvl = mode(imageDataSmoothed, 'all');
%     %%
%     thLvl = bgLvl + 0.07 * bgLvl;
% 
%     mask = imageDataSmoothed > thLvl;
% 
%     mask = imopen(mask, strel('sphere', 4));
%     mask = imfill(mask, 4, 'holes');
% 
%     dd = -bwdist(~mask);
%     dd(~mask) = Inf;
% 
%     dd = imhmin(dd, 1);
% 
%     LL = watershed(dd);
% 
%     mask(LL == 0) = 0;
%     % 
%     % %Create a mask of the cells using a smoothed image to remove unwanted
%     % %hot spots
%     % cellMask = imbinarize(medfilt3(image, [ip.Results.Smoothing, ip.Results.Smoothing, ip.Results.Smoothing]),...
%     %      'adaptive', 'Sensitivity', ip.Results.CellSensitivity);
%     % 
%     % %Clean up the cell mask
%     % cellMask = imopen(cellMask, strel('sphere', 4));
%     % cellMask = imclose(cellMask, strel('sphere', 5));
%     % 
%     % cellMask = imfill(cellMask, 26, 'holes');
%     % 
%     % %Watershed to separate clusters of objects
%     % dd = -bwdist(~cellMask);
%     % dd(~cellMask) = Inf;
%     % 
%     % dd = imhmin(dd, 7, 26);
%     % 
%     % LL = watershed(dd, 26);
%     % LL(~cellMask) = 0;
%     % 
%     % cellMask(LL == 0) = false;
%     % 
%     % cellMask = bwareaopen(cellMask, 5000, 26);
%     % % volshow(cellMask)
%     % 
%     % %--Measure cell data--
%     % cellData = regionprops3(cellMask, image, 'Volume', 'VoxelValues', 'VoxelIdxList');
%     % 
%     % for iCell = 1:height(cellData)
%     % 
%     %     %currCellMask = cellData(iCell, :).Image{:};
%     % 
%     %     currVoxelValues = cellData(iCell, :).VoxelValues;
%     %     meanVoxelValue = mean(currVoxelValues{:});
%     %     stdVoxelValue = std(double(currVoxelValues{:}), 0, 'all');
%     % 
%     %     spotTh = meanVoxelValue + 1.5 * stdVoxelValue;
%     %     vacuoleTh = meanVoxelValue - stdVoxelValue;
%     % 
%     %     %currSpotMask = cellData(iCell, :).VoxelValues > spotTh;
%     %     %currCellMask((cellData(iCell, :).VoxelValues{:}) < vacuoleTh) = false;
%     %     % 
%     %     % volshow(currCellMask);
%     %     % pause
%     % 
%     %     try
%     %         isVacuole = (cellData(iCell, :).VoxelValues{:}) < vacuoleTh;
%     %         currVoxelIdxList = cellData(iCell, :).VoxelIdxList{:};
%     % 
%     %         if any(isVacuole)
%     %             cellMask(currVoxelIdxList(isVacuole)) = false;
%     %         end
%     %     catch
%     %         keyboard
%     %     end
%     % end
%     % 
%     % 
% 
% 
% 
%     % 
%     % %Make the cell mask
%     % for iZ = 1:reader.sizeZ
%     % 
%     %     I = getPlane(reader, iZ, 1, 1);
%     % 
%     %     image(:, :, iZ) = I;
%     % 
%     %     %Find cells
%     %     currMask = imbinarize(I, 'adaptive', 'Sensitivity', 0.1);
%     %     currMask = imopen(currMask, strel('disk', 3));
%     %     currMask = imfill(currMask, 'holes');
%     % 
%     %     cellMask(:, :, iZ) = currMask;
%     % 
%     %     % %Find spots
%     %     % df1 = imgaussfilt(image(:, :, iZ), sigma1);
%     %     % df2 = imgaussfilt(image(:, :, iZ), sigma2);
%     %     % 
%     %     % DoG = df1 - df2;
%     %     % 
%     %     % currSpotMask = DoG > 200;
%     %     % 
%     %     % tmp_spotData = regionprops(currSpotMask, 'Circularity', 'PixelIdxList');
%     %     % 
%     %     % for iSpot = 1:numel(tmp_spotData)
%     %     % 
%     %     %     if tmp_spotData(iSpot).Circularity < 0.8
%     %     % 
%     %     %         currSpotMask(tmp_spotData(iSpot).PixelIdxList) = false;
%     %     % 
%     %     %     end
%     %     % 
%     %     % end
%     %     % 
%     %     % spotMask(:, :, iZ) = currSpotMask;
%     % 
%     % 
%     %     %%imshow(mask)
%     % end
% 
%     % %% Spot finding
%     % sigma1 = 6;
%     % sigma2 = 3;
%     % 
%     % df1 = imgaussfilt3(image, sigma1);
%     % df2 = imgaussfilt3(image, sigma2);
%     % 
%     % DoG = df2 - df1;
%     % 
%     % spotMask = DoG > 250;
% 
% 
%     %%
%     % largeObjs = bwareaopen(spotMask, 500, 26);
%     % 
%     % spotMask(largeObjs) = false;
% 
%     %volshow(spotMask)
% 
%     %%
% 
% 
%     %% Visualize by generating zStack of images
%     for iZ = 1:size(image, 3)
% 
%         imgOut = showoverlay(image(:, :, iZ), bwperim(cellMask(:, :, iZ)), 'Color', [0 1 0]);
%         %imgOut = showoverlay(imgOut, bwperim(spotMask(:, :, iZ)), 'Color', [1 0 1]);
% 
%         if iZ == 1
%             imwrite(imgOut, fullfile(outputDir, [outputFN, '.tiff']), 'Compression', 'none')
%         else
%             imwrite(imgOut, fullfile(outputDir, [outputFN, '.tiff']), 'Compression', 'none', 'WriteMode', 'append')
%         end
% 
%     end
% 
%     % 
%     % 
%     % celldata = regionprops(cellMask, image, 'PixelIdxList', 'PixelList', 'PixelValues');
%     % spotData = regionprops(spotMask, image, 'Centroid', 'MeanIntensity', 'PixelValues', 'PixelIdxList');
%     % 
%     % spotCentroids = cat(1, spotData.Centroid);
%     % 
%     % %%
%     % for iCell = 1:numel(celldata)
%     % 
%     %     celldata(iCell).NumSpots = 0;
%     %     celldata(iCell).TotalSpotInt = 0;
%     %     celldata(iCell).SpotPixelIdxList = {};
%     %     celldata(iCell).DiffusePixelIdxList = celldata(iCell).PixelIdxList;
%     % 
%     %     for iSpot = 1:numel(spotData)
%     %         %Check if spot is in cell
%     %         if ismember(round(spotCentroids(iSpot, :)), celldata(iCell).PixelList, 'rows')
%     % 
%     %             %Store spot information as well
%     % 
%     %             celldata(iCell).NumSpots = celldata(iCell).NumSpots + 1;
%     %             celldata(iCell).TotalSpotInt = celldata(iCell).TotalSpotInt + ...
%     %                 sum(spotData(iSpot).PixelValues);
%     % 
%     %             if isempty(celldata(iCell).SpotPixelIdxList)
%     %                 celldata(iCell).SpotPixelIdxList = {spotData(iSpot).PixelIdxList};
%     %             else
%     %                 celldata(iCell).SpotPixelIdxList{end + 1} = spotData(iSpot).PixelIdxList;
%     %             end
%     % 
%     %             %Remove the spot indices from the indices used to calculate the
%     %             %diffuse signal - this is probably not very efficient right now
%     %             for ii = 1:numel(spotData(iSpot).PixelIdxList)
%     %                 idxMatch = celldata(iCell).DiffusePixelIdxList == spotData(iSpot).PixelIdxList(ii);
%     %                 celldata(iCell).DiffusePixelIdxList(idxMatch) = [];
%     %             end
%     % 
%     %         end
%     %     end
%     % 
%     %     %Calculate total intensity in cell vs total intensity in spot - need to
%     %     %exclude regions in spots
%     %     celldata(iCell).TotalIntDiffuse = sum(image(celldata(iCell).DiffusePixelIdxList));
%     %     celldata(iCell).TotalIntSpot = 0;
%     %     for iSpots = 1:numel(celldata(iCell).SpotPixelIdxList)
%     %         celldata(iCell).TotalIntSpot = ...
%     %             celldata(iCell).TotalIntSpot + sum(image(celldata(iCell).SpotPixelIdxList{iSpots}));
%     %     end
%     % end
%     % 
%     % %% Save data
%     % 
%     % 
%     % save(fullfile(outputDir, [outputFN, '.mat']), 'celldata', 'cellMask', 'spotMask')
% 
%     fprintf('\b DONE\n')
% end
% 
% 
% %%
% 
% 
% 
% 





