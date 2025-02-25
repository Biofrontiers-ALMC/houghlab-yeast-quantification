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

%Parse the inputs
if ischar(input) || isStringScalar(input)
    input = {input};
end

ip = inputParser;
ip.addParameter('CellSensitivity', 0.2);
parse(ip, varargin{:});

%inputDir = 'Z:\Microscopy\Yeast\Sup35\20250211_JA_Sup35MstartWT_4hrinduction';
%outputDir = 'Z:\Microscopy\Yeast\Sup35\20250214 Analysis JWT\MATLAB';
%outputDir = 'D:\Work\ALMC\houghlab-yeast-quantification\data\processed\20250225';

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

%% Process files

%Size of spot filter
sigma1 = 1/(1 + sqrt(2)) * 6;
sigma2 = 1/(1 + sqrt(2)) * 12;

if ~exist(outputDir, 'dir')
    mkdir(outputDir)
end

for iFile = 1:numel(fileList)

    %Print statement showing which file we're working on
    fprintf('[%s] Processing %s (file %d of %d)...\n', ...
        datetime, fileList{iFile})

    %Create a BioformatsImage object to read ND2 files
    reader = BioformatsImage(fileList{iFile});

    %Generate an output filename
    [~, outputFN] = fileparts(reader.filename);

    %Initialize variables to store 3D masks and images
    cellMask = zeros(reader.height, reader.width, reader.sizeZ, 'logical');
    spotMask = zeros(reader.height, reader.width, reader.sizeZ, 'logical');
    image = zeros(reader.height, reader.width, reader.sizeZ, 'uint16');
    
    %Read in the image
    for iZ = 1:reader.sizeZ
     
        I = getPlane(reader, iZ, 1, 1);
     
        image(:, :, iZ) = I;
    end
    


    
    % 
    % %Make the cell mask
    % for iZ = 1:reader.sizeZ
    % 
    %     I = getPlane(reader, iZ, 1, 1);
    % 
    %     image(:, :, iZ) = I;
    % 
    %     %Find cells
    %     currMask = imbinarize(I, 'adaptive', 'Sensitivity', 0.1);
    %     currMask = imopen(currMask, strel('disk', 3));
    %     currMask = imfill(currMask, 'holes');
    % 
    %     cellMask(:, :, iZ) = currMask;
    % 
    %     % %Find spots
    %     % df1 = imgaussfilt(image(:, :, iZ), sigma1);
    %     % df2 = imgaussfilt(image(:, :, iZ), sigma2);
    %     % 
    %     % DoG = df1 - df2;
    %     % 
    %     % currSpotMask = DoG > 200;
    %     % 
    %     % tmp_spotData = regionprops(currSpotMask, 'Circularity', 'PixelIdxList');
    %     % 
    %     % for iSpot = 1:numel(tmp_spotData)
    %     % 
    %     %     if tmp_spotData(iSpot).Circularity < 0.8
    %     % 
    %     %         currSpotMask(tmp_spotData(iSpot).PixelIdxList) = false;
    %     % 
    %     %     end
    %     % 
    %     % end
    %     % 
    %     % spotMask(:, :, iZ) = currSpotMask;
    % 
    % 
    %     %%imshow(mask)
    % end

    %Clean up the cell mask
    cellMask = imopen(cellMask, strel('sphere', 4));

    %%
    dd = -bwdist(~cellMask);
    dd(~cellMask) = Inf;

    dd = imhmin(dd, 0.8);

    LL = watershed(dd);
    LL(~cellMask) = 0;

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
    cellMask(LL == 0) = false;
    %volshow(mask3)


    %% Visualize by generating zStack of images
    for iZ = 1:size(image, 3)

        imgOut = showoverlay(image(:, :, iZ), bwperim(cellMask(:, :, iZ)), 'Color', [0 1 0]);
        imgOut = showoverlay(imgOut, bwperim(spotMask(:, :, iZ)), 'Color', [1 0 1]);

        if iZ == 1
            imwrite(imgOut, fullfile(outputDir, [outputFN, '.tiff']), 'Compression', 'none')
        else
            imwrite(imgOut, fullfile(outputDir, [outputFN, '.tiff']), 'Compression', 'none', 'WriteMode', 'append')
        end

    end


    celldata = regionprops(cellMask, image, 'PixelIdxList', 'PixelList', 'PixelValues');
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


    save(fullfile(outputDir, [outputFN, '.mat']), 'celldata', 'cellMask', 'spotMask')

    fprintf('\b DONE\n')
end


%%









