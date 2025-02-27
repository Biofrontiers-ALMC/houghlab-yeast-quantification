clearvars
clc

%file = 'D:\Work\ALMC\houghlab-yeast-quantification\data\yJM1837_Sup35MstartWT_richmedia002.nd2';
file = '../data/yJM1837_Sup35WT_4_5hrs_EDmedia021.nd2';

reader = BioformatsImage(file);

imageData = zeros(reader.height, reader.width, reader.sizeZ, 'uint16');

for iZ = 1:reader.sizeZ

    imageData(:, :, iZ) = getPlane(reader, iZ, 1, 1);

end

%%

%Try and calculate an acceptable threshold
%currPlane = medfilt2(imageData(:, :, 1), [3 3]);

imageDataSmoothed = medfilt3(imageData, [3 3 3]);

bgLvl = mode(imageDataSmoothed, 'all');
%%
thLvl = bgLvl + 0.07 * bgLvl;

mask = imageDataSmoothed > thLvl;

mask = imopen(mask, strel('sphere', 4));
mask = imfill(mask, 4, 'holes');

dd = -bwdist(~mask);
dd(~mask) = Inf;

dd = imhmin(dd, 1);

LL = watershed(dd);

mask(LL == 0) = 0;

mask = imerode(mask, strel('sphere', 3));
mask = imclearborder(mask, 4);

%showoverlay(imadjust(imageData(:, :, 16)), bwperim(mask(:, :, 16)));

%% Find vacuoles
%meanCellInt = mean(imageData(mask), 'all');

%vacMask = imageDataSmoothed <= (0.9 * meanCellInt);

vacMask = false(size(mask));
spotMask = false(size(mask));

cellData = regionprops3(mask, imageDataSmoothed, 'MeanIntensity', 'VoxelIdxList', 'VoxelValues');

for iCell = 1:height(cellData)

    %Calculate the vacuole threshold
    vacTh = 1.2 * cellData(iCell, :).MeanIntensity;

    %Find voxels which are below this threshold
    isVacuole = cellData(iCell, :).VoxelValues{:} < vacTh;

    %Update the mask
    idxList = cellData(iCell, :).VoxelIdxList{:};
    idxList = idxList(isVacuole);
    vacMask(idxList) = true;

    %Recalculate the mean cell intensity for spot detection
    vxVal = cellData(iCell, :).VoxelValues{:};
    netMeanIntensity = mean(vxVal(~isVacuole));
    spotTh = 4 * cellData(iCell, :).MeanIntensity; %Should modify to exclude vacuoles?

    %Find voxels which are below this threshold
    isSpot = cellData(iCell, :).VoxelValues{:} > spotTh;

    %Update the mask
    idxList = cellData(iCell, :).VoxelIdxList{:};
    idxList = idxList(isSpot);
    spotMask(idxList) = true;    

end
%showoverlay(imadjust(imageData(:, :, 16)), vacMask(:, :, 16));
%%
% plane = 18;
% showoverlay(imageData(:, :, plane), bwperim(spotMask(:, :, plane)))
% 
% %% Retry with spot detection
% 
% %Size of spot filter
% sigma1 = 1/(1 + sqrt(2)) * 5;
% sigma2 = 1/(1 + sqrt(2)) * 12;
% 
% dogSpotMask = false(size(mask));
% 
% for iZ = 1:size(imageData, 3)
% 
%     %Find spots
%     df1 = imgaussfilt(imageDataSmoothed(:, :, iZ), sigma1);
%     df2 = imgaussfilt(imageDataSmoothed(:, :, iZ), sigma2);
% 
%     DoG = df1 - df2;
% 
%     dogSpotMask(:, :, iZ) = DoG > 200;
% 
% end
% 
% %%
% plane = 30;
% showoverlay(imageData(:, :, plane), bwperim(dogSpotMask(:, :, plane)))



%%
cellMask = bwlabeln(mask);
cellMask(vacMask) = 0;

outputDir = '../processed/20250227';
outputFN = 'tmp';

%Renormalize the image data for output image
imageDataNorm = double(imageData);
imageDataNorm = (imageData - min(imageDataNorm, [], 'all')) / (max(imageDataNorm, [], 'all') - min(imageDataNorm, [], 'all'));

for iZ = 1:size(imageData, 3)

    imgOut = showoverlay(imageDataNorm(:, :, iZ), bwperim(cellMask(:, :, iZ)), 'Color', [0 1 0]);
    imgOut = showoverlay(imgOut, bwperim(spotMask(:, :, iZ)), 'Color', [1 0 1]);

    if iZ == 1
        imwrite(imgOut, fullfile(outputDir, [outputFN, '.tiff']), 'Compression', 'none')
    else
        imwrite(imgOut, fullfile(outputDir, [outputFN, '.tiff']), 'Compression', 'none', 'WriteMode', 'append')
    end

end



%imshow(imadjust(currPlane), [])

