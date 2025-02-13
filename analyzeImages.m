clearvars
clc

reader = BioformatsImage('../data/yJM1837_Sup35WT_4_5hrs_EDmedia019.nd2');

mask = zeros(reader.height, reader.width, reader.sizeZ, 'logical');
spotMask = zeros(reader.height, reader.width, reader.sizeZ, 'logical');
image = zeros(reader.height, reader.width, reader.sizeZ, 'uint16');

for iZ = 1:reader.sizeZ

    I = getPlane(reader, iZ, 1, 1);

    image(:, :, iZ) = I;
    %%
    currMask = imbinarize(I, 'adaptive');
    currMask = imopen(currMask, strel('disk', 3));
    currMask = imfill(currMask, 'holes');

    mask(:, :, iZ) = currMask;

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
%% Spot finding
sigma1 = 2.8;
sigma2 = 2;

df1 = imgaussfilt3(image, sigma1);
df2 = imgaussfilt3(image, sigma2);

DoG = df1 - df2;

spotMask = DoG > 50;

%%
largeObjs = bwareaopen(spotMask, 100, 26);

spotMask(largeObjs) = false;

volshow(spotMask)

%%
mask3(LL == 0) = false;
volshow(mask3)

%%
celldata = regionprops(mask3, image, 'PixelIdxList', 'PixelList', 'PixelValues');
spotData = regionprops(spotMask, image, 'Centroid', 'MeanIntensity', 'PixelValues', 'PixelIdxList');

spotCentroids = cat(1, spotData.Centroid);

%%
for iCell = 1:numel(celldata)

    celldata(iCell).NumSpots = 0;
    celldata(iCell).TotalSpotInt = 0;
    celldata(iCell).SpotPixelIdxList = {};

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
        end
    end

    %Calculate total intensity in cell vs total intensity in spot - need to
    %exclude regions in spots
    


    
end

%%




%% Visualize in 





