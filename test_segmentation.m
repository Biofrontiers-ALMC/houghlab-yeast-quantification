clearvars
clc

file = 'D:\Work\ALMC\houghlab-yeast-quantification\data\yJM1837_Sup35MstartWT_richmedia002.nd2';

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
showoverlay(imadjust(imageData(:, :, 16)), bwperim(mask(:, :, 16)));

%% Find vacuoles

cellData = regionprops3(mask, imageData, 'MeanIntensity', 'VoxelIdxList', 'VoxelValues');

for iCell = 1:height(cellData)
    
    
end


%imshow(imadjust(currPlane), [])

