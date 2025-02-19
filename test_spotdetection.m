clearvars
clc

file = 'Z:\Microscopy\Yeast\Sup35\20250211_JA_Sup35MstartWT_4hrinduction\yJM1837_Sup35WT_4_5hrs_EDmedia017.nd2';
reader = BioformatsImage(file);

image = zeros(reader.height, reader.width, reader.sizeZ, 'uint16');
spotMask = zeros(reader.height, reader.width, reader.sizeZ, 'logical');

sigma1 = 1/(1 + sqrt(2)) * 6;
sigma2 = 1/(1 + sqrt(2)) * 12;


    % sigma1 = 6;
    % sigma2 = 3;

for iZ = 1:reader.sizeZ

    image(:, :, iZ) = getPlane(reader, iZ, 1, 1);

    df1 = imgaussfilt(image(:, :, iZ), sigma1);
    df2 = imgaussfilt(image(:, :, iZ), sigma2);

    DoG = df1 - df2;

    cspotmask = DoG > 200;

    tmp_spotData = regionprops(cspotmask, 'Circularity', 'PixelIdxList');

    for iSpot = 1:numel(tmp_spotData)

        if tmp_spotData(iSpot).Circularity < 0.8

            cspotmask(tmp_spotData(iSpot).PixelIdxList) = false;

        end

    end

    spotMask(:, :, iZ) = cspotmask;

end

%% Spot finding
% sigma1 = 6;
% sigma2 = 3;

% df1 = imgaussfilt3(image, sigma1);
% df2 = imgaussfilt3(image, sigma2);
% 
% DoG = df1 - df2;
% 
% spotMask = DoG > 200;

for iZ = 1:reader.sizeZ

    imgOut = showoverlay(image(:, :, iZ), bwperim(spotMask(:, :, iZ)), 'Color', [1 0 1]);

    if iZ == 1
        imwrite(imgOut, 'test.tiff', 'Compression', 'none')
    else
        imwrite(imgOut, 'test.tiff', 'Compression', 'none', 'WriteMode', 'append')
    end

end

%% For single plane
return


I = getPlane(reader, iZ, 1, 1);
sigma1 = 5;
sigma2 = 2;

df1 = imgaussfilt(I, sigma1);
df2 = imgaussfilt(I, sigma2);

DoG = df2 - df1;

spotMask = DoG > 250;

figure(1);
subplot(1, 2, 1)
showoverlay(I, spotMask)
subplot(1, 2, 2)
imshow(I, [])