clearvars
close all
clc

folder = 'Z:\Microscopy\Yeast\Sup35\20250305 Analysis\20250211_JA_Sup35MstartWT_4hrinduction';
outputDir = 'Z:\Microscopy\Yeast\Sup35\20250306 Plots';
files = dir(fullfile(folder, '*.mat'));

storeNetCellMeanIntensity = [];
storeDiffuseMeanIntensity = [];
storeSpotMeanIntensity = [];

for iF = 1:numel(files)

    load(fullfile(folder, files(iF).name));

    %Collate the data we want to plot
    storeNetCellMeanIntensity((end + 1):(end + numel(finalCellData))) = ...
        cat(2, finalCellData.NetCellMeanIntensity);

    storeDiffuseMeanIntensity((end + 1):(end + numel(finalCellData))) = ...
        cat(2, finalCellData.DiffuseMeanIntensity);

     storeSpotMeanIntensity((end + 1):(end + numel(finalCellData))) = ...
        cat(2, finalCellData.SpotMeanIntensity);

end

%%
scatter(storeNetCellMeanIntensity, storeDiffuseMeanIntensity)
hold on
scatter(storeNetCellMeanIntensity, storeSpotMeanIntensity)
hold off
xlabel('Cell Mean Intensity')
ylabel('Diffuse/Spot Mean Intensity')
legend('Diffuse', 'Spot')

[~, dataSetName] = fileparts(folder);
title(dataSetName, 'Interpreter', 'none');

saveas(gcf, fullfile(outputDir, [dataSetName, '.fig']));
saveas(gcf, fullfile(outputDir, [dataSetName, '.jpg']));