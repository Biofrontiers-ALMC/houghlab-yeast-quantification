clearvars
clc

folder = 'Z:\Microscopy\Yeast\Sup35\20250305 Analysis\20241120_JA_Sup35Mstart_2hrs_induction';

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
