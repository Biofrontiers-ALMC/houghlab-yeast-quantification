clearvars
close all
clc

folder = 'Y:\Microscopy\Yeast\Sup35\20250321 Analysis\20250211_JA_Sup35MstartWT_4hrinduction\richmedia';
outputDir = 'Y:\Microscopy\Yeast\Sup35\20250321 Analysis\Plots\20250211_JA_Sup35MstartWT_4hrinduction';

[~, fpath] = fileparts(folder);
if strcmpi(fpath, 'edmedia')
    isEDmedia = true;
else
    isEDmedia = false;
end


if ~exist(outputDir, 'dir')
    mkdir(outputDir)
end


files = dir(fullfile(folder, '*.mat'));

storeNetCellMeanIntensity = [];
storeDiffuseMeanIntensity = [];
storeSpotMeanIntensity = [];
storeCreationTime = [];

for iF = 1:numel(files)

    load(fullfile(folder, files(iF).name));

    %There was an error in earlier versions so reacquire the modification
    %date
    % if isempty(storeCreationDate)
    %     storeCreationDate = {finalCellData.CreationDate};
    % else
    %     % storeCreationDate{(end + 1):(end + numel(finalCellData))} = ...
    %     %     deal({finalCellData.CreationDate});
    %     storeCreationDate = [storeCreationDate, {finalCellData.CreationDate}];
    % end

    d = System.IO.File.GetLastWriteTime(finalCellData(1).Filename);

    ModifiedDateTime = datetime(d.Year, d.Month, d.Day, d.Hour, d.Minute, d.Second);
    PosixTime = convertTo(ModifiedDateTime, 'posixTime');

    storeCreationTime((end + 1):(end + numel(finalCellData))) = ...
        repmat(PosixTime, 1, numel(finalCellData));

    %Collate the data we want to plot
    storeNetCellMeanIntensity((end + 1):(end + numel(finalCellData))) = ...
        cat(2, finalCellData.NetCellMeanIntensity);

    storeDiffuseMeanIntensity((end + 1):(end + numel(finalCellData))) = ...
        cat(2, finalCellData.DiffuseMeanIntensity);

     storeSpotMeanIntensity((end + 1):(end + numel(finalCellData))) = ...
        cat(2, finalCellData.SpotMeanIntensity);

end

if isEDmedia
    %keyboard

    uniqueDT = unique(storeCreationTime);

    diffT = uniqueDT - min(uniqueDT);
    idxStart = find(diffT > 600, 1, 'first');

    storeNetCellMeanIntensity(storeCreationTime < uniqueDT(idxStart)) = [];

    storeDiffuseMeanIntensity(storeCreationTime < uniqueDT(idxStart)) = [];

    storeSpotMeanIntensity(storeCreationTime < uniqueDT(idxStart)) = [];

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