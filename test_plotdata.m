clearvars
clc

load('Z:\Microscopy\Yeast\Sup35\20250214 Analysis JWT\MATLAB\yJM1837_Sup35WT_4_5hrs_EDmedia015.mat');

numCells = numel(celldata);

totalProtein = zeros(numCells, 1);
totalSpot = totalProtein;
totalDiffuse = totalProtein;

%Make a bar plot showing total protein concentration in this image
for iCell = 1:numel(celldata)

    totalProtein(iCell) = (celldata(iCell).TotalIntDiffuse + celldata(iCell).TotalIntSpot)/...
        numel(celldata(iCell).PixelIdxList);
    
    totalDiffuse(iCell) = celldata(iCell).TotalIntDiffuse/numel(celldata(iCell).DiffusePixelIdxList);

    totalSpot(iCell) = celldata(iCell).TotalIntSpot;
    spotVol = 0;
    for ii = 1:numel(celldata(iCell).SpotPixelIdxList)

        spotVol = spotVol + numel(celldata(iCell).SpotPixelIdxList{ii});

    end

    totalSpotConc(iCell) = totalSpot(iCell) / spotVol;

end

%%
figure;

scatter(totalProtein, totalDiffuse)
hold on
scatter(totalProtein, totalSpotConc)

% 
% stackedData = cat(2, totalDiffuse, totalSpot);
% 
% bar(stackedData, 'stacked')
% ylabel('Total intensity')
% xlabel('Cell')
legend('Diffuse', 'Punctate')