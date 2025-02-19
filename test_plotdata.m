clearvars
clc

load('Z:\Microscopy\Yeast\Sup35\20250214 Analysis JWT\MATLAB\yJM1837_Sup35WT_4_5hrs_EDmedia015.mat');

numCells = numel(celldata);

totalProtein = zeros(numCells, 1);
totalSpot = totalProtein;
totalDiffuse = totalProtein;

%Make a bar plot showing total protein concentration in this image
for iCell = 1:numel(celldata)

    totalProtein(iCell) = celldata(iCell).TotalIntDiffuse + celldata(iCell).TotalIntSpot;
    totalSpot(iCell) = celldata(iCell).TotalIntSpot;
    totalDiffuse(iCell) = celldata(iCell).TotalIntDiffuse;

end

%%
stackedData = cat(2, totalDiffuse, totalSpot);

bar(stackedData, 'stacked')
ylabel('Total intensity')
xlabel('Cell')
legend('Diffuse', 'Punctate')