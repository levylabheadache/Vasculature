clear; clc; close all;
dataDir = 'D:\2photon\'; % 'D:\2photon\Simone\'; %'C:\2photon';
dataSet = 'Vasculature'; %'Afferents'; %  'Neutrophil_Simone'; %'Pollen'; %'Astrocyte'; %    'NGC'; % 'Neutrophil'; %
% Parse data table
dataTablePath = 'R:\Levy Lab\2photon\ImagingDatasets.xlsx'; % 'R:\Levy Lab\2photon\ImagingDatasetsSimone2.xlsx'; %'D:\MATLAB\ImagingDatasets.xlsx'; % 'D:\MATLAB\NGCdata.xlsx';  Simone
dataTable = readcell(dataTablePath, 'sheet',dataSet);  % 'NGC', ''
colNames = dataTable(1,:); dataTable(1,:) = [];
dataCol = struct('mouse',find(contains(colNames, 'Mouse')), 'date',find(contains(colNames, 'Date')), 'FOV',find(contains(colNames, 'FOV')), 'vasc',find(contains(colNames, 'vascChan')),...
    'volume',find(contains(colNames, 'Volume')), 'run',find(contains(colNames, 'Runs')), 'Ztop',find(contains(colNames, 'Zbot')), 'Zbot',find(contains(colNames, 'Ztop')), 'csd',find(contains(colNames, 'CSD')), 'ref',find(contains(colNames, 'Ref')), 'done',find(contains(colNames, 'Done')));
Nexpt = size(dataTable, 1);
dataTable(:,dataCol.date) = cellfun(@num2str, dataTable(:,dataCol.date), 'UniformOutput',false);

% Initialize variables
expt = cell(1,Nexpt); runInfo = cell(1,Nexpt); Tscan = cell(1,Nexpt); loco = cell(1,Nexpt); %vessels(x) = cell(1,Nexpt);
vesselROI = cell(1,Nexpt); NvesselROI = cell(1,Nexpt); tifStackMax = cell(1,Nexpt);
%vessels = repmat(struct('notes',[], 'data',[]), 1, Nexpt);
roiTemplate = struct( 'boxPosition',[], 'xSize',[], 'ySize',NaN, 'vesselLine',[], 'projectionAngle',[], 'projection',NaN, 'modalFixedDiameter',[], 'vesselType','' );

% Choose which subset to  process
xPresent = 3; % flip(100:102); %45:47; % [66:69]; %6;  62,64,
Npresent = numel(xPresent);
overwrite = false;
FS = 12;
for x = xPresent  %30 %x2D % x2Dcsd % x3D %% 51
    % Experiment/run level metadata
    [expt{x}, runInfo{x}, regParams, projParam] = ParseDataTable(dataTable, x, dataCol, dataDir); % , regParams, projParam
    expt{x}.vascChan = 'red';
    
    catInfo(x) = ConcatenateRunInfo(expt{x}, runInfo{x}, 'overwrite',false);  % 'suffix','sbxcat',
    [Tscan{x}, runInfo{x}] = GetTime(runInfo{x});
    projParam = GenerateExptProjections(expt{x}, catInfo(x), Tscan{x});
    [vesselROI{x}, NvesselROI{x}, tifStackMax{x}] = AnalyzeVasculature(expt{x}, projParam);

    % figure of vessel diameter changes over time
    for Z = find(NvesselROI{x})
        tempDiam = [vesselROI{x}{Z}.diameter];
        colorMat = distinguishable_colors(NvesselROI{x}(Z));
        close all;
        figure;
        subplot(1,2,1);
        imshow(tifStackMax{x}{Z}, []);
        for roi = 1:NvesselROI{x}(Z), drawpolygon('Position', vesselROI{x}{Z}(roi).boxPosition.xy, 'color',colorMat(roi,:));  end
        title( sprintf('%s Z projection %i', expt{x}.name, Z), 'Interpreter','none' )

        subplot(1,2,2)
        plot(vertcat(Tscan{x}{:}), vertcat(tempDiam.um_gauss)') %plot(vesselROI{Z}.diameter.um, 'k')
        xlabel('Time (sec)'); ylabel('Diameter (\mum)');
        axis tight; axis square;
        pause;
    end
    vesselROIpool = [vesselROI{x}{:}];
    diamPool = [vesselROIpool.diameter];
    allDiam = cat(1, diamPool.um_gauss)';
    allDiamZ = zscore(allDiam, [], 1);
    imagesc(allDiamZ')
end