%% Use GLM to assess contribution of different variables

% Clear any previous variables in the Workspace and Command Window to start fresh
clear; clc; close all;

% TODO -- Set the directory of where animal folders are located
dataDir =  'D:\2photon\Simone\Simone_Macrophages\'; % 'D:\2photon\Simone\Simone_Vasculature\'; 'D:\2photon\Simone\Simone_Macrophages\'; %  

% PARSE DATA TABLE 

% TODO --- Set excel sheet
dataSet = 'Vasculature'; %'Macrophage'; 'AffCSD'; 'Pollen'; 'Vasculature'; %  'Astrocyte'; %  'Anatomy'; %  'Neutrophil_Simone'; % 'Afferents'
[regParam, projParam] = DefaultProcessingParams(dataSet); % get default parameters for processing various types of data

regParam.method = 'translation'; %rigid 
regParam.name = 'translation'; %rigid  

% TODO --- Set data spreadsheet directory
dataTablePath = 'R:\Levy Lab\2photon\ImagingDatasets_Simone.xlsx'; % 'R:\Levy Lab\2photon\ImagingDatasetsSimone2.xlsx';
dataTable = readcell(dataTablePath, 'sheet',dataSet);  % 'NGC', ''
colNames = dataTable(1,:); dataTable(1,:) = [];
dataCol = struct('mouse',find(contains(colNames, 'Mouse')), 'date',find(contains(colNames, 'Date')), 'FOV',find(contains(colNames, 'FOV')), 'vascChan',find(contains(colNames, 'VascChan')),...
    'volume',find(contains(colNames, 'Volume')), 'run',find(contains(colNames, 'Runs')), 'Ztop',find(contains(colNames, 'Zbot')), 'Zbot',find(contains(colNames, 'Ztop')), 'csd',find(contains(colNames, 'CSD')), ...
    'ref',find(contains(colNames, 'Ref')), 'edges',find(contains(colNames, 'Edge')), 'Zproj',find(contains(colNames, 'Zproj')), 'done',find(contains(colNames, 'Done')));
Nexpt = size(dataTable, 1);
dataTable(:,dataCol.date) = cellfun(@num2str, dataTable(:,dataCol.date), 'UniformOutput',false);

% Initialize variables
locoDiamDeform_pred = cell(1,Nexpt); 
locoDiamDeform_resp = cell(1,Nexpt); 
locoDiamDeform_opts = cell(1,Nexpt); 
locoDiamDeform_result = cell(1,Nexpt); 
locoDiamDeform_summary = cell(1,Nexpt);
GLMname = 'locoDeformDiam';
%GLMrate = 15.49/30;

for x = xPresent % x3D % 
    % GLMparallel options
    locoDiamDeform_opts{x}.name = sprintf('%s_%s', expt{x}.name, GLMname); %strcat(expt{x}.name, , '_preCSDglm');
    locoDiamDeform_opts{x}.rShow = NaN;
    locoDiamDeform_opts{x}.figDir = ''; % figDir;
    locoDiamDeform_opts{x}.alpha = 0.01;  % The regularization parameter, default is 0.01
    locoDiamDeform_opts{x}.standardize = true; 
    locoDiamDeform_opts{x}.trainFrac = 0.75; % 1; %
    locoDiamDeform_opts{x}.Ncycle = 20;
    locoDiamDeform_opts{x}.distribution = 'gaussian'; % 'poisson'; %  
    locoDiamDeform_opts{x}.CVfold = 10;
    locoDiamDeform_opts{x}.nlamda = 1000;
    locoDiamDeform_opts{x}.maxit = 5*10^5;
    locoDiamDeform_opts{x}.minDev = 0.05; 
    locoDiamDeform_opts{x}.minDevFrac = 0.1;
    locoDiamDeform_opts{x}.maxP = 0.05;
    locoDiamDeform_opts{x}.Nshuff = 0;  
    locoDiamDeform_opts{x}.minShuff = 15; 
    locoDiamDeform_opts{x}.window = [-10,10]; % [0,0]; % [-0.5, 0.5]; % 
    locoDiamDeform_opts{x}.lopo = true; %false; %
    locoDiamDeform_opts{x}.frameRate = expt{x}.scanRate;  % GLMrate; %
    locoDiamDeform_opts{x}.binSize = 1; %expt{x}.scanRate/GLMrate;
    locoDiamDeform_opts{x}.minShuffFrame = round( locoDiamDeform_opts{x}.frameRate*locoDiamDeform_opts{x}.minShuff );
    windowFrame = round(locoDiamDeform_opts{x}.window*locoDiamDeform_opts{x}.frameRate);
    locoDiamDeform_opts{x}.shiftFrame = windowFrame(1):windowFrame(2);
    locoDiamDeform_opts{x}.maxShift = max( abs(windowFrame) );
    locoDiamDeform_opts{x}.Nshift = numel( locoDiamDeform_opts{x}.shiftFrame );  %Nshift = preCSDOpts(x).Nshift;
    locoDiamDeform_opts{x}.lags = locoDiamDeform_opts{x}.shiftFrame/locoDiamDeform_opts{x}.frameRate;
    locoDiamDeform_opts{x}.xVar = 'Time';

    % PREDICTORS
    % locomotion predictors
    tempVelocityCat = BinDownMean( vertcat(loco{x}(expt{x}.preRuns).Vdown), locoDiamDeform_opts{x}.binSize );
    tempAccelCat = BinDownMean( abs(vertcat(loco{x}(expt{x}.preRuns).Adown)), locoDiamDeform_opts{x}.binSize ); 
    tempStateCat = BinDownMean( vertcat(loco{x}(expt{x}.preRuns).stateDown), locoDiamDeform_opts{x}.binSize );

    % deformation predictors
    zUse = [projParam.z{:}];
    transMag = vertcat(deform{x}(expt{x}.preRuns).transMag);
    transPred = prctile(transMag(:,zUse),80,2);
    transPred(normalize(transPred) > 5) = NaN;
    transPred = BinDownMean( transPred, locoDiamDeform_opts{x}.binSize );
    % Scaling
    scaleMag = vertcat(deform{x}(expt{x}.preRuns).scaleMag);
    scalePred = prctile(scaleMag(:,zUse),80,2);
    scalePred(normalize(scalePred) > 5) = NaN;
    scalePred = BinDownMean( scalePred, locoDiamDeform_opts{x}.binSize );

    shearMag = vertcat(deform{x}(expt{x}.preRuns).shearMag);
    shearPred = prctile(shearMag(:,zUse),80,2);
    shearPred(normalize(shearPred) > 5) = NaN;
    shearPred = BinDownMean( shearPred, locoDiamDeform_opts{x}.binSize );
    
    transSpdMag = vertcat(deform{x}(expt{x}.preRuns).DtransMag);
    transSpdPred = prctile(transSpdMag(:,zUse),80,2);
    transSpdPred(normalize(transSpdPred) > 5) = NaN;
    transSpdPred = BinDownMean( transSpdPred, locoDiamDeform_opts{x}.binSize );
    
    stretchMag = vertcat(deform{x}(expt{x}.preRuns).stretchMag);
    stretchPred = prctile(stretchMag(:,zUse),80,2);
    stretchPred(normalize(stretchPred) > 5) = NaN;
    stretchPred = BinDownMean( stretchPred, locoDiamDeform_opts{x}.binSize );
    
    shearRateMag = vertcat(deform{x}(expt{x}.preRuns).DshearMag);
    shearRatePred = prctile(shearRateMag(:,zUse),80,2);
    shearRatePred(normalize(shearRatePred) > 5) = NaN;
    shearRatePred = BinDownMean( shearRatePred, locoDiamDeform_opts{x}.binSize );
    
    tempShift = vertcat(deform{x}(expt{x}.preRuns).shiftZ);
    shiftPred = prctile(tempShift(:,zUse),80,2);
    shiftPred(abs(normalize(shiftPred)) > 5) = NaN;
    shiftPred = BinDownMean(shiftPred, locoDiamDeform_opts{x}.binSize);
    
    locoDiamDeform_pred{x} = struct('data',[], 'name',[], 'N',NaN, 'TB',[], 'lopo',[], 'fam',[]); 
    locoDiamDeform_pred{x}.data = [transPred, scalePred, shearPred, transSpdPred, stretchPred, shearRatePred, shiftPred, tempVelocityCat, tempAccelCat, tempStateCat]; % 
    locoDiamDeform_pred{x}.name = {'|Translation|', '|Scale|', '|Shear|', 'TransSpd', 'Stretch', 'Shear Rate','Z Shift', 'Velocity', '|Accel|', 'State'};
    locoDiamDeform_pred{x}.N = size(locoDiamDeform_pred{x}.data,2);
    for p = flip(1:locoDiamDeform_pred{x}.N), locoDiamDeform_pred{x}.lopo.name{p} = ['No ',locoDiamDeform_pred{x}.name{p}]; end
    
    % Set up leave-one-family-out
    locoDiamDeform_pred{x}.fam.col = {1:7, 8:10}; %{1:4, 5:7}; %{1:2, 3:4, 5:6, 7:8, 9:10, 11:12};  % {1:12};%{1, 2:3, 4:5, 6:7, 8, 9}; 
    locoDiamDeform_pred{x}.fam.N = numel(locoDiamDeform_pred{x}.fam.col); 
    locoDiamDeform_pred{x}.fam.name = {'Deform','Loco'}; 

    % Define response
    vesselROIpool = [vesselROI{x}{:}];
    diamPool = [vesselROIpool.diameter];
    allDiam = cat(1, diamPool.um_gauss)';
    allDiamZ = zscore(allDiam, [], 1);
    diamResp = BinDownMean( allDiam, locoDiamDeform_opts{x}.binSize ); % allDiamZ
    locoDiamDeform_resp{x}.data = diamResp; 
    locoDiamDeform_resp{x}.N = size(locoDiamDeform_resp{x}.data, 2); 
    locoDiamDeform_resp{x}.name = sprintfc('Diam %i', 1:locoDiamDeform_resp{x}.N);

    % Remove scans with missing data 
    nanFrame = find(any(isnan([locoDiamDeform_pred{x}.data, locoDiamDeform_resp{x}.data]),2)); % find( isnan(sum(pred(x).data,2)) ); 
    fprintf('\nRemoving %i NaN-containing frames', numel(nanFrame));
    locoDiamDeform_pred{x}.data(nanFrame,:) = []; locoDiamDeform_resp{x}.data(nanFrame,:) = [];

    % Run the GLM
    locoDiamDeform_opts{x}.load = true; % false; % 
    locoDiamDeform_opts{x}.saveRoot = expt{x}.dir; %''; %
    [locoDiamDeform_result{x}, locoDiamDeform_summary{x}, ~, locoDiamDeform_pred{x}, locoDiamDeform_resp{x}] = GLMparallel(locoDiamDeform_pred{x}, locoDiamDeform_resp{x}, locoDiamDeform_opts{x}); 
    %locoDiamDeform_summary{x} = SummarizeGLM(locoDiamDeform_result{x}, locoDiamDeform_pred{x}, locoDiamDeform_resp{x}, locoDiamDeform_opts{x});
end
%%
for x = xPresent
    locoDiamDeform_opts{x}.rShow = 1:sum(NvesselROI{x}); %1:locoDiamDeform_resp{x}.N; % 1:LocoDeform_resp{x}.N; %NaN;
    locoDiamDeform_opts{x}.xVar = 'Time';
    ViewGLM(locoDiamDeform_pred{x}, locoDiamDeform_resp{x}, locoDiamDeform_opts{x}, locoDiamDeform_result{x}, locoDiamDeform_summary{x}); %GLMresultFig = 
end

%%