clear; clc; close all;
% Get experiment metadata
dataTablePath = 'C:\Users\ablaeser\Documents\MATLAB\Vasculature\vascData.xlsx'; % 'D:\MATLAB\LevyLab\Dura\3D\3DduraData.xlsx';
[expt, Nexpt] = GetVascExpt(dataTablePath); % Parse data index and get data


% Initialize variables
runInfo = cell(1,Nexpt); Tscan = cell(1,Nexpt); loco = cell(1,Nexpt); % deform = cell(1,Nexpt);

% Registration parameters
regParams.method = 'translation';
regParams.refChan = 1;  %1; %1 = red, 2 = green
regParams.refRun = 1;
regParams.refScan = []; %[];
regParams.chunkSize = 1;
regParams.histmatch = false; % true; % 
regParams.binXY = 2;
regParams.binT = 1;
regParams.prereg = false;
regParams.highpass = 0;
regParams.lowpass = 0;
regParams.medFilter = [0,0,0];
regParams.minInt = 1000;
regParams.edges = []; %[80,80,20,20]; 

xPresent = 3; % Determine which experiments to analyze
for x = xPresent 
    % Get basic run-level metadata and data, and maybe perform registration on individual runs
    for run = flip(expt(x).runs)
        runInfo{x}(run) = MakeInfoStruct( expt(x), run ); % expt(x).mouse, expt(x).date, r
        %loco{x}(r) = GetLocoData( DuraDataPath(expt(x).mouse, expt(x).date, r, 'quad'), runInfo{x}(r) ); % 
        %RegisterCat3D( runInfo{x}(run), 'preaff',true, 'fix',true, 'overwrite',false, 'writeZ',false, 'refChan',1, 'minInt',1000);  % regParams, 
    end  
    RegisterVascData(expt(x), runInfo{x}(expt(x).runs), 'refScans',2:runInfo{x}(expt(x).runs(1)).Nscan, 'fix',true);
    %[Tscan{x}, runInfo{x}] = GetTime(runInfo{x}(expt(x).runs));

    % Run registration on concatenated data
    % plot( vertcat( loco{x}.speedDown ) ) % Find periods of stillness
    [catInfo(x), ~] = ConcatenateRunInfo(expt(x), runInfo{x}(expt(x).runs), 'suffix','sbxcat', 'overwrite',true); % Get concatenated metadata % expt(x)
    ConcatenateExptRuns(expt(x), runInfo{x}(expt(x).runs), catInfo(x), 'overwrite',false); % , 'setEdge',[190, 100, 44, 122]

    % Final registration
    RegisterCat3D( catInfo(x), regParams, 'overwrite',false, 'writeZ',false, 'Zint',3, 'chunk',30, 'preaff',false, 'minInt',1000); % expt(x).mouse, expt(x).date, expt(x).runs , 'cat',false

    % Get deformation data
    %[~, deform{x}, ~] = GetDeformCat3D( catInfo(x), deformLim, 'show',false);  %  deformCat, affParams
end
