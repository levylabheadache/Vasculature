function RegisterVascData(expt, runInfo, varargin) % [outputArg1,outputArg2] = 
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'runInfo', @isstruct )
addParameter( IP, 'margin', 1, @isnumeric ) % grab this many planes above/below the central superficial/penetrating planes
addParameter( IP, 'edge', [80,80,40,20], @isnumeric ) % [60,60,40,40]
addParameter( IP, 'minEdge', [60,60,20,20], @isnumeric )
addParameter( IP, 'chunk', 20, @isnumeric ) %don't go over 20
addParameter( IP, 'refChan', 1, @isnumeric ) % for scanbox, 1 = green, 2 = red. -1 = both
addParameter( IP, 'refScans', 1, @isnumeric )
addParameter( IP, 'minInt', 1500, @isnumeric )
addParameter( IP, 'fix', false, @islogical )
addParameter( IP, 'overwrite', false, @islogical )
parse( IP, expt, runInfo, varargin{:} ); % mouse, exptDate,
margin = IP.Results.margin; 
fixSbx = IP.Results.fix;
initEdge = IP.Results.edge;
minEdge = IP.Results.minEdge;
minInt = IP.Results.minInt;
refChan = IP.Results.refChan;  %1; %1 = red, 2 = green
refScans = IP.Results.refScans; 
overwrite = IP.Results.overwrite;
fSep = '\';
[optimizer, metric] = imregconfig('monomodal');
Nrow = runInfo(1).sz(1);
Ncol = runInfo(1).sz(2);
Nplane = runInfo(1).Nplane;
totScan = sum([runInfo.Nscan]);

zSuper = expt.superficial-margin:expt.superficial+margin;
zSuper(zSuper < 1 | zSuper > Nplane) = [];
Nsuper = numel(zSuper);
zPen = expt.penetrating-margin:expt.penetrating+margin;
zPen(zPen < 1 | zPen > Nplane) = [];
Npen = numel(zPen);
planeData = cell(expt.Nruns, Nplane);
for runs = 1:expt.Nruns
    ZtifDir = strcat(runInfo(runs).dir, 'Ztifs', fSep); mkdir( ZtifDir );
    sbxInputPath = runInfo(runs).path; %sbxPath{runs};
    [fDir, fName] = fileparts(sbxInputPath);
    pathTemplate = strcat( fDir, fSep, fName );
    rawProjPath = strcat(pathTemplate,'_rawProj.tif');
    
    % Fix the sequencing issue if needed 
    if fixSbx 
        sbxFixPath = [pathTemplate, '.sbxfix'];
        if (~exist(sbxFixPath,'file') || overwrite)
            fprintf('\n   Correcting sbx z order... ');
            runInfo(runs) = FixSBX(sbxInputPath, runInfo(runs));
        end
        sbxInputPath = sbxFixPath;
    end
    
    % Write tifs from the raw data
    if (~exist(rawProjPath,'file') || overwrite)
        fprintf('\n   Writing raw projection stack... ');
        if runInfo(runs).Nplane > 1 
            WriteSbxProjection(sbxInputPath, runInfo(runs), rawProjPath, 'verbose',true, 'chan',refChan, 'dir',ZtifDir, 'name',runInfo(runs).exptName); % , 'binT',binT  [ZtifDir,'Raw\']
        else
            error('RegisterVascData is for volume data only');
        end
    end
    
    % Load the imaging data from the relevant planes
    for z = unique([zSuper, zPen])
        planeData{runs, z} = ...
            WriteSbxPlaneTif(sbxInputPath, runInfo(runs), z, 'chan',refChan, 'dir',ZtifDir, 'name',runInfo(runs).exptName);
    end
end

% SUPERFICIAL DATA
superDataRaw = zeros(Nrow, Ncol, totScan, Nsuper, 'uint16'); 
for Z = 1:Nsuper
    superDataRaw(:,:,:,Z) = cat(3, planeData{:,zSuper(Z)});
end
superProj = squeeze(uint16(mean(superDataRaw, 3)));
saveastiff(superProj, [expt.dir, expt.name, '_superStackUnaligned.tif']); % 'D:\2photon\CGRP02\201110_CGRP02\superProj_FOV2.tif'
%superMeanRaw = uint16( mean(superDataRaw, 4) ); 
%saveastiff(superMeanRaw, 'D:\2photon\CGRP02\201110_CGRP02\superMeanRaw.tif'); 

% Align other planes to the central one
Zcent = ceil(median(1:Nsuper));
superDataAlign = zeros(Nrow, Ncol, totScan, Nsuper, 'uint16'); 
for Z = flip(1:Nsuper)
    superAlignTform(Z) = imregtform(superProj(initEdge(3)+1:end-initEdge(4),initEdge(1)+1:end-initEdge(2),Z), superProj(initEdge(3)+1:end-initEdge(4),initEdge(1)+1:end-initEdge(2),Zcent), 'translation', optimizer, metric); % get shift relative to central plane
    % Apply shift to each frame
    if Z ~= Zcent
        superDataAlign(:,:,:,Z) = imwarp(superDataRaw(:,:,:,Z), superAlignTform(Z), 'OutputView',imref2d([Nrow,Ncol])); 
    else
        superDataAlign(:,:,:,Z) = superDataRaw(:,:,:,Z);
    end
end
% Determine edges of final aligned data
superProjAlign = uint16( squeeze(mean(superDataAlign, 3)) );
saveastiff(superProjAlign, [expt.dir, expt.name, '_superStackAligned.tif']);
superEdge = GetEdges(superProjAlign(:,:,1), 'minInt',minInt, 'show',false);
superEdge = max([superEdge; minEdge], [], 1);
superMeanAlign = uint16( mean(superDataAlign(superEdge(3)+1:end-superEdge(4),superEdge(1)+1:end-superEdge(2),:,:), 4) );
saveastiff(superMeanAlign, [expt.dir, expt.name, '_superMovie_unreg.tif']);

% Register each frame of each aligned submovie
superRef = uint16(median(superMeanAlign(:,:,refScans), 3));
saveastiff(superRef, [expt.dir, expt.name, '_superRef.tif']); % _median
tic
superMovieReg = zeros(size(superMeanAlign), 'uint16');
parfor s = 1:totScan %101:200 %flip()
    superMovieReg(:,:,s) = imregister( superMeanAlign(:,:,s), superRef, 'translation', optimizer, metric);
end
toc
saveastiff(superMovieReg, [expt.dir, expt.name, '_superMovie_reg.tif']);
toc

% PENETRATING DATA
penDataRaw = zeros(Nrow, Ncol, totScan, Npen, 'uint16'); 
for Z = 1:Npen
    penDataRaw(:,:,:,Z) = cat(3, planeData{:,zPen(Z)});
end
penProj = squeeze(uint16(mean(penDataRaw, 3)));
saveastiff(penProj, [expt.dir, expt.name, '_penStackUnaligned.tif']); 

% Align other planes to the central one
Zcent = ceil(median(1:Npen));
penDataAlign = zeros(Nrow, Ncol, totScan, Npen, 'uint16'); 
for Z = flip(1:Npen)
    penAlignTform(Z) = imregtform(penProj(initEdge(3)+1:end-initEdge(4),initEdge(1)+1:end-initEdge(2),Z), penProj(initEdge(3)+1:end-initEdge(4),initEdge(1)+1:end-initEdge(2),Zcent), 'translation', optimizer, metric); % get shift relative to central plane
    % Apply shift to each frame
    if Z ~= Zcent
        penDataAlign(:,:,:,Z) = imwarp(penDataRaw(:,:,:,Z), penAlignTform(Z), 'OutputView',imref2d([Nrow,Ncol])); 
    else
        penDataAlign(:,:,:,Z) = penDataRaw(:,:,:,Z);
    end
end
% Determine edges of final aligned data
penProjAlign = uint16( squeeze(mean(penDataAlign, 3)) );
saveastiff(penProjAlign, [expt.dir, expt.name, '_penStackAligned.tif']);
penEdge = GetEdges(penProjAlign(:,:,1), 'minInt',minInt, 'show',false);
penEdge = max([penEdge; minEdge], [], 1);
penMeanAlign = uint16( mean(penDataAlign(penEdge(3)+1:end-penEdge(4),penEdge(1)+1:end-penEdge(2),:,:), 4) );
saveastiff(penMeanAlign, [expt.dir, expt.name, '_penMovie_unreg.tif']);

% Register each frame of each aligned submovie
penRef = uint16(median(penMeanAlign(:,:,refScans), 3));
saveastiff(penRef, [expt.dir, expt.name, '_penRef.tif']); % _median
tic
penMovieReg = zeros(size(penMeanAlign), 'uint16');
parfor s = 1:totScan 
    penMovieReg(:,:,s) = imregister( penMeanAlign(:,:,s), penRef, 'translation', optimizer, metric);
end
toc
saveastiff(penMovieReg, [expt.dir, expt.name, '_penMovie_reg.tif']);
toc
end