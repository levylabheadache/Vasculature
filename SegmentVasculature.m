function [vesselROI, NvesselROI, tifStackMax] = SegmentVasculature(expt, projParam, varargin)

IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired(IP, 'projParam', @isstruct)
addParameter( IP, 'overwrite', false, @islogical )
addParameter( IP, 'Zredo', [], @isnumeric ) %planes to be resegmented
addParameter( IP, 'review', false, @islogical )
parse( IP, expt, projParam, varargin{:} ); % mouse, exptDate,
overwrite = IP.Results.overwrite;
Zredo = IP.Results.Zredo;
review = IP.Results.review;

% Determine which channel to segment
if ~isfield(expt, 'vascChan')
    warning('vascChan not specified, set to red by default')
    vascChan = 'red'; 
else
    vascChan = expt.vascChan; 
end
vascChanInd = find(contains({'red','green'}, vascChan));

% Load segmentation, if it already exists
vessPath = sprintf('%s%s_vasc.mat', expt.dir, expt.name);
if exist(vessPath, 'file') && ~overwrite
    fprintf('\nLoading %s', vessPath);
    load(vessPath, 'vesselROI','NvesselROI', 'tifStackMax');
else
    tifStackMax = cell(1,projParam.Nz); 
    vesselROI = cell(1,projParam.Nz); 
    NvesselROI = zeros(1,projParam.Nz);
end

% Draw vessel ROIs and estimate diameters
Zseg = [];
if ~isempty(Zredo)
    Zseg = Zredo;
elseif overwrite || ~exist(vessPath, 'file') %projParam.z
    Zseg = 1:projParam.Nz;
end
if expt.Nplane > 1
    for Z = Zseg
        if expt.Nruns == 1
            catPath = projParam.path.run.reg.z{1,vascChanInd,Z}; % projParam.path.run.reg.z(runs,:,Z)
        else
            catPath = projParam.path.cat.reg.z{Z,vascChanInd}; %#ok<FNDSB>
        end
        tifStack = loadtiff(catPath);
        tifStackMax{Z} = max(tifStack, [], 3);
        [vesselROI{Z}, NvesselROI(Z)] = MakeVesselROI( tifStackMax{Z} );

        for roi = 1:NvesselROI(Z)  
            vesselROI{Z}(roi).z = projParam.z{Z};  
        end
     
    end
else
    if Zseg == 1
        if expt.Nruns == 1
            catPath = projParam.path.run.raw.z{1,vascChanInd,1}; % projParam.path.run.reg.z(runs,:,Z)
        else
            catPath = projParam.path.cat.reg.z{1,vascChanInd,1}; %#ok<FNDSB>
        end
        disp(catPath);
        tifStack = loadtiff(catPath);
        tifStackMax{1} = max(tifStack, [], 3);
        [vesselROI{1}, NvesselROI(1)] = MakeVesselROI( tifStackMax{1} );

        for roi = 1:NvesselROI  
            vesselROI{1}(roi).z = projParam.z{1};  
        end
        
        vesselROI{1} = GetVesselProfile(vesselROI{1}, tifStackMax{1}, 'max');
        vesselROI{1} = GetVesselProfile(vesselROI{1}, tifStack);

        % Save the results
        fprintf('\nSaving %s', vessPath);
        save(vessPath, 'vesselROI','NvesselROI', 'projParam','tifStackMax');
    end
end

if review
    for Z = 1:projParam.Nz
        colorMat = distinguishable_colors(NvesselROI(Z));
        close all;
        figure;
        imshow(tifStackMax{Z}, []);
        impixelinfo;
        for roi = 1:NvesselROI(Z), drawpolygon('Position', vesselROI{Z}(roi).boxPosition.xy, 'color',colorMat(roi,:));  end
        title( sprintf('%s Z projection %i', expt.name, Z), 'Interpreter','none' )
        pause;
    end
end
end


%savefig = {}
