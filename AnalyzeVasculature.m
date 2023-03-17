function [vesselROI, NvesselROI, tifStackMax] = SegmentVasculature(expt, projParam, overwrite)
if nargin < 3, overwrite = false; end
vessPath = sprintf('%s%s_vasc.mat', expt.dir, expt.name);
if ~isfield(expt, 'vascChan')
    warning('vascChan not specified, set to red by default')
    vascChan = 'red'; 
else
    vascChan = expt.vascChan; 
end
if ~exist(vessPath, 'file') || overwrite
    vascChanInd = find(contains({'red','green'}, vascChan));
    % Draw vessel ROIs and estimate diameters
    tifStackMax = cell(1,projParam.Nz); vesselROI = cell(1,projParam.Nz); NvesselROI = zeros(1,projParam.Nz);
    for Z = 1:projParam.Nz
        if expt.Nruns == 1
            catPath = projParam.path.run.reg.z{1,vascChanInd,Z}; % projParam.path.run.reg.z(runs,:,Z)
        else
            catPath = projParam.path.cat.reg.z{Z,vascChanInd}; %#ok<FNDSB> 
        end
        tifStack = loadtiff(catPath);
        tifStackMax{Z} = max(tifStack, [], 3);
        [vesselROI{Z}, NvesselROI(Z)] = MakeVesselROI( tifStackMax{Z} );
        for roi = 1:NvesselROI(Z),  vesselROI{Z}(roi).z = projParam.z{Z};  end
        vesselROI{Z} = GetVesselProfile(vesselROI{Z}, tifStack);
        vesselROI{Z} = GetVesselProfile(vesselROI{Z}, tifStackMax{Z}, 'max');
        %find(all(~vesselROI(roi).projection, 2))
        %vesselROI{Z} = EstimateVesselDiameter(vesselROI{Z}, projParam.umPerPixel_scale, projParam.rate_bin, 'motion',false, 'smooth',round(4/expt.umPerPixel), 'lp',0.25, 'show',true); %
        % Save the results
        fprintf('\nSaving %s', vessPath);
        save(vessPath, 'vesselROI','NvesselROI', 'projParam','tifStackMax');
    end

    % Estimate diameter by curve fitting
    for Z = 1:projParam.Nz
        vesselROI{Z} = EstimateVesselDiameter(vesselROI{Z}, projParam.umPerPixel_scale, projParam.rate_bin, 'motion',false, 'smooth',round(4/expt.umPerPixel), 'lp',0.25, 'show',true); %
        % Save the results
        fprintf('\nSaving %s', vessPath);
        save(vessPath, 'vesselROI','NvesselROI', 'projParam','tifStackMax');
    end
else
    fprintf('\nLoading %s', vessPath);
    load(vessPath, 'vesselROI','NvesselROI', 'tifStackMax');
end