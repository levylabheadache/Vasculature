function vesselROI = GetVesselDiameter(expt, projParam, vesselROI, varargin)
IP = inputParser;
addRequired(IP, 'expt', @isstruct )
addRequired(IP, 'projParam', @isstruct)
addRequired(IP, 'vesselROI', @iscell)
addParameter( IP, 'smooth', 0, @isnumeric ) % specify in microns
addParameter( IP, 'lp', 0.25, @isnumeric ) % specify in Hz
addParameter( IP, 'overwrite', false, @islogical )
addParameter( IP, 'show', true, @islogical )
parse( IP, expt, projParam, vesselROI, varargin{:} ); 
smooth_um = IP.Results.smooth;
smooth_pix = round(smooth_um/projParam.umPerPixel_scale);
lpFreq = IP.Results.lp;
overwrite = IP.Results.overwrite;
show = IP.Results.show;
vessPath = sprintf('%s%s_vasc.mat', expt.dir, expt.name);
% Estimate diameter by curve fitting
for Z = 1:projParam.Nz
    if ~isempty(vesselROI{Z}) && (~isfield(vesselROI{Z}, 'diameter') || overwrite)
        vesselROI{Z} = EstimateVesselDiameter(vesselROI{Z}, projParam.umPerPixel_scale, projParam.rate_bin, 'smooth',smooth_pix, 'lp',lpFreq, 'show',true); % % round(4/expt.umPerPixel) 'motion',false,
        % Save the results
        fprintf('\nSaving %s', vessPath);
        save(vessPath, 'vesselROI', '-append');
    end
end
end