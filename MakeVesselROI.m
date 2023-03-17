function [vesselROI, NvesselROI] = MakeVesselROI( tifStackMax  )
% GUI for segmenting vessel ROIs, to measure their diameter
roiTemplate = struct( 'boxPosition',[], 'xSize',[], 'ySize',NaN, 'vesselLine',[], 'projectionAngle',[], 'projection',NaN, 'vesselType','' ); % 'modalFixedDiameter',[], 
vesselROI = repmat(roiTemplate, 0, 1);  % initialize  vessels(x).notes.vesselROI

[ySize, xSize] = size(tifStackMax);
FS = 12;
tifStackMaxLims = prctile(tifStackMax(:), [5,99.5]);
% Draw vessel ROI and axis line
close all;
figure('WindowState','maximized', 'color','k');
imshow(tifStackMax, tifStackMaxLims)
colormap('gray');
axis image

for roi = 1:50 % NthreshROI % while ishandle( proj.mean ) && c <= NthreshROI %
    title('Zoom in on a segment, then unpause', 'FontSize',FS, 'color','w'); pause;
    title(sprintf('roi = %d: Press any button to continue segmentation OR close the figure to stop',roi), 'FontSize',FS, 'Color','w');
    try
        waitforbuttonpress;
    catch
        roi = roi - 1; %#ok<FXSET>
        fprintf('\n  Figure closed. Found %d good ROI', roi );
        break;
    end

    % Box
    title('Draw a box around a small segment of the blood vessel with clear boundaries, then unpause')
    tempBox = drawpolygon; pause;
    vesselROI(roi).boxPosition.xy = tempBox.Position; 
    vesselROI(roi).xSize = xSize;
    vesselROI(roi).ySize = ySize;

    % Line
    title('Place the line along the axis of that small segment of the blood vessel, then unpause')
    diamAxis = drawline(); pause; 
    vesselROI(roi).vesselLine.position.xy = diamAxis.Position; 

    % Select vessel type
    flag = false;
    while ~flag
        vesselROI(roi).vesselType = input('Is the vessel an artery (A), vein (V), or dural (D) vessel?: ', 's'); disp(' '); % vessels(x).notes.vesselType
        answer = vesselROI(roi).vesselType; %vessels(x).notes.vesselType;
        if strcmpi(answer, 'A') || strcmpi(answer, 'D') || strcmpi(answer, 'V')
            flag = true;
        end
    end
    vesselROI(roi).projectionAngle = atand(diff(vesselROI(roi).vesselLine.position.xy(:, 1))/diff(vesselROI(roi).vesselLine.position.xy(:, 2)));
end
NvesselROI = numel(vesselROI);
end