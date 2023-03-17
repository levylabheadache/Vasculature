function [MScanData] = DiamCalcSurfaceVesselExample
%________________________________________________________________________________________________________________________
% Written/Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyzes the change in vessel diameter over time for a surface vessel (artery, vein, etc)
%________________________________________________________________________________________________________________________
%
%   Inputs: Click play - you will be prompted to select a single .tif stack. The stack used in this example is avaible
%           for download at https://drive.google.com/drive/folders/1YKzuRAuJtlcBx3zLxZods0PeUsiUrTlm?usp=sharing
%
%           The objective size for the provided stack is 16X, but is left as an input as this function is intended to be
%           shared and edited.
%
%           The vessel type for the provided stack is an artery (A).
%
%           You will then be prompted to draw two boxes, the first is an ROI around the diameter of the vessel. The second
%           is a line along the vessel center-axis. See README for more details.
%
%   Outputs: MScanData.mat struct, movie (5X speed), fig showing vessel diameter over time
%
%   Last Revised: May 9th, 2019
%________________________________________________________________________________________________________________________

clear
clc
%cd 'D:\2photon\CGRPAi6-1\220801_FOV1\Projections\Runs' % need to change diretory to targe folder first
disp('Select a single .tif stack from the current directory'); disp(' ')
tifStack = uigetfile('*.tif', 'Multiselect', 'off');
disp(['Loading: ' tifStack '...']); disp(' ')
movieInfo = imfinfo(tifStack);   % Pull information from graphics file

%take file info and extrac magnification and frame rate
MScanData.notes.header.fileName = movieInfo(1).Filename;
MScanData.notes.header.frameWidth = num2str(movieInfo(1).Width);
MScanData.notes.header.frameHeight = num2str(movieInfo(1).Height);
MScanData.notes.header.numberOfFrames = length(movieInfo);
MScanData.notes.xSize = str2double(MScanData.notes.header.frameWidth);
MScanData.notes.ySize = str2double(MScanData.notes.header.frameHeight);

% Read header and take further action based on header information
textHold = strread(movieInfo(1).ImageDescription, '%s', 'delimiter', '\n'); %#ok<*DSTRRD>
%magStart = strfind(textHold{20}, ': ');
%MScanData.notes.header.magnification = textHold{20}(magStart + 2:end);
MScanData.notes.header.magnification = 16*2.4;
%rotationStart = strfind(textHold{19}, ': ');
MScanData.notes.header.rotation = 0; % textHold{19}(rotationStart + 2:end);
%frameRateStart = strfind(textHold{24}, ': ');
%MScanData.notes.header.frameRate = 0.5; %(textHold{24}(frameRateStart + 2:end - 3));
MScanData.notes.frameRate = 0.5; % 1/str2num(MScanData.notes.header.frameRate); %#ok<*ST2NM>
MScanData.notes.startframe = 1;
MScanData.notes.endframe = MScanData.notes.header.numberOfFrames;

% Select objective size
flag = false;
while ~flag
    disp('Which objective was used during this stack? Objective IDs:')
    disp('      [1] Olympus UMPlanFL N 10X')
    disp('      [2] Nikon CFI75 LWD 16X')
    disp('      [3] Olympus UMPlanFL N 20X (0.50 NA)')
    disp('      [4] Olympus UMPlanFL N 20X (1.00 NA)')
    disp('      [5] Olympus UMPlanFL N 40X')
    answer = input('Enter objective ID number [1 2 3 4 5]: ', 's'); disp(' ')
    if strcmp(answer, '1') || strcmp(answer, '2') || strcmp(answer, '3') || strcmp(answer, '4') || strcmp(answer, '40X')
        flag = true;
    end
end

% Handle response
switch answer
    case '1'
        micronsPerPixel = 1.2953;
    case '2'
        micronsPerPixel = 0.825;
    case '3'
        micronsPerPixel = 0.5595;       
    case '4'
        micronsPerPixel = 0.64;    
    case '5'
        micronsPerPixel = 0.3619; 
end

% Select vessel type
flag = false;
while ~flag
    MScanData.notes.vesselType = input('Is the vessel an artery (A), vein (V), or dural (D) vessel?: ', 's'); disp(' ')
    answer = MScanData.notes.vesselType;
    if strcmp(answer, 'A') || strcmp(answer, 'D') || strcmp(answer, 'V')
        flag = true;
    end
end

MScanData.notes.micronsPerPixel = micronsPerPixel;
MScanData.notes.header.timePerLine = 1/(MScanData.notes.frameRate*str2num(MScanData.notes.header.frameHeight));
xFactor = micronsPerPixel/MScanData.notes.header.magnification; %(str2num(MScanData.notes.header.magnification(1:end - 1)));

% Draw vessel ROI and axis line
image = imread(tifStack, 'TIFF', 'Index', 1);
figure;
imagesc(double(image))
title('Two-photon vessel example')
colormap('gray');
axis image
xlabel('pixels')
ylabel('pixels')

yString = 'y';
theInput = 'n';
xSize = size(image, 2);
ySize = size(image, 1);
area = impoly(gca, [1 1; 1 20; 20 20; 20 1]); %#ok<*IMPOLY>

while strcmp(yString, theInput) ~= 1
    theInput = input('Is the diameter of the box ok? (y/n): ', 's'); disp(' ')
end

if strcmp(yString, theInput)
    get_API = iptgetapi(area);
    MScanData.notes.vesselROI.boxPosition.xy = get_API.getPosition();
    MScanData.notes.vesselROI.xSize = xSize;
    MScanData.notes.vesselROI.ySize = ySize;
    theInput = 'n';
end

diamAxis = imline(gca, round(xSize*[.25 .75]), round(ySize*[.25 .75])); %#ok<*IMLINE>
while strcmp(yString, theInput) ~= 1
    theInput = input('Is the line along the diameter axis ok? (y/n): ', 's'); disp(' ')
end

if strcmp(yString, theInput)
    get_API = iptgetapi(diamAxis);
    MScanData.notes.vesselROI.vesselLine.position.xy = get_API.getPosition();
end
MScanData.notes.xFactor = xFactor;

disp('Analyzing vessel projections from defined polygons...'); disp(' ');
[MScanData] = GetDiameterFromMovie(MScanData, tifStack);

% Calc FWHM
try
    [MScanData] = FWHM_MovieProjection(MScanData, [MScanData.notes.startframe MScanData.notes.endframe]);
catch
    disp([MScanData.notes.imageID ' FWHM calculation failed!']); disp(' ')
end

% Attempt to remove SMALL x-y motion artifacts, does not work well in vertical Z-plane
% 1 dural/vein, >40% changes spline, artery: >60% spline
% 2 dural/vein, >30% changes interpolate, artery: >50% interpolate
if strcmp(MScanData.notes.vesselType, 'D') || strcmp(MScanData.notes.vesselType, 'V')
    MScanData.data.vesselDiameter = RemoveMotion(MScanData.data.tempVesselDiameter, MScanData.notes.vesselROI.modalFixedDiameter, 2, 0.3);
else
    MScanData.data.vesselDiameter = RemoveMotion(MScanData.data.tempVesselDiameter, MScanData.notes.vesselROI.modalFixedDiameter, 2, 0.5);
end

% Create .avi movie of the vessel. File name is hardcoded, can be changed to input for filename.
disp('Generating movie (sped-up 5X)...'); disp(' ')
outputVideo = VideoWriter('VesselDiamCalcExample.avi');
outputVideo.FrameRate = MScanData.notes.frameRate*5;
open(outputVideo);
for a = 1:MScanData.notes.header.numberOfFrames
    tifImg = double(imread(tifStack, 'TIFF', 'Index', a));
    writeVideo(outputVideo, mat2gray(tifImg));
end
close(outputVideo)

% Lowpass filter for vessel diameter in figure
[B, A] = butter(4, 1/(MScanData.notes.frameRate/2), 'low');
filtVesselDiam = filtfilt(B, A, MScanData.data.vesselDiameter);

% figure of vessel diameter changes over time
figure;
subplot(2,1,1)
plot((1:MScanData.notes.header.numberOfFrames) / MScanData.notes.frameRate, ...
    MScanData.data.vesselDiameter, 'k')
title('Vessel diameter over time')
xlabel('Time (sec)')
ylabel('Diameter (\mum)')
axis tight

subplot(2,1,2)
plot((1:MScanData.notes.header.numberOfFrames) / MScanData.notes.frameRate, filtVesselDiam, 'k')
title('With 1 Hz lowpass filter')
xlabel('Time (sec)')
ylabel('Diameter (\mum)')
axis tight

implay('VesselDiamCalcExample.avi')   % show movie

disp('Diameter calculation for surface vessel example  - complete.'); disp(' ')

end










