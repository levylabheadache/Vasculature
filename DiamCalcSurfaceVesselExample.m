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
magStart = strfind(textHold{20}, ': ');
MScanData.notes.header.magnification = textHold{20}(magStart + 2:end);
rotationStart = strfind(textHold{19}, ': ');
MScanData.notes.header.rotation = textHold{19}(rotationStart + 2:end);
frameRateStart = strfind(textHold{24}, ': ');
MScanData.notes.header.frameRate=(textHold{24}(frameRateStart + 2:end - 3));
MScanData.notes.frameRate = 1/str2num(MScanData.notes.header.frameRate); %#ok<*ST2NM>
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
xFactor = micronsPerPixel/(str2num(MScanData.notes.header.magnification(1:end - 1)));

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

%% Opens the tiff file and gets the  vessel projections from the defined polygons
function [MScanData] = GetDiameterFromMovie(MScanData, fileID)
MScanData.notes.firstFrame = imread(fileID, 'TIFF', 'Index', 1);
fftFirstFrame = fft2(double(MScanData.notes.firstFrame));
X = repmat(1:MScanData.notes.xSize, MScanData.notes.ySize, 1);
Y = repmat((1:MScanData.notes.ySize)', 1, MScanData.notes.xSize);
MScanData.notes.vesselROI.projectionAngle = atand(diff(MScanData.notes.vesselROI.vesselLine.position.xy(:, 1))/diff(MScanData.notes.vesselROI.vesselLine.position.xy(:, 2)));
atand(diff(MScanData.notes.vesselROI.vesselLine.position.xy(:, 1))/diff(MScanData.notes.vesselROI.vesselLine.position.xy(:, 2)));

for theFrame = MScanData.notes.startframe:MScanData.notes.endframe
    rawFrame = imread(fileID, 'TIFF', 'Index', theFrame);
    fftRawFrame = fft2(double(rawFrame));
    
    [MScanData.notes.pixelShift(:, theFrame), ~] = DftRegistration(fftFirstFrame, fftRawFrame, 1);
    
    inpolyFrame = inpolygon(X + MScanData.notes.pixelShift(3, theFrame), Y + MScanData.notes.pixelShift(4, theFrame), MScanData.notes.vesselROI.boxPosition.xy(:, 1), MScanData.notes.vesselROI.boxPosition.xy(:, 2));
    boundedrawFrame = rawFrame.*uint16(inpolyFrame);
    MScanData.notes.vesselROI.projection(theFrame, :) = radon(boundedrawFrame, MScanData.notes.vesselROI.projectionAngle);
end

end

%% Calculate diameter using FWHM and get the baseline diameter
function [MScanData] = FWHM_MovieProjection(MScanData, theFrames)

for f = min(theFrames):max(theFrames)
    % Add in a 5 pixel median filter
    MScanData.data.rawVesselDiameter(f) = CalcFWHM(medfilt1(MScanData.notes.vesselROI.projection(f, :), 5));
end

MScanData.data.tempVesselDiameter = MScanData.data.rawVesselDiameter*MScanData.notes.xFactor;
[holdHist, d] = hist(MScanData.data.tempVesselDiameter, 0:.25:100);
[~, maxD] = max(holdHist);
MScanData.notes.vesselROI.modalFixedDiameter = d(maxD);

end


%% Calc full-width at half-max
function width = CalcFWHM(data,smoothing,threshold)
data = double(data(:));     % make sure this is column, and cast to double

% smooth data, if appropriate
if nargin < 2
    % smoothing not passed in, set to default (none)
    smoothing = 1;
end

if smoothing > 1
    data = conv2(data,rectwin(smoothing) ./ smoothing,'valid');
end

% subtract out baseline
%data = data - min(data);   %jd - don't subract out min, in case threshold was set externally

if nargin < 3
    %threshold = max(data)/2;
    offset = min(data);                           % find the baseline
    threshold = max(data - offset) / 2 + offset;  % threshold is half max, taking offset into account
end

aboveI = find(data > threshold);    % all the indices where the data is above half max

if isempty(aboveI)
    % nothing was above threshold!
    width = 0;
    return
end

firstI = aboveI(1);                 % index of the first point above threshold
lastI = aboveI(end);                % index of the last point above threshold

if (firstI-1 < 1) || (lastI+1) > length(data)
    % interpolation would result in error, set width to zero and just return ...
    width = 0;
    return
end

% use linear interpolation to get a more accurate picture of where the max was
% find value difference between the point and the threshold value,
% and scale this by the difference between integer points ...
point1offset = (threshold-data(firstI-1)) / (data(firstI)-data(firstI-1));
point2offset = (threshold-data(lastI)) / (data(lastI+1)-data(lastI));

point1 = firstI-1 + point1offset;
point2 = lastI + point2offset;

width = point2-point1;

end

%% Remove motion artifacts
function newVesselDiameter = RemoveMotion(vesselDiameter, baseline, diffThresh, rawThresh)
indx1 = find(diff(vesselDiameter) > diffThresh);
indx2 = find(abs((vesselDiameter - baseline)/baseline) > (rawThresh));
indx = union(indx1 + 1, indx2);   % indx: points need to be interpolated
indx0 = 1:length(vesselDiameter);
indx0(indx) = [];   % indx0: good points
count = 1;

if isempty(indx) ~= 1
    if indx0(1) ~= 1
        indx0 = [1:indx0(1) - 1, indx0];
    end
    
    for a = 1:length(indx0) - 1
        step = indx0(a + 1) - indx0(a);
        if step == 1
            newVesselDiameter(count) = vesselDiameter(count); %#ok<*AGROW>
        end
        
        if (step ~= 1)
            newVesselDiameter(count) = vesselDiameter(count);
            newVesselDiameter(count + 1:count + step - 1) = (vesselDiameter(indx0(a + 1)) + vesselDiameter(indx0(a)))/2;
        end
        
        count = count + step;
    end
    
    newVesselDiameter(count) = vesselDiameter(indx0(end));
    if indx(end) == length(vesselDiameter)
        newVesselDiameter(indx0(end) + 1:length(vesselDiameter)) = vesselDiameter(indx0(end));
    end
else
    newVesselDiameter = vesselDiameter;
end

end

%% Registers image via cross-correlation
function [output, Greg] = DftRegistration(buf1ft, buf2ft, usfac)
%________________________________________________________________________________________________________________________
% Utilized in analysis by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Code unchanged with the exception of this title block for record keeping
%
%   Last Opened: February 23rd, 2019
%________________________________________________________________________________________________________________________
%
% function [output Greg] = dftregistration(buf1ft,buf2ft,usfac);
% Efficient subpixel image registration by crosscorrelation. This code
% gives the same precision as the FFT upsampled cross correlation in a
% small fraction of the computation time and with reduced memory
% requirements. It obtains an initial estimate of the crosscorrelation peak
% by an FFT and then refines the shift estimation by upsampling the DFT
% only in a small neighborhood of that estimate by means of a
% matrix-multiply DFT. With this procedure all the image points are used to
% compute the upsampled crosscorrelation.
% Manuel Guizar - Dec 13, 2007

% Portions of this code were taken from code written by Ann M. Kowalczyk
% and James R. Fienup.
% J.R. Fienup and A.M. Kowalczyk, "Phase retrieval for a complex-valued
% object by using a low-resolution image," J. Opt. Soc. Am. A 7, 450-458
% (1990).

% Citation for this algorithm:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup,
% "Efficient subpixel image registration algorithms," Opt. Lett. 33,
% 156-158 (2008).

% Inputs
% buf1ft    Fourier transform of reference image,
%           DC in (1,1)   [DO NOT FFTSHIFT]
% buf2ft    Fourier transform of image to register,
%           DC in (1,1) [DO NOT FFTSHIFT]
% usfac     Upsampling factor (integer). Images will be registered to
%           within 1/usfac of a pixel. For example usfac = 20 means the
%           images will be registered within 1/20 of a pixel. (default = 1)

% Outputs
% output =  [error,diffphase,net_row_shift,net_col_shift]
% error     Translation invariant normalized RMS error between f and g
% diffphase     Global phase difference between the two images (should be
%               zero if images are non-negative).
% net_row_shift net_col_shift   Pixel shifts between images
% Greg      (Optional) Fourier transform of registered version of buf2ft,
%           the global phase difference is compensated for.

% Default usfac to 1
if exist('usfac')~=1, usfac=1; end %#ok<*EXIST>

% Compute error for no pixel shift
if usfac == 0
    CCmax = sum(sum(buf1ft.*conj(buf2ft)));
    rfzero = sum(abs(buf1ft(:)).^2);
    rgzero = sum(abs(buf2ft(:)).^2);
    error = 1.0 - CCmax.*conj(CCmax)/(rgzero*rfzero);
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax));
    output=[error,diffphase];
    
    % Whole-pixel shift - Compute crosscorrelation by an IFFT and locate the
    % peak
elseif usfac == 1
    [m,n]=size(buf1ft);
    CC = ifft2(buf1ft.*conj(buf2ft));
    [max1,loc1] = max(CC);
    [~,loc2] = max(max1);
    rloc=loc1(loc2);
    cloc=loc2;
    CCmax=CC(rloc,cloc);
    rfzero = sum(abs(buf1ft(:)).^2)/(m*n);
    rgzero = sum(abs(buf2ft(:)).^2)/(m*n);
    error = 1.0 - CCmax.*conj(CCmax)/(rgzero(1,1)*rfzero(1,1));
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax));
    md2 = fix(m/2);
    nd2 = fix(n/2);
    if rloc > md2
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end
    
    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end
    output=[error,diffphase,row_shift,col_shift];
    
    % Partial-pixel shift
else
    
    % First upsample by a factor of 2 to obtain initial estimate
    % Embed Fourier data in a 2x larger array
    [m,n]=size(buf1ft);
    mlarge=m*2;
    nlarge=n*2;
    CC=zeros(mlarge,nlarge);
    CC(m+1-fix(m/2):m+1+fix((m-1)/2),n+1-fix(n/2):n+1+fix((n-1)/2)) = ...
        fftshift(buf1ft).*conj(fftshift(buf2ft));
    
    % Compute crosscorrelation and locate the peak
    CC = ifft2(ifftshift(CC)); % Calculate cross-correlation
    [max1,loc1] = max(CC);
    [~,loc2] = max(max1);
    rloc=loc1(loc2);cloc=loc2;
    CCmax=CC(rloc,cloc);
    
    % Obtain shift in original pixel grid from the position of the
    % crosscorrelation peak
    [m,n] = size(CC); md2 = fix(m/2); nd2 = fix(n/2);
    if rloc > md2
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end
    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end
    row_shift=row_shift/2;
    col_shift=col_shift/2;
    
    % If upsampling > 2, then refine estimate with matrix multiply DFT
    if usfac > 2
        %%% DFT computation %%%
        % Initial shift estimate in upsampled grid
        row_shift = round(row_shift*usfac)/usfac;
        col_shift = round(col_shift*usfac)/usfac;
        dftshift = fix(ceil(usfac*1.5)/2); %% Center of output array at dftshift+1
        % Matrix multiply DFT around the current shift estimate
        CC = conj(dftups(buf2ft.*conj(buf1ft),ceil(usfac*1.5),ceil(usfac*1.5),usfac,...
            dftshift-row_shift*usfac,dftshift-col_shift*usfac))/(md2*nd2*usfac^2);
        % Locate maximum and map back to original pixel grid
        [max1,loc1] = max(CC);
        [~,loc2] = max(max1);
        rloc = loc1(loc2); cloc = loc2;
        CCmax = CC(rloc,cloc);
        rg00 = dftups(buf1ft.*conj(buf1ft),1,1,usfac)/(md2*nd2*usfac^2);
        rf00 = dftups(buf2ft.*conj(buf2ft),1,1,usfac)/(md2*nd2*usfac^2);
        rloc = rloc - dftshift - 1;
        cloc = cloc - dftshift - 1;
        row_shift = row_shift + rloc/usfac;
        col_shift = col_shift + cloc/usfac;
        
        % If upsampling = 2, no additional pixel shift refinement
    else
        rg00 = sum(sum( buf1ft.*conj(buf1ft) ))/m/n;
        rf00 = sum(sum( buf2ft.*conj(buf2ft) ))/m/n;
    end
    error = 1.0 - CCmax.*conj(CCmax)/(rg00*rf00);
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax));
    % If its only one row or column the shift along that dimension has no
    % effect. We set to zero.
    if md2 == 1
        row_shift = 0;
    end
    if nd2 == 1
        col_shift = 0;
    end
    output=[error,diffphase,row_shift,col_shift];
end

% Compute registered version of buf2ft
if (nargout > 1)&&(usfac > 0)
    [nr,nc]=size(buf2ft);
    Nr = ifftshift((-fix(nr/2):ceil(nr/2)-1));
    Nc = ifftshift((-fix(nc/2):ceil(nc/2)-1));
    [Nc,Nr] = meshgrid(Nc,Nr);
    Greg = buf2ft.*exp(1i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
    Greg = Greg*exp(1i*diffphase);
elseif (nargout > 1)&&(usfac == 0)
    Greg = buf2ft*exp(1i*diffphase);
end
return

    function out=dftups(in,nor,noc,usfac,roff,coff)
        % function out=dftups(in,nor,noc,usfac,roff,coff);
        % Upsampled DFT by matrix multiplies, can compute an upsampled DFT in just
        % a small region.
        % usfac         Upsampling factor (default usfac = 1)
        % [nor,noc]     Number of pixels in the output upsampled DFT, in
        %               units of upsampled pixels (default = size(in))
        % roff, coff    Row and column offsets, allow to shift the output array to
        %               a region of interest on the DFT (default = 0)
        % Recieves DC in upper left corner, image center must be in (1,1)
        % Manuel Guizar - Dec 13, 2007
        % Modified from dftus, by J.R. Fienup 7/31/06
        
        % This code is intended to provide the same result as if the following
        % operations were performed
        %   - Embed the array "in" in an array that is usfac times larger in each
        %     dimension. ifftshift to bring the center of the image to (1,1).
        %   - Take the FFT of the larger array
        %   - Extract an [nor, noc] region of the result. Starting with the
        %     [roff+1 coff+1] element.
        
        % It achieves this result by computing the DFT in the output array without
        % the need to zeropad. Much faster and memory efficient than the
        % zero-padded FFT approach if [nor noc] are much smaller than [nr*usfac nc*usfac]
        
        [nr,nc]=size(in);
        % Set defaults
        if exist('roff')~=1, roff=0; end
        if exist('coff')~=1, coff=0; end
        if exist('usfac')~=1, usfac=1; end
        if exist('noc')~=1, noc=nc; end
        if exist('nor')~=1, nor=nr; end
        % Compute kernels and obtain DFT by matrix products
        kernc=exp((-1i*2*pi/(nc*usfac))*( ifftshift((0:nc-1)).' - floor(nc/2) )*( (0:noc-1) - coff ));
        kernr=exp((-i1*2*pi/(nr*usfac))*( (0:nor-1).' - roff )*( ifftshift((0:nr-1)) - floor(nr/2)  ));
        out=kernr*in*kernc;
        return
        
    end

end
