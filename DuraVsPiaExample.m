%% NaVAi6-4 220321 FOV1 Run1
crop_rect = [100,20,400,400];
redPct = [5, 99.5];
greenPct = [5,99.7];

z_dura = 5:8;
z_pia = 12:15;
% Load projection stacks
stack_green = loadtiff('D:\2photon\NaVAi6-4\220321_FOV1\NaVAi6-4_220321_FOV1_run1\NaVAi6-4_220321_001_dft_meanProj_green.tif');
stack_red = loadtiff('D:\2photon\NaVAi6-4\220321_FOV1\NaVAi6-4_220321_FOV1_run1\NaVAi6-4_220321_001_dft_meanProj_red.tif');

% 
dura_green = mean(stack_green(:,:,z_dura), 3);
dura_red = mean(stack_red(:,:,z_dura), 3);
dura_RGB = zeros([size(dura_green),3 ], 'uint8');
dura_RGB(:,:,1) = uint8(rescale(dura_red, 0, 2^8-1, 'inputMin',prctile(dura_red(:), redPct(1)), 'inputMax',prctile(dura_red(:), redPct(2)))); % min(stackChan{chan}(:))
dura_RGB(:,:,2) = uint8(rescale(dura_green, 0, 2^8-1, 'inputMin',prctile(dura_green(:), greenPct(1)), 'inputMax',prctile(dura_green(:), greenPct(2))));
dura_RGB =  imcrop(dura_RGB, crop_rect);


pia_green = mean(stack_green(:,:,z_pia), 3);
pia_red = mean(stack_red(:,:,z_pia), 3);
pia_RGB = zeros([size(pia_green),3], 'uint8');
pia_RGB(:,:,1) = uint8(rescale(pia_red, 0, 2^8-1, 'inputMin',prctile(pia_red(:), redPct(1)), 'inputMax',prctile(pia_red(:), redPct(2)))); % min(stackChan{chan}(:))
pia_RGB(:,:,2) = uint8(rescale(pia_green, 0, 2^8-1, 'inputMin',prctile(pia_green(:), greenPct(1)), 'inputMax',prctile(pia_green(:), greenPct(2))));
pia_RGB =  imcrop(pia_RGB, crop_rect);

figDir = 'D:\MATLAB\Figures\VascDeformGrant\';
close all
DuralProj = figure;
%subplot(1,2,1)
imshow(dura_RGB);
MakeScaleBar( [20/expt{x}.umPerPixel,0], {get(gca,'Xlim'), get(gca,'Ylim')}, [0.05,0.95], [0,0], 'label',false, 'color','w', 'width',3 ); % 
figPath = sprintf('%sNaVAi6-4_220321_FOV1_run1_VascFibers_dura.pdf', figDir);
exportgraphics(DuralProj, figPath, 'Resolution',300); fprintf('\nSaved %s', figPath);

PialProj = figure;
%subplot(1,2,2)
imshow(pia_RGB);
figPath = sprintf('%sNaVAi6-4_220321_FOV1_run1_VascFibers_pia.pdf', figDir);
exportgraphics(PialProj, figPath, 'Resolution',300); fprintf('\nSaved %s', figPath);


%% CGRPAi6-1_220831_FOV1P/D
crop_rect = [250,50,400,400];
redPct = [5, 95];
greenPct = [5,99.3];

% Dural
dura_green = loadtiff('D:\2photon\CGRPAi6-1\220831_FOV1D\CGRPAi6-1_220831_FOV1D_reg_meanProj_green.tif');
dura_red = loadtiff('D:\2photon\CGRPAi6-1\220831_FOV1D\CGRPAi6-1_220831_FOV1D_reg_meanProj_red.tif');
dura_RGB = zeros([size(dura_green),3 ], 'uint8');
dura_RGB(:,:,1) = uint8(rescale(dura_red, 0, 2^8-1, 'inputMin',prctile(dura_red(:), redPct(1)), 'inputMax',prctile(dura_red(:), redPct(2)))); % min(stackChan{chan}(:))
dura_RGB(:,:,2) = uint8(rescale(dura_green, 0, 2^8-1, 'inputMin',prctile(dura_green(:), greenPct(1)), 'inputMax',prctile(dura_green(:), greenPct(2))));
dura_RGB =  imcrop(dura_RGB, crop_rect);

% Pial
pia_green = loadtiff('D:\2photon\CGRPAi6-1\220831_FOV1P\CGRPAi6-1_220831_FOV1P_reg_meanProj_green.tif');
pia_red = loadtiff('D:\2photon\CGRPAi6-1\220831_FOV1P\CGRPAi6-1_220831_FOV1P_reg_meanProj_red.tif');
pia_RGB = zeros([size(pia_green),3], 'uint8');
pia_RGB(:,:,1) = uint8(rescale(pia_red, 0, 2^8-1, 'inputMin',prctile(pia_red(:), redPct(1)), 'inputMax',prctile(pia_red(:), redPct(2)))); % min(stackChan{chan}(:))
pia_RGB(:,:,2) = uint8(rescale(pia_green, 0, 2^8-1, 'inputMin',prctile(pia_green(:), greenPct(1)), 'inputMax',prctile(pia_green(:), greenPct(2))));
pia_RGB =  imcrop(pia_RGB, crop_rect);

figDir = 'D:\MATLAB\Figures\VascDeformGrant\';
close all
DuralProj = figure;
%subplot(1,2,1)
imshow(dura_RGB);
MakeScaleBar( [50/expt{x}.umPerPixel,0], {get(gca,'Xlim'), get(gca,'Ylim')}, [0.05,0.90], [0,0], 'label',false, 'color','w', 'width',3 ); % 
figPath = sprintf('%sCGRPAi6-1_220831_FOV1D_duralVascFibers.pdf', figDir);
exportgraphics(DuralProj, figPath, 'Resolution',300); fprintf('\nSaved %s', figPath);

PialProj = figure;
%subplot(1,2,2)
imshow(pia_RGB);
%saveastiff(pia_RGB, )
figPath = sprintf('%sCGRPAi6-1_220831_FOV1P_pialVascFibers.pdf', figDir);
exportgraphics(PialProj, figPath, 'Resolution',300); fprintf('\nSaved %s', figPath);

%% CGRPAi6-1_220831_FOV2P/D
crop_rect = [60,10,450,450];
redPct = [5, 99.5];
greenPct = [5,99.7];

% Dural
dura_green = loadtiff('D:\2photon\CGRPAi6-1\220831_FOV2D\CGRPAi6-1_220831_FOV2D_run1\CGRPAi6-1_220831_006_meanProj_green.tif');
dura_red = loadtiff('D:\2photon\CGRPAi6-1\220831_FOV2D\CGRPAi6-1_220831_FOV2D_run1\CGRPAi6-1_220831_006_meanProj_red.tif');
dura_RGB = zeros([size(dura_green),3 ], 'uint8');
dura_RGB(:,:,1) = uint8(rescale(dura_red, 0, 2^8-1, 'inputMin',prctile(dura_red(:), redPct(1)), 'inputMax',prctile(dura_red(:), redPct(2)))); % min(stackChan{chan}(:))
dura_RGB(:,:,2) = uint8(rescale(dura_green, 0, 2^8-1, 'inputMin',prctile(dura_green(:), greenPct(1)), 'inputMax',prctile(dura_green(:), greenPct(2))));
dura_RGB =  imcrop(dura_RGB, crop_rect);

% Pial
pia_green = loadtiff('D:\2photon\CGRPAi6-1\220831_FOV2P\CGRPAi6-1_220831_FOV2P_run1\CGRPAi6-1_220831_007_meanProj_green.tif');
pia_red = loadtiff('D:\2photon\CGRPAi6-1\220831_FOV2P\CGRPAi6-1_220831_FOV2P_run1\CGRPAi6-1_220831_007_meanProj_red.tif');
pia_RGB = zeros([size(pia_green),3], 'uint8');
pia_RGB(:,:,1) = uint8(rescale(pia_red, 0, 2^8-1, 'inputMin',prctile(pia_red(:), redPct(1)), 'inputMax',prctile(pia_red(:), redPct(2)))); % min(stackChan{chan}(:))
pia_RGB(:,:,2) = uint8(rescale(pia_green, 0, 2^8-1, 'inputMin',prctile(pia_green(:), greenPct(1)), 'inputMax',prctile(pia_green(:), greenPct(2))));
pia_RGB =  imcrop(pia_RGB, crop_rect);

figDir = 'D:\MATLAB\Figures\VascDeformGrant\';
close all
DuralProj = figure;
%subplot(1,2,1)
imshow(dura_RGB);
MakeScaleBar( [50/expt{x}.umPerPixel,0], {get(gca,'Xlim'), get(gca,'Ylim')}, [0.78,0.95], [0,0], 'label',false, 'color','w', 'width',3 ); % 
figPath = sprintf('%sCGRPAi6-1_220831_FOV2D_duralVascFibers.pdf', figDir);
exportgraphics(DuralProj, figPath, 'Resolution',300); fprintf('\nSaved %s', figPath);

PialProj = figure;
%subplot(1,2,2)
imshow(pia_RGB);
%saveastiff(pia_RGB, )
figPath = sprintf('%sCGRPAi6-1_220831_FOV2P_pialVascFibers.pdf', figDir);
exportgraphics(PialProj, figPath, 'Resolution',300); fprintf('\nSaved %s', figPath);

