% This script compares the segmentation of FSL FAST (data-driven) 
% with the FreeSurfer approach (which relies on some human priors) to
% see if the FreeSurfer approach introduces any systematic bias.

%% Load data
load('segmentationData.mat');

%% Compare FreeSurfer segmentation with FSL segmentation
% Total brain volume
y = fast.totalBrainVol; % FSL FAST
x = aseg.BrainSegVol; % FreeSurfer
f = figure('Color', 'white', 'position', [200,200,290,290], 'Renderer', 'painters');
s = loglog(x, y, '.k', 'MarkerSize', 25);
xlabel('FreeSurfer volume (mm^3)'); % FS total brain volume (excl brainstem)
ylabel('FSL volume (mm^3)'); % FSL GM + WM volume
title('Total brain volume');
lims = [min([x;y]) .* 0.9, max([x;y]) .* 1.1];
xlim(lims);
ylim(lims);
hold on; l = plot(lims, lims);
l.LineStyle = '--';
l.Color = [0.8,0.8,0.8];
legend({'Data', 'Y = X'}, 'Location', 'SouthEast');
mdl = fitlm(log(x),log(y)); % stats
ci = coefCI(mdl);

% Gray matter volume
y = fast.grayMatterVol; % FSL FAST
x = aseg.TotalGrayVol; % FreeSurfer
f = figure('Color', 'white', 'position', [200,200,290,290], 'Renderer', 'painters');
s = loglog(x, y, '.k', 'MarkerSize', 25);
xlabel('FreeSurfer volume (mm^3)'); % FS gray matter volume (excl brainstem)
ylabel('FSL volume (mm^3)'); % FSL GM volume
title('Gray matter volume');
lims = [min([x;y]) .* 0.9, max([x;y]) .* 1.1];
xlim(lims);
ylim(lims);
hold on; l = plot(lims, lims);
l.LineStyle = '--';
l.Color = [0.8,0.8,0.8];
legend({'Data', 'Y = X'}, 'Location', 'SouthEast');
mdl = fitlm(log(x),log(y)); % stats
ci = coefCI(mdl);

% White matter volume
y = fast.whiteMatterVol; % FSL FAST
x = aseg.CerebralWhiteMatterVol + ...
    aseg.Right_Cerebellum_White_Matter + ...
    aseg.Left_Cerebellum_White_Matter; % FreeSurfer
f = figure('Color', 'white', 'position', [200,200,290,290], 'Renderer', 'painters');
s = loglog(x, y, '.k', 'MarkerSize', 25);
xlabel('FreeSurfer volume (mm^3)'); % FS white matter matter volume (excl brainstem)
ylabel('FSL volume (mm^3)'); % FSL WM volume
title('White matter volume');
lims = [min([x;y]) .* 0.9, max([x;y]) .* 1.1];
xlim(lims);
ylim(lims);
hold on; l = plot(lims, lims);
l.LineStyle = '--';
l.Color = [0.8,0.8,0.8];
legend({'Data', 'Y = X'}, 'Location', 'SouthEast');
mdl = fitlm(log(x),log(y)); % stats
ci = coefCI(mdl);

%% Compare some allometric scaling exponents between FSL and FreeSurfer
% White matter vs total gray matter

% FreeSurfer
x = aseg.TotalGrayVol;
y = aseg.CerebralWhiteMatterVol + ...
    aseg.Right_Cerebellum_White_Matter + ...
    aseg.Left_Cerebellum_White_Matter;
mdlFS = fitlm(log(x),log(y));
ciFS = coefCI(mdlFS);
fFS = figureLogLog(y, x, mdlFS.Coefficients.Estimate(2), mdlFS.Coefficients.Estimate(1), ...
    'White matter volume (mm^2)', 'Gray matter volume (mm^3)', ...
    'addLine', 1, 'confidenceBands', mdlFS);

% FSL
x = fast.grayMatterVol;
y = fast.whiteMatterVol;
mdlFSL = fitlm(log(x),log(y));
ciFSL = coefCI(mdlFSL);
fFSL = figureLogLog(y, x, mdlFSL.Coefficients.Estimate(2), mdlFSL.Coefficients.Estimate(1), ...
    'White matter volume (mm^2)', 'Gray matter volume (mm^3)', ...
    'addLine', 1, 'confidenceBands', mdlFSL);
