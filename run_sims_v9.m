% Script to run various simulations, via the examples contained within the
% MEDTEG.m class. The outputs from the present script form the "results" of
% the associated publication:
%
%   XXXXXXXX
%
% Dependencies:
%   Requires MEDTEG.m
%   Requires the fig-matlab package (https://github.com/petejonze/fig-matlab)
%   May require various MathWorks packages (e.g., stats) -- untested
%
% History:
%   0.0.1	PRJ	06/09/2024 : initial piloting
%   0.0.2	PRJ	09/09/2024 : further development and modifications to MEDTEG.m and plotting
%   0.0.3	PRJ	12/09/2024 : first full working version (Sim1 & Sim2)
%   0.0.4	PRJ	13/09/2024 : starting to add Sim3 and Sim4
%   0.0.5	PRJ	15/09/2024 : working versions of Sim3 and Sim4
%   0.0.6	PRJ	17/09/2024 : corrected and simplified Sim1
%   0.0.7	PRJ	19/09/2024 : incremenetal tweaks and aesthetic improvements
%   0.0.8	PRJ	20/09/2024 : first full manuscript draft
%   0.0.9	PRJ	19/12/2024 : R1 (with simple max-gradient algorithm for comparison)


% init (clear/close everything)
clc
close all;
clear all %#ok


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% user params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For print-quality manuscript [~N hours for all 4 Sims]
% N_SIM = [1000 1000 1000 1000]; % Sims1-4
% N_IND_SIMS_TO_PLOT = 5;
% N_BOOTSTRP = 20000;

% For almost-final writeup [1.1 hours for all 4 Sims]  [comment out otherwise]
% N_SIM = [200 200 200 200]; % Sims1-4
% N_IND_SIMS_TO_PLOT = 5;
% N_BOOTSTRP = 2000;
% SAVE_PLOTS = true;

% For initial drafting [~25 mins for all 4 Sims]  [comment out otherwise]
N_SIM = [50 50 50 50]; % Sims1-4
N_IND_SIMS_TO_PLOT = 5;
N_BOOTSTRP = 1000;
SAVE_PLOTS = true;

% for rapid prototyping/debugging [~2.5 mins for all 4 Sims] [comment out otherwise]
% N_SIM = [5 5 5 5]; % Sims1-4
% N_IND_SIMS_TO_PLOT = 1;
% N_BOOTSTRP = 100;
% SAVE_PLOTS = false;


% ************* SHOULD BE NO NEED TO EDIT BELOW THIS POINT *************

% initialise output options
EXPORT_DIR      = [];
EXPORT_FORMAT   = 'pdf';
pkgVer          = 0.3;
[EXPORT_DIR, EXPORT_FORMAT] = fig_init(pkgVer,EXPORT_FORMAT,EXPORT_DIR);

% start timer (so we can see how long the script takes to run)
tic();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIM 1: Small scotoma -- example 1 in MEDTEG.runExamples()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n\n ***** SIM 1 ***** \n\n')

% =========================================================================
% 1.1 run simulations
% =========================================================================

    % params
    N_SIMS = N_SIM(1); 
    N_EXTRA_POINTS = 5; % to match in MEDTEG.m
    
    % init results data storage
    nExtraLocations = 0:5;
    nLocsInScotoma = nan(N_SIMS, N_EXTRA_POINTS+1); % incl. 0
    HillOfVision_volume_est = nan(N_SIMS, N_EXTRA_POINTS+1);
    HillOfVision_volume_true= nan(N_SIMS,1);
    sum_squared_residuals = nan(N_SIMS,N_EXTRA_POINTS+1);


    % create a "ground truth" to validate against
    id = [ % 24-2, right eye format (OD)
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN % +27
        NaN, NaN, NaN,   1,   2,   3,   4, NaN, NaN, NaN % +21
        NaN, NaN,   5,   6,   7,   8,   9,  10, NaN, NaN % +15
        NaN,  11,  12,  13,  14,  15,  16,  17,  18, NaN % +9
        19,   20,  21,  22,  23,  24,  25,  26,  27, NaN % +3
        28,   29,  30,  31,  32,  33,  34,  35,  36, NaN % -3
        NaN,  37,  38,  39,  40,  41,  42,  43,  44, NaN % -9
        NaN, NaN,  45,  46,  47,  48,  49,  50, NaN, NaN % -15
        NaN, NaN, NaN,  51,  52,  53,  54, NaN, NaN, NaN % -21
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN % -27
        ];  %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
    [x, y] = meshgrid(linspace(-27, 27, 10), linspace(27, -27, 10));
    x = x(~isnan(id));
    y = y(~isnan(id));

    % create a "ground truth" to validate against
    N = 100;
    [x1,y1] = meshgrid(linspace(min(x), max(x), N), linspace(min(y), max(y), N));

    % remove outside convex hull
    k = convhull([x y]);
    chverts = [x(k) y(k)];
    idx = ~inhull([x1(:) y1(:)], chverts, [], 0.1);
    x1(idx) = NaN;
    y1(idx) = NaN;

    % vectorised version, without NaNs
    x_hidensity = x1(~isnan(x1));
    y_hidensity = y1(~isnan(y1));

    % run
    for i = 1:N_SIMS
        close all; % prevent build up of windows  

        try
            [x_deg, y_deg, simObs_true_DLS, F_simObs_DLS, dat] = MEDTEG.runExample(1);
        catch ME
            warning(ME.message);
            continue
        end


        % measure volume (see volumetric.m)
        % NB: should technically convert to steradians, as areas "cannot be represented perfectly without distortion of any two-dimensional projection of a curved surface.... The graphied visual field... is an azimuthal map projection of the inner surface of the perimetry bowl. This means that distances and directions from the center of the graph are correctly represeted but area, circumferential distances, and directional relationship other than from the center are not correctly represented."
        % - Weleber & Tobler, 1986 (American Journal of Ophthalmology Volume 101, Issue 4, April 1986, Pages 461468
        % However, as shown in Table III of the above, the distortion is relatively minor within the central 30-degrees.
        if (exist('scatterquad2','file')~=2)
            error('For extended, volumetric analysis, add scatterquad2.m to your MATLAB path (available from Mathworks FileExchange');
        end
    
        % unpack dat matrix  [x y norm_DLS observer_true_DLS observer_est_DLS isNewlyAddedPoint]
        x = dat.x_deg;
        y = dat.y_deg;
        norm_DLS = dat.norm_DLS;
        observer_true_DLS = dat.observer_true_DLS;
        estimates_dB = dat.observer_est_DLS;
        error_dB = estimates_dB - observer_true_DLS;
        isNew = dat.isNewlyAddedPoint;
    
        % create a "ground truth" to validate against
        true_dB_hidensity = F_simObs_DLS(x_hidensity, y_hidensity);
        HillOfVision_volume_true(i) = scatterquad2(x_hidensity, y_hidensity, true_dB_hidensity);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute error (volumetric deviation from ground truth)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % compute HoV with 0...n extra points
        HillOfVision_volume_est(i,1) = scatterquad2(x(isNew==0), y(isNew==0), estimates_dB(isNew==0));
        nidx = find(isNew); % new point indices
        % defensive validation
        if length(nidx) ~= N_EXTRA_POINTS
            error('??????');
        end
        % compute
        for j = 1:N_EXTRA_POINTS
            idx = 1:nidx(j);
            HillOfVision_volume_est(i,j+1) = scatterquad2(x(idx), y(idx), estimates_dB(idx));
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute error 2 (sum sequared error)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % compute sum squared residuals with 0...n extra points
        F_DLS = scatteredInterpolant(x(isNew==0), y(isNew==0), estimates_dB(isNew==0));
        DLS_estimates = F_DLS(x_hidensity, y_hidensity);
        sum_squared_residuals(i,1) = sum((true_dB_hidensity - DLS_estimates).^2);

        % compute
        nidx = find(isNew); % new point indices
        for j = 1:length(nidx)
            idx = 1:nidx(j);
            F_DLS = scatteredInterpolant(x(idx), y(idx), estimates_dB(idx));
            DLS_estimates = F_DLS(x_hidensity, y_hidensity);
            sum_squared_residuals(i,j+1) = sum((true_dB_hidensity - DLS_estimates).^2);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute N test locations within scotoma(s)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        TD = observer_true_DLS - norm_DLS;
        inScotoma = TD < -8;
        inBlindSpot = x > 13 & x < 17 & y > -4 & y < 4; % crude but good enough
        idxOriginal = inScotoma & ~inBlindSpot & isNew==0;
        idxAll = inScotoma & ~inBlindSpot;
    
        % [x y TD TD<-10 inBlindSpot isNew]
        % fprintf('N points within regions of damage (within -6 dB or more loss region, excl. physiologic blindspot):\n')
        % fprintf('   original locations only: %i\n', sum(idxOriginal));
        % fprintf('   w/ additional locations: %i\n', sum(idxAll));
    
        nLocsInScotoma(i,1) = sum(idxOriginal);
        for j = 1:length(nidx)
            nLocsInScotoma(i,j+1) = nLocsInScotoma(i,j) + idxAll(nidx(j));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot first N individual examples
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % skip if plotted enough already
        if (i > N_IND_SIMS_TO_PLOT)
            continue; % skip
        end

        % open a new figure window
        figDim      = [15 20];
        nPlots      = [1 2];
        isTight     = true;
        isSquare    = true;
        styleFlag   = [];
        axesLims    = [];
        hFig        = fig_make(figDim, nPlots, isTight, isSquare, styleFlag, axesLims);

        % extract data
        idx = dat.isNewlyAddedPoint==1;
        existingTestLocations_XY = [dat.x_deg(~idx) dat.y_deg(~idx)];
        existingTestLocations_DLS_dB = dat.observer_est_DLS(~idx);
        newTestLocations_XY = [dat.x_deg(idx) dat.y_deg(idx)];
        newTestLocations_DLS_dB = dat.observer_est_DLS(idx);

        % Plot 1
        fig_subplot(1);
        hold on
        plot(existingTestLocations_XY(:,1), existingTestLocations_XY(:,2), 'ko', 'markerfacecolor',[.8 .8 .8])
        plot([15 15],[-3 3],'kx');
        for ii = 1:length(newTestLocations_DLS_dB)
            txt = sprintf('%i', ii);
            text(newTestLocations_XY(ii,1), newTestLocations_XY(ii,2), txt, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 7, 'Color','k', 'Interpreter','tex');
        end
        axis square
        plot(xlim(), [0 0], 'k--', [0 0], ylim(), 'k--'); % crosshair

        % title
        hTxt = textLoc('\bf Test Locations / Order', 'NorthWest', 'Interpreter','tex');
        fig_nudge(hTxt, -.075, .075);

        % format the axes
        xTick       = -20:10:20;
        xTickLabels = [];
        yTick       = xTick;
        yTickLabels = [];
        xAxisTitle  = [];
        yAxisTitle  = [];
        xlims       = [-35 30];
        ylims       = [-30 30];
        fontSize    = [];
        fig_axesFormat(gca, xTick,xTickLabels, yTick,yTickLabels, xAxisTitle,yAxisTitle, xlims,ylims, fontSize);

        % Plot 2
        drawnow();
        fig_subplot(2);
        plot(existingTestLocations_XY(:,1), existingTestLocations_XY(:,2), '.w')
        hold on
        for ii = 1:size(existingTestLocations_XY,1)
            DLS_dB = existingTestLocations_DLS_dB(ii);
            txt = sprintf('\\bf%1.0f', DLS_dB);
            text(existingTestLocations_XY(ii,1), existingTestLocations_XY(ii,2), txt, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 7, 'Interpreter','tex');
        end
        for ii = 1:size(newTestLocations_XY,1)
            DLS_dB = newTestLocations_DLS_dB(ii);
            txt = sprintf('\\bf%1.0f', DLS_dB);
            text(newTestLocations_XY(ii,1), newTestLocations_XY(ii,2), txt, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 7, 'Color','r', 'Interpreter','tex');
        end
        axis square
        plot(xlim(), [0 0], 'k--', [0 0], ylim(), 'k--'); % crosshair
        set(gca, 'xlim',[-30 30], 'ylim',[-30 30]);

        % title
        hTxt = textLoc('\bf DLS Estimates', 'NorthWest', 'Interpreter','tex');
        fig_nudge(hTxt, -.075, .075);

        % format the axes
        fig_axesFormat(gca, xTick,xTickLabels, yTick,yTickLabels, xAxisTitle,yAxisTitle, xlims,ylims, fontSize);

        %%format the figure
        xTitle      = ['X' char(176)];
        yTitle      = ['Y' char(176)];
        mainTitle   = [];
        fontSize    = 16;
        [hXTitle,hYTitle] = fig_figFormat(hFig, xTitle,yTitle,mainTitle, fontSize, false);
        fig_nudge(hXTitle, 0, .3);
        fig_nudge(hYTitle, .3, 0);

        % save
        fn = sprintf('sim1_example_%i', i);
        fig_save(hFig, fn, EXPORT_DIR, EXPORT_FORMAT);
        close(hFig);
    end


% =======================================================================
%% 1.2 Plot simulated observer
% =======================================================================
% close all

    % plot 30-2
    DLS = simObs_true_DLS;
    isRightEye = true;
    doAnnotation = false;
    age = 18;
    maskmatrix = [];
    hFig1 = VfPlot.plotWithNums(DLS, isRightEye, doAnnotation, age, maskmatrix, x_deg, y_deg);
    
    
    % plot 24-2
    ids = [ % 24-2, right eye format (OD)
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN % +27
        NaN, NaN, NaN,   1,   2,   3,   4, NaN, NaN, NaN % +21
        NaN, NaN,   5,   6,   7,   8,   9,  10, NaN, NaN % +15
        NaN,  11,  12,  13,  14,  15,  16,  17,  18, NaN % +9
        19,   20,  21,  22,  23,  24,  25,  26,  27, NaN % +3
        28,   29,  30,  31,  32,  33,  34,  35,  36, NaN % -3
        NaN,  37,  38,  39,  40,  41,  42,  43,  44, NaN % -9
        NaN, NaN,  45,  46,  47,  48,  49,  50, NaN, NaN % -15
        NaN, NaN, NaN,  51,  52,  53,  54, NaN, NaN, NaN % -21
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN % -27
        ];  %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
    
    DLS = simObs_true_DLS();
    DLS(isnan(ids)) = NaN;
    DLS([1 end],:) = [];
    %
    x = x_deg;
    x(isnan(ids)) = NaN;
    x([1 end],:) = [];
    %
    y = y_deg;
    y(isnan(ids)) = NaN;
    y([1 end],:) = [];
    %
    isRightEye = true;
    doAnnotation = false;
    age = 18;
    maskmatrix = [];
    hFig2 = VfPlot.plotWithNums(DLS, isRightEye, doAnnotation, age, maskmatrix, x, y);
    
    % save
    fig_save(hFig1, 'sim1_simObserver_30-2', EXPORT_DIR, EXPORT_FORMAT);
    fig_save(hFig2, 'sim1_simObserver_24-2', EXPORT_DIR, EXPORT_FORMAT);


% =========================================================================
% 1.3 analysis
% =========================================================================
% close all
    
    % Plotting ---------------------------------------------------------

    % open a new figure window.
    figDim = [10 10];
    nPlots = 1;
    isTight = [];
    isSquare = true;
    styleFlag = [];
    axesLims = [];
    hFig = fig_make(figDim, nPlots, isTight, isSquare, styleFlag, axesLims);

    % plot
    dat = nLocsInScotoma;
    errorbar(nExtraLocations, mean(dat), diff(bootci(N_BOOTSTRP,{@mean,dat})), 'ko');

    % run repeated measures ANOVA
    data = array2table(nLocsInScotoma, 'VariableNames', {'n0', 'n1', 'n2', 'n3', 'n4', 'n5'});
    withinDesign = table([1 2 3 4 5 6]', 'VariableNames', {'nExtraLocs'});% Create a within-subjects design table indicating the measurements are repeated
    rmModel = fitrm(data, 'n0-n5~1', 'WithinDesign', withinDesign);
    rmANOVA = ranova(rmModel);
    % add stats
    txt = fStr(rmANOVA.DF(1), rmANOVA.DF(2), rmANOVA.F(1), rmANOVA.pValue(1));
    hTxt = textLoc(txt, 'NorthWest');
    fig_nudge(hTxt, -0.05, 0.05);

    % format all the axes
    xTick   = 0:5;
    xTickLabels = [];
    yTick = 3:8;
    yTickLabels = [];
    xTitle = '$N$ Additional Test Locations';
    yTitle = '$N$ test points within scotoma';
    xlims = [-.5 5.5];
    ylims = [2.5 8.5];
    fig_axesFormat(gca, xTick,xTickLabels, yTick,yTickLabels, xTitle,yTitle, xlims,ylims);
    
    % format the figure
    xTitle = [];
    yTitle = [];
    mainTitle = [];
    fontSize = 16;
    fig_figFormat(hFig, xTitle,yTitle,mainTitle, fontSize, false);

    %% save
    fig_save(hFig, 'sim1_results', EXPORT_DIR, EXPORT_FORMAT);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIM 2: Extensive loss in 1 hemifield -- example 2 in MEDTEG.runExamples()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ratio of points in the preserved/damaged hemifield, as function of n
% points (vline where the fixed seed points end)
fprintf('\n\n ***** SIM 2 ***** \n\n')

% =========================================================================
% 2.1 run simulations
% =========================================================================

    % params
    N_SIMS = N_SIM(2);
    
    % init results data storage
    ratio = nan(N_SIMS, 1);
    ratio_newLocs = nan(N_SIMS, 1);
    nInSpared_newLocs = nan(N_SIMS, 1);
    
    % run
    for i = 1:N_SIMS
        close all

        try
            [x_deg, y_deg, simObs_true_DLS, F_simObs_DLS, dat] = MEDTEG.runExample(2); % dat: [x y norm_DLS observer_true_DLS observer_est_DLS isNewlyAddedPoint]
        catch ME
            warning(ME.message);
            continue
        end


        DLS_upper = nanmean(simObs_true_DLS(y_deg>0)); %#ok
        DLS_lower = nanmean(simObs_true_DLS(y_deg<0)); %#ok
        isUpperHemifieldLoss = DLS_upper < DLS_lower;
    
        
        idx = dat.isNewlyAddedPoint==1;
        y = dat.y_deg(idx);
    
        nInUpperOriginal = sum(dat.y_deg(~idx)>0);
        nInLowerOriginal = sum(dat.y_deg(~idx)>0);
    
        for j = 1:length(y)
            nInUpper = [sum(y(1:j)>0) sum(y(1:j)>0)+nInUpperOriginal]; % [ new_only, all]
            nInLower = [sum(y(1:j)<0) sum(y(1:j)<0)+nInLowerOriginal]; % [ new_only, all]
    
            if isUpperHemifieldLoss
                nInAffected = nInUpper;
                nInSpared = nInLower;
            else
                nInAffected = nInLower;
                nInSpared = nInUpper;
            end
    
            ratio_newLocs(i,j) = nInSpared(1)/nInAffected(1);
            nInSpared_newLocs(i,j) = nInSpared(1);
            ratio(i,j) = nInSpared(2)/nInAffected(2);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot first N individual examples
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % skip if plotted enough already
        if (i > N_IND_SIMS_TO_PLOT)
            continue; % skip
        end

        % open a new figure window
        figDim      = [15 20];
        nPlots      = [1 2];
        isTight     = true;
        isSquare    = true;
        styleFlag   = [];
        axesLims    = [];
        hFig        = fig_make(figDim, nPlots, isTight, isSquare, styleFlag, axesLims);

        % extract data
        idx = dat.isNewlyAddedPoint==1;
        existingTestLocations_XY = [dat.x_deg(~idx) dat.y_deg(~idx)];
        existingTestLocations_DLS_dB = dat.observer_est_DLS(~idx);
        newTestLocations_XY = [dat.x_deg(idx) dat.y_deg(idx)];
        newTestLocations_DLS_dB = dat.observer_est_DLS(idx);

        % Plot 1
        fig_subplot(1);
        hold on
        plot(existingTestLocations_XY(:,1), existingTestLocations_XY(:,2), 'ko', 'markerfacecolor',[.8 .8 .8])
        plot([15 15],[-3 3],'kx');
        for ii = 1:length(newTestLocations_DLS_dB)
            txt = sprintf('%i', ii);
            text(newTestLocations_XY(ii,1), newTestLocations_XY(ii,2), txt, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 7, 'Color','k', 'Interpreter','tex');
        end
        axis square
        plot(xlim(), [0 0], 'k--', [0 0], ylim(), 'k--'); % crosshair

        % title
        hTxt = textLoc('\bf Test Locations / Order', 'NorthWest', 'Interpreter','tex');
        fig_nudge(hTxt, -.075, .075);

        % format the axes
        xTick       = -20:10:20;
        xTickLabels = [];
        yTick       = xTick;
        yTickLabels = [];
        xAxisTitle  = [];
        yAxisTitle  = [];
        xlims       = [-35 30];
        ylims       = [-30 30];
        fontSize    = [];
        fig_axesFormat(gca, xTick,xTickLabels, yTick,yTickLabels, xAxisTitle,yAxisTitle, xlims,ylims, fontSize);

        % Plot 2
        drawnow();
        fig_subplot(2);
        plot(existingTestLocations_XY(:,1), existingTestLocations_XY(:,2), '.w')
        hold on
        for ii = 1:size(existingTestLocations_XY,1)
            DLS_dB = existingTestLocations_DLS_dB(ii);
            txt = sprintf('\\bf%1.0f', DLS_dB);
            text(existingTestLocations_XY(ii,1), existingTestLocations_XY(ii,2), txt, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 7, 'Interpreter','tex');
        end
        for ii = 1:size(newTestLocations_XY,1)
            DLS_dB = newTestLocations_DLS_dB(ii);
            txt = sprintf('\\bf%1.0f', DLS_dB);
            text(newTestLocations_XY(ii,1), newTestLocations_XY(ii,2), txt, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 7, 'Color','r', 'Interpreter','tex');
        end
        axis square
        plot(xlim(), [0 0], 'k--', [0 0], ylim(), 'k--'); % crosshair
        set(gca, 'xlim',[-30 30], 'ylim',[-30 30]);

        % title
        hTxt = textLoc('\bf DLS Estimates', 'NorthWest', 'Interpreter','tex');
        fig_nudge(hTxt, -.075, .075);

        % format the axes
        fig_axesFormat(gca, xTick,xTickLabels, yTick,yTickLabels, xAxisTitle,yAxisTitle, xlims,ylims, fontSize);

        %%format the figure
        xTitle      = ['X' char(176)];
        yTitle      = ['Y' char(176)];
        mainTitle   = [];
        fontSize    = 16;
        [hXTitle,hYTitle] = fig_figFormat(hFig, xTitle,yTitle,mainTitle, fontSize, false);
        fig_nudge(hXTitle, 0, .3);
        fig_nudge(hYTitle, .3, 0);

        % save
        fn = sprintf('sim2_example_%i', i);
        fig_save(hFig, fn, EXPORT_DIR, EXPORT_FORMAT);
        close(hFig);
    end


% =========================================================================
% 2.2 Plot simulated observer
% =========================================================================
% close all

    % plot 30-2
    DLS = simObs_true_DLS;
    isRightEye = true;
    doAnnotation = false;
    age = 18;
    maskmatrix = [];
    hFig1 = VfPlot.plotWithNums(DLS, isRightEye, doAnnotation, age, maskmatrix, x_deg, y_deg);
    
    
    % plot 24-2
    ids = [ % 24-2, right eye format (OD)
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN % +27
        NaN, NaN, NaN,   1,   2,   3,   4, NaN, NaN, NaN % +21
        NaN, NaN,   5,   6,   7,   8,   9,  10, NaN, NaN % +15
        NaN,  11,  12,  13,  14,  15,  16,  17,  18, NaN % +9
        19,   20,  21,  22,  23,  24,  25,  26,  27, NaN % +3
        28,   29,  30,  31,  32,  33,  34,  35,  36, NaN % -3
        NaN,  37,  38,  39,  40,  41,  42,  43,  44, NaN % -9
        NaN, NaN,  45,  46,  47,  48,  49,  50, NaN, NaN % -15
        NaN, NaN, NaN,  51,  52,  53,  54, NaN, NaN, NaN % -21
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN % -27
        ];  %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
    
    DLS = simObs_true_DLS();
    DLS(isnan(ids)) = NaN;
    DLS([1 end],:) = [];
    %
    x = x_deg;
    x(isnan(ids)) = NaN;
    x([1 end],:) = [];
    %
    y = y_deg;
    y(isnan(ids)) = NaN;
    y([1 end],:) = [];
    %
    isRightEye = true;
    doAnnotation = false;
    age = 18;
    maskmatrix = [];
    hFig2 = VfPlot.plotWithNums(DLS, isRightEye, doAnnotation, age, maskmatrix, x, y);
    
    % save
    fig_save(hFig1, 'sim2_simObserver_30-2', EXPORT_DIR, EXPORT_FORMAT);
    fig_save(hFig2, 'sim2_simObserver_24-2', EXPORT_DIR, EXPORT_FORMAT);

% =========================================================================
% 2.3 analysis
% =========================================================================
% close all

    % open a new figure window
    figDim      = [9.5 11];
    nPlots      = 1;
    isTight     = true;
    isSquare    = true;
    styleFlag   = [];
    axesLims    = [];
    hFig        = fig_make(figDim, nPlots, isTight, isSquare, styleFlag, axesLims);

    % plot data
    errorbar(1:14, mean(ratio), diff(bootci(N_BOOTSTRP,{@mean,ratio})), 'ko', 'markerfacecolor',[.8 .8 1]);

    % format the axes
    xTick       = 0:2:15;
    xTickLabels = [];
    yTick       = 1:.15:1.6;
    yTickLabels = [];
    xAxisTitle  = '$N$ Additional Test Locations';
    yAxisTitle  = {'Ratio of total test locations', 'in each hemifield $\frac{N_{healthy}}{N_{impaired}}$'}; % 'Ratio $\frac{Healthy \;\; Hemifield}{Impaired \;\; Hemifield}$';
    xlims       = [-1 15];
    ylims       = [.9 1.6];
    fontSize    = [];
    [~, ~, ~, hxAxisTitle] = fig_axesFormat(gca, xTick,xTickLabels, yTick,yTickLabels, xAxisTitle,yAxisTitle, xlims,ylims, fontSize);
    fig_nudge(hxAxisTitle, 0, 0.1);

    % add annotation
    hline(1);

    % add secondary x axis: note, currently comes *after* figFormat(!)
    xTickLabels = xTick + 10;
    xAxisTitle = '$N$ Location Total (incl. 10 seed points)';
    fontSize = [10 9];
    hTickLbl = fig_addSecondAxis(gca, 'top', xTick, xTickLabels, xAxisTitle, fontSize);
    fig_nudge(hTickLbl, 0, -0.02);

    % add secondary y axis
    % NB: ratio == (nInSpared_newLocs+5)./(repmat(1:14,N_SIMS,1)-nInSpared_newLocs+5)
    


    %%format the figure
    xTitle      = [];
    yTitle      = [];
    mainTitle   = [];
    fontSize    = 16;
    fig_figFormat(hFig, xTitle,yTitle,mainTitle, fontSize, false);

    % save
    fig_save(hFig, 'sim2_results', EXPORT_DIR, EXPORT_FORMAT);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIM 3: Integrating structural data -- example 3 in MEDTEG.runExamples()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n\n ***** SIM 3 ***** \n\n')

% =========================================================================
% 3.1 run simulations
% =========================================================================

    % params
    N_SIMS = N_SIM(3);
    condNames = {'No noise; veridical', 'w/ noise; veridical', 'No noise; x-flipped', 'w/ noise; x-flipped'};
    noiseMode = {'none', 'high', 'none', 'high'};
    veridicalMode = [true true false false];
    nConditions = 4;
    
    % init results data storage
    nNewLocsInScotoma = nan(N_SIMS, nConditions); % [nonoise-veridical lownoise-veridical highnoise-veridical nonoise-false lownoise-false highnoise-false]

    % run
    for i = 1:N_SIMS
        close all

        for j = 1:nConditions
            % run
            try
                [x_deg, y_deg, simObs_true_DLS, F_simObs_DLS, dat] =  MEDTEG.runExample(3, struct('structNoiseLevel',noiseMode{j}, 'structIsVerical',veridicalMode(j))); % dat: [x y norm_DLS observer_true_DLS observer_est_DLS isNewlyAddedPoint]
            catch ME
                warning(ME.message);
                continue
            end

            % store result
            TD = dat.observer_true_DLS - dat.norm_DLS;
            inScotoma = TD < -8;
            inBlindSpot = dat.x_deg > 13 & dat.x_deg < 17 & dat.y_deg > -4 & dat.y_deg < 4; % crude but good enough
            isNew = dat.isNewlyAddedPoint;
            idx = inScotoma & ~inBlindSpot & dat.isNewlyAddedPoint==1;
            nNewLocsInScotoma(i, j) = sum(idx);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Lattice plot first N individual examples
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % skip if plotted enough already
            if (i > N_IND_SIMS_TO_PLOT)
                continue; % skip
            end

            % open a new figure window
            if (j == 1)
                figDim      = [15 20];
                nPlots      = [2 2];
                isTight     = true;
                isSquare    = true;
                styleFlag   = [];
                axesLims    = [];
                hFig        = fig_make(figDim, nPlots, isTight, isSquare, styleFlag, axesLims);
            else
                figure(hFig);
            end

            % extract data
            idx = dat.isNewlyAddedPoint==1;
            existingTestLocations_XY = [dat.x_deg(~idx) dat.y_deg(~idx)];
            existingTestLocations_DLS_dB = dat.observer_est_DLS(~idx);
            newTestLocations_XY = [dat.x_deg(idx) dat.y_deg(idx)];
            newTestLocations_DLS_dB = dat.observer_est_DLS(idx);

            % Plot DLS
            fig_subplot(j);
            plot(existingTestLocations_XY(:,1), existingTestLocations_XY(:,2), '.w')
            hold on
            for ii = 1:size(existingTestLocations_XY,1)
                DLS_dB = existingTestLocations_DLS_dB(ii);
                txt = sprintf('\\bf%1.0f', DLS_dB);
                text(existingTestLocations_XY(ii,1), existingTestLocations_XY(ii,2), txt, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 7, 'Interpreter','tex');
            end
            for ii = 1:size(newTestLocations_XY,1)
                DLS_dB = newTestLocations_DLS_dB(ii);
                txt = sprintf('\\bf%1.0f', DLS_dB);
                text(newTestLocations_XY(ii,1), newTestLocations_XY(ii,2), txt, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 7, 'Color','r', 'Interpreter','tex');
            end
            axis square
            set(gca, 'xlim',[-30 30], 'ylim',[-30 30]);

            % title
            txt = condNames{j}; % ['\bf' condNames{j}]
            hTxt = textLoc(txt, 'NorthWest', 'Interpreter','tex');
            fig_nudge(hTxt, -.075, .075);

            % format the axes
            xTick       = -20:10:20;
            xTickLabels = [];
            yTick       = xTick;
            yTickLabels = [];
            xAxisTitle  = [];
            yAxisTitle  = [];
            xlims       = [-30 30];
            ylims       = [-30 30];
            fontSize    = [14 10];
            fig_axesFormat(gca, xTick,xTickLabels, yTick,yTickLabels, xAxisTitle,yAxisTitle, xlims,ylims, fontSize);

            % add crosshairs
            plot(xlim(), [0 0], 'k--', [0 0], ylim(), 'k--'); % crosshair

            if (j == nConditions)
                %%format the figure
                xTitle      = ['X' char(176)];
                yTitle      = ['Y' char(176)];
                mainTitle   = [];
                fontSize    = [14 10];
                [hXTitle,hYTitle] = fig_figFormat(hFig, xTitle,yTitle,mainTitle, fontSize, false);
                fig_nudge(hXTitle, 0, .3);
                fig_nudge(hYTitle, .3, 0);

                % save
                fn = sprintf('sim3_example%i', i);
                fig_save(hFig, fn, EXPORT_DIR, EXPORT_FORMAT);

                % close
                close(hFig);
            end
        end
    end

% =========================================================================
% 3.2 Plot simulated observer
% =========================================================================
% close all

    % plot 30-2
    DLS = simObs_true_DLS;
    isRightEye = true;
    doAnnotation = false;
    age = 18;
    maskmatrix = [];
    hFig1 = VfPlot.plotWithNums(DLS, isRightEye, doAnnotation, age, maskmatrix, x_deg, y_deg);
    
    
    % plot 24-2
    ids = [ % 24-2, right eye format (OD)
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN % +27
        NaN, NaN, NaN,   1,   2,   3,   4, NaN, NaN, NaN % +21
        NaN, NaN,   5,   6,   7,   8,   9,  10, NaN, NaN % +15
        NaN,  11,  12,  13,  14,  15,  16,  17,  18, NaN % +9
        19,   20,  21,  22,  23,  24,  25,  26,  27, NaN % +3
        28,   29,  30,  31,  32,  33,  34,  35,  36, NaN % -3
        NaN,  37,  38,  39,  40,  41,  42,  43,  44, NaN % -9
        NaN, NaN,  45,  46,  47,  48,  49,  50, NaN, NaN % -15
        NaN, NaN, NaN,  51,  52,  53,  54, NaN, NaN, NaN % -21
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN % -27
        ];  %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
    
    DLS = simObs_true_DLS();
    DLS(isnan(ids)) = NaN;
    DLS([1 end],:) = [];
    %
    x = x_deg;
    x(isnan(ids)) = NaN;
    x([1 end],:) = [];
    %
    y = y_deg;
    y(isnan(ids)) = NaN;
    y([1 end],:) = [];
    %
    isRightEye = true;
    doAnnotation = false;
    age = 18;
    maskmatrix = [];
    hFig2 = VfPlot.plotWithNums(DLS, isRightEye, doAnnotation, age, maskmatrix, x, y);
    
    % save
    fig_save(hFig1, 'sim3_simObserver_30-2', EXPORT_DIR, EXPORT_FORMAT);
    fig_save(hFig2, 'sim3_simObserver_24-2', EXPORT_DIR, EXPORT_FORMAT);

% =========================================================================
% 3.3 analysis
% =========================================================================
% close all

    % open a new figure window
    figDim      = [10 9];
    nPlots      = 1;
    isTight     = true;
    isSquare    = true;
    styleFlag   = [];
    axesLims    = [];
    hFig        = fig_make(figDim, nPlots, isTight, isSquare, styleFlag, axesLims);

    % plot data
    colors = [.8 .8 1; .6 .6 1; 1 .8 .8; 1 .6 .6];
    hold on 
    xnudge = [0 0 .5 .5];
    for i = 1:4
        bar(i+xnudge(i), nanmean(nNewLocsInScotoma(:,i)), 'facecolor',colors(i,:)) %#ok
    end
    errorbar((1:4)+xnudge, nanmean(nNewLocsInScotoma), diff(bootci(N_BOOTSTRP,{@nanmean, nNewLocsInScotoma})), 'ko', 'markerfacecolor',[.8 .8 .8]); %#ok

    % format the axes
    xTick       = (1:4)+xnudge;
    xTickLabels = {'no\; noise','noise','no\; noise','noise'};
    yTick       = 0:2:12;
    yTickLabels = [];
    xAxisTitle  = '\textit{veridical} \qquad \textit{x-flipped}\quad';
    yAxisTitle  = [];
    xlims       = [];
    ylims       = [];
    fontSize    = [];
    formatData    = [];
    xMinorTick    = [];
    yMinorTick    = [];
    mrTickLgth    = [];
    xRotation    = 45;
    yRotation    = [];
    [~, ~, ~, hxAxisTitle] = fig_axesFormat(gca, xTick,xTickLabels, yTick,yTickLabels, xAxisTitle,yAxisTitle, xlims,ylims, fontSize, formatData, xMinorTick,yMinorTick,mrTickLgth, xRotation,yRotation);
    fig_nudge(hxAxisTitle, 0, .4);

    %%format the figure
    xTitle      = 'Structural Data';
    yTitle      = '$N$ new points within scotoma';
    mainTitle   = [];
    fontSize    = 16;
    [hXTitle,~,~] = fig_figFormat(hFig, xTitle,yTitle,mainTitle, fontSize, false);
    fig_nudge(hXTitle, 0, .07);

    
    % add secondary axes: note, currently comes *after* figFormat(!)
    drawnow();
    yTickLabels = round(yTick ./ 14 * 100);
    yAxisTitle = '%';
    fontSize = [10 9];
    [hTickLbl, hAxisLabel] = fig_addSecondAxis(gca, 'right', yTick, yTickLabels, yAxisTitle, fontSize);
    fig_nudge(hAxisLabel, -0.1, 0);
    fig_expandFigureWindow(hFig, 0.5, 0);

    % save
    fig_save(hFig, 'sim3_results', EXPORT_DIR, EXPORT_FORMAT);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIM 4: Early termination (inattentiveness) -- example 4 in MEDTEG.runExamples()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n\n ***** SIM 4 ***** \n\n')

% =========================================================================
% 4.1 run simulations
% =========================================================================

    % params
    N_SIMS = N_SIM(4);
    
    % init results data storage
    SSE_diffFromGroundTruth = nan(N_SIMS, 3); % [nLocs 24-2 MEDTEG]

    % create a "ground truth" to validate against
    id = [ % 24-2, right eye format (OD)
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN % +27
        NaN, NaN, NaN,   1,   2,   3,   4, NaN, NaN, NaN % +21
        NaN, NaN,   5,   6,   7,   8,   9,  10, NaN, NaN % +15
        NaN,  11,  12,  13,  14,  15,  16,  17,  18, NaN % +9
        19,   20,  21,  22,  23,  24,  25,  26,  27, NaN % +3
        28,   29,  30,  31,  32,  33,  34,  35,  36, NaN % -3
        NaN,  37,  38,  39,  40,  41,  42,  43,  44, NaN % -9
        NaN, NaN,  45,  46,  47,  48,  49,  50, NaN, NaN % -15
        NaN, NaN, NaN,  51,  52,  53,  54, NaN, NaN, NaN % -21
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN % -27
        ];  %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
    [x, y] = meshgrid(linspace(-27, 27, 10), linspace(27, -27, 10));
    x = x(~isnan(id));
    y = y(~isnan(id));
    k = convhull([x y]);
    chverts = [x(k) y(k)];
    % create a "ground truth" to validate against
    N = 100;
    [x1,y1] = meshgrid(linspace(min(x), max(x), N), linspace(min(y), max(y), N));


    % remove outside convex hull
    % vals in matrix form to be retained for use in surf() later
    idx = ~inhull([x1(:) y1(:)], chverts, [], 0.1);
    x1(idx) = NaN;
    y1(idx) = NaN;

    % remove around blindspot
    inBlindSpot = x1 > 12 & x1 < 18 & y1 > -5 & y1 < 5; % crude but good enough
    x1(inBlindSpot) = NaN;
    y1(inBlindSpot) = NaN;

    % vectorised version, without NaNs
    x_hidensity = x1(~isnan(x1));
    y_hidensity = y1(~isnan(y1));

    % run
    for i = 1:N_SIMS
        % run
        try
            [x_deg, y_deg, simObs_true_DLS, F_simObs_DLS, dat] =  MEDTEG.runExample(4); % dat: [x y norm_DLS observer_true_DLS observer_est_DLS isNewlyAddedPoint isStandard]
        catch ME
            warning(ME.message);
            continue
        end

        % extract key metrics
        nLocs = height(dat)/2;

        % ground truth
        DLS_true = F_simObs_DLS(x_hidensity, y_hidensity);

        % 24-2
        idx = 1:nLocs;
        x_242 = dat.x_deg(idx);
        y_242 = dat.y_deg(idx);
        v_242 = dat.observer_est_DLS(idx);
        F_242_DLS = scatteredInterpolant(x_242, y_242, v_242);
        %
        DLS_242_estimates = F_242_DLS(x_hidensity, y_hidensity);
        sum_squared_residuals_242 = sum((DLS_true - DLS_242_estimates).^2);

        % MEDTEG
        idx = (nLocs+1):(nLocs*2);
        x_MEDTEG = dat.x_deg(idx);
        y_MEDTEG = dat.y_deg(idx);
        v_MEDTEG = dat.observer_est_DLS(idx);
        F_MEDTEG_DLS = scatteredInterpolant(x_MEDTEG, y_MEDTEG, v_MEDTEG);
        %
        DLS_MEDTEG_estimates = F_MEDTEG_DLS(x_hidensity, y_hidensity);
        sum_squared_residuals_MEDTEG = sum((DLS_true - DLS_MEDTEG_estimates).^2);

        % store result
        SSE_diffFromGroundTruth(i,:) = [nLocs, sum_squared_residuals_242, sum_squared_residuals_MEDTEG];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Lattice plot first N individual examples
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % skip if plotted enough already
        if (i > N_IND_SIMS_TO_PLOT)
            continue; % skip
        end

        % open a new figure window
        figDim      = [20 20];
        nPlots      = [2 2];
        isTight     = true;
        isSquare    = true;
        styleFlag   = [];
        axesLims    = [];
        hFig        = fig_make(figDim, nPlots, isTight, isSquare, styleFlag, axesLims);

        % =================================================================
        % Plot DLS (24-2)
        fig_subplot(1,1);
        hold on
        for ii = 1:length(v_242)
            txt = sprintf('\\bf%1.0f', v_242(ii));
            text(x_242(ii), y_242(ii), txt, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 7, 'Interpreter','tex');
        end

        % format the axes
        axis square
        xTick       = -20:10:20;
        xTickLabels = [];
        yTick       = xTick;
        yTickLabels = [];
        xAxisTitle  = [];
        yAxisTitle  = [];
        xlims       = [-30 30];
        ylims       = [-30 30];
        fontSize    = [];
        fig_axesFormat(gca, xTick,xTickLabels, yTick,yTickLabels, xAxisTitle,yAxisTitle, xlims,ylims, fontSize);

        % title
        hTxt = textLoc('24-2', 'NorthWest', 'Interpreter','tex');
        fig_nudge(hTxt, -.05, .05);

        % add crosshairs
        plot(xlim(), [0 0], 'k--', [0 0], ylim(), 'k--'); % crosshair

        % plot surface
        drawnow();
        fig_subplot(2,1);
        surf(x1, y1, F_242_DLS(x1, y1));
        fig_axesFormat(gca, xTick,xTickLabels, yTick,yTickLabels, xAxisTitle,yAxisTitle, xlims,ylims, fontSize);

        % label with error
        txt = sprintf('SSR = %1.1f', SSE_diffFromGroundTruth(i, 2));
        hTxt = textLoc(txt, 'NorthWest');
        fig_nudge(hTxt, -0.05, 0.05);

        % =================================================================
        % Plot DLS (MEDTEG)
        fig_subplot(1,2);
        hold on
        for ii = 1:length(v_MEDTEG)
            txt = sprintf('\\bf%1.0f', v_MEDTEG(ii));
            hTxt = text(x_MEDTEG(ii), y_MEDTEG(ii), txt, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 7, 'Interpreter','tex');
            if (ii > 10)
                set(hTxt, 'Color','r');
            end
        end

        % format the axes
        axis square
        fig_axesFormat(gca, xTick,xTickLabels, yTick,yTickLabels, xAxisTitle,yAxisTitle, xlims,ylims, fontSize);

        % title
        hTxt = textLoc('MEDTEG', 'NorthWest', 'Interpreter','tex');
        fig_nudge(hTxt, -.05, .05);

        % add crosshairs
        plot(xlim(), [0 0], 'k--', [0 0], ylim(), 'k--'); % crosshair

        % plot surface
        fig_subplot(2,2);
        surf(x1, y1, F_MEDTEG_DLS(x1, y1));
        fig_axesFormat(gca, xTick,xTickLabels, yTick,yTickLabels, xAxisTitle,yAxisTitle, xlims,ylims, fontSize);

        % label with error
        txt = sprintf('SSR = %1.1f', SSE_diffFromGroundTruth(i, 3));
        hTxt = textLoc(txt, 'NorthWest');
        fig_nudge(hTxt, -0.05, 0.05);

        % =================================================================
        %%format the figure
        xTitle      = ['X' char(176)];
        yTitle      = ['Y' char(176)];
        mainTitle   = [];
        fontSize    = 16;
        [hXTitle,hYTitle] = fig_figFormat(hFig, xTitle,yTitle,mainTitle, fontSize, false);
        fig_nudge(hXTitle, 0, .3);
        fig_nudge(hYTitle, .3, 0);

        % save
        fn = sprintf('sim4_example%i', i);
        fig_save(hFig, fn, EXPORT_DIR, EXPORT_FORMAT);

        % close
        close(hFig);
    end

% =========================================================================
% 4.2 Plot simulated observer
% =========================================================================
% close all

    % plot 30-2
    DLS = simObs_true_DLS;
    isRightEye = true;
    doAnnotation = false;
    age = 18;
    maskmatrix = [];
    hFig1 = VfPlot.plotWithNums(DLS, isRightEye, doAnnotation, age, maskmatrix, x_deg, y_deg);
    
    
    % plot 24-2
    ids = [ % 24-2, right eye format (OD)
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN % +27
        NaN, NaN, NaN,   1,   2,   3,   4, NaN, NaN, NaN % +21
        NaN, NaN,   5,   6,   7,   8,   9,  10, NaN, NaN % +15
        NaN,  11,  12,  13,  14,  15,  16,  17,  18, NaN % +9
        19,   20,  21,  22,  23,  24,  25,  26,  27, NaN % +3
        28,   29,  30,  31,  32,  33,  34,  35,  36, NaN % -3
        NaN,  37,  38,  39,  40,  41,  42,  43,  44, NaN % -9
        NaN, NaN,  45,  46,  47,  48,  49,  50, NaN, NaN % -15
        NaN, NaN, NaN,  51,  52,  53,  54, NaN, NaN, NaN % -21
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN % -27
        ];  %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
    
    DLS = simObs_true_DLS();
    DLS(isnan(ids)) = NaN;
    DLS([1 end],:) = [];
    %
    x = x_deg;
    x(isnan(ids)) = NaN;
    x([1 end],:) = [];
    %
    y = y_deg;
    y(isnan(ids)) = NaN;
    y([1 end],:) = [];
    %
    isRightEye = true;
    doAnnotation = false;
    age = 18;
    maskmatrix = [];
    hFig2 = VfPlot.plotWithNums(DLS, isRightEye, doAnnotation, age, maskmatrix, x, y);
    
    % save
    fig_save(hFig1, 'sim4_simObserver_30-2', EXPORT_DIR, EXPORT_FORMAT);
    fig_save(hFig2, 'sim4_simObserver_24-2', EXPORT_DIR, EXPORT_FORMAT);

% =========================================================================
% 4.3 analysis
% =========================================================================
% close all

    % open a new figure window
    figDim      = [10 20];
    nPlots      = [1 2];
    isTight     = true;
    isSquare    = true;
    styleFlag   = [];
    axesLims    = [];
    hFig        = fig_make(figDim, nPlots, isTight, isSquare, styleFlag, axesLims);

    % make 2nd subplot half-width
    fig_subplot(2);
    set(gca,'Position',get(gca,'Position').*[1 1 .5 1]);

    % =====================================================================
    fig_subplot(1);

    % plot raw
    hold on
    plot(SSE_diffFromGroundTruth(:,1), SSE_diffFromGroundTruth(:,2), 'bo');
    plot(SSE_diffFromGroundTruth(:,1), SSE_diffFromGroundTruth(:,3), 'rs');

    % plot fit
    idx = ~isnan(SSE_diffFromGroundTruth(:,1));
    fitobj = fit(SSE_diffFromGroundTruth(idx,1), SSE_diffFromGroundTruth(idx,2), 'log10' );
    xFit = linspace(min(SSE_diffFromGroundTruth(:,1)), max(SSE_diffFromGroundTruth(:,1)), 200);
    yFit = fitobj(xFit);
    plot(xFit, yFit, 'b-', 'linewidth',2);
    %
    fitobj = fit(SSE_diffFromGroundTruth(idx,1), SSE_diffFromGroundTruth(idx,3), 'log10');
    yFit = fitobj(xFit);
    plot(xFit, yFit, 'r--', 'linewidth',2);

    % format the axes
    xTick       = 15:10:45;
    xTickLabels = [];
    yTick       = 0:50000:350000;
    yTickLabels = [];
    xAxisTitle  = '$N$ Locations Total';
    yAxisTitle  = 'Sum Squared Residuals (SSR)';
    xlims       = [];
    ylims       = [0 yTick(end)];
    fontSize    = [];
    fig_axesFormat(gca, xTick,xTickLabels, yTick,yTickLabels, xAxisTitle,yAxisTitle, xlims,ylims, fontSize);
    
    % add legend
    hDat = plot(NaN, NaN, 'bo-', NaN, NaN, 'rs--');
    legend(hDat, '24-2', 'MEDTEG')

    % =====================================================================
    fig_subplot(2);

    % plot data
    hold on
    %bar([1 2], nanmean(X));
    idx = ~any(isnan(SSE_diffFromGroundTruth),2);
    errorbar(1, mean(SSE_diffFromGroundTruth(idx,2)), diff(bootci(N_BOOTSTRP,{@mean, SSE_diffFromGroundTruth(idx,2)})), 'bo', 'linewidth',2 , 'markerfacecolor',[.8 .8 .8]);
    errorbar(2, mean(SSE_diffFromGroundTruth(idx,3)), diff(bootci(N_BOOTSTRP,{@mean, SSE_diffFromGroundTruth(idx,3)})), 'rs', 'linewidth',2, 'markerfacecolor',[.8 .8 .8]);

    
    % format the axes
    xTick       = [1 2];
    xTickLabels = {'24-2\;','MEDTEG\;'};
    yTick       = [];
    yTickLabels = [];
    xAxisTitle  = [];
    yAxisTitle  = [];
    xlims       = [0.33 2.66];
    %ylims       = [];
    formatData  = [];
    xMinorTick  = [];
    yMinorTick  = [];
    mrTickLgth  = [];
    xRotation   = 45;
    yRotation   = [];
    [hXTickLbl, hYTickLbl, c_axes, hxAxisTitle] = fig_axesFormat(gca, xTick,xTickLabels, yTick,yTickLabels, xAxisTitle,yAxisTitle, xlims,ylims, fontSize, formatData, xMinorTick,yMinorTick,mrTickLgth, xRotation,yRotation);
    fig_nudge(hXTickLbl, .2, 0);

    % highlight significance
    [H,P] = ttest2(SSE_diffFromGroundTruth(idx,2), SSE_diffFromGroundTruth(idx,3));
    sigstar([1 2], P, false, true)

    % =====================================================================
    %%format the figure
    xTitle      = [];
    yTitle      = [];
    mainTitle   = [];
    fontSize    = 16;
    [hXTitle,hYTitle,hTitle] = fig_figFormat(hFig, xTitle,yTitle,mainTitle, fontSize, false);

    % save
    fig_save(hFig, 'sim4_results', EXPORT_DIR, EXPORT_FORMAT);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Finish Up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % save workspace, in case we want to look at the data posthoc
    fn = sprintf('run_sims_v5_workspace_%s', datestr(now(), 30)); %#ok
    save(fn);

    toc() % report time elapsed£