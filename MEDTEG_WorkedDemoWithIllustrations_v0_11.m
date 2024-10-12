% MEDTEG (Minimum Entropy Dynamic TEst Grid)
%
%   Demo / walk-through of an algorithm that dynamically adds points to a
%   perimetric test grid (e.g., to provide a personalized or more
%   comprehensive Visual Field assessment).
%
%   Note that this code is intended to demo / illustrate its logic. It is
%   designed to be stepped through once from start to finish, and many of
%   the lines of code pertain purely to the generation of graphical
%   illustrations (including illustrations used in the accompanying
%   peer-reviewed paper; see below for reference).
%   
%   For people wanting to *use* MEDTEG in a project, see the separate
%   MEDTEG.m fle for an OOP implementation. That object does not include
%   the code used here to plot/illustrate the functionality, and may be
%   more easy to follow for technically minded readers. This MEDTEG class
%   also includes further examples of simulated usage.
%
% To cite MEDTEG and/or this code, and for further info on the principles
% underlying the algorithm, see:
%	Jones, P.R. (in prep). <<<TITLE>>>
%
% For the latest version of this code, or to contribute improvements:
%   <<<<GITHUB>>>>
%
% @Requires:
%   Files from MATLAB FileExchange:
%       textLoc
%       supLabel
%       inhull
%       BoundedVoronoiDiagram
%   Custom files:
%       QuestPlus.m
%
% @History:
%   v0.0    PRJ     Oct 2020        Initial prototype (C#)
%   v0.1    PRJ     06/08/2024      Initial prototype (MATLAB)
%   v0.2    PRJ     07/08/2024      Full working version using H maximization (no QUEST+)
%   v0.3    PRJ     08/08/2024      Full working version using Quest+ and EH minimization
%   v0.4    PRJ     09/08/2024      Clean up code/documentation. Add pre-specified simulated observers
%   v0.5    PRJ     19/08/2024      Modified to choose location based on the biggest reduction in entropy (i.e., EDH maximization). Plus general clean up and improvements.
%   v0.6    PRJ     20/08/2024      Misc performance improvements, bug fixes, and code tidying
%   v0.7    PRJ     20/08/2024      Further cleaning & hardening
%   v0.8    PRJ     02/09/2024      Aesthetic tweaks for manuscript
%   v0.9    PRJ     06/09/2024      Modify sibson weights so as to respect the horizontal meridian
%   v0.10   PRJ     09/09/2024      Further tweaks to code/comments for consistency with MEDTEG.m
%   v0.11   PRJ     22/09/2024      Improvements to formatting & comments
%
% This code is licenced under CC BY-SA 3.0. This permits the code to be
% used for free, even in commercial applications, as long as: (1) the
% present work is cited, and (2) the derivative work is licensed under the
% same terms.
%
% Copyright (2024) P R Jones <petejonze@gmail.com>
% *********************************************************************
%

% start by clearing everything in memory
clear all %#ok

% user defined params (try changing these)
mode_auto                   = false;            % false to manually step through by pressing SPACE key
do_plots                    = true;             % false to suppress the generation of graphical illustrations
nPointsToAddToGrid          = 4;                % will run MEDTEG for this many cycles (each time adding one new point to the grid)
grid_type                   = 'seeds';          % 'microperimetry'  '24-2', 'seeds'                                    -- as used in getAllTestLocations()
weight_type                 = 'prefer_central'; % 'macula_only', 'prefer_central', 'flat'                   -- as used in getWeighting()
sim_observer_type           = 'hemianopia';     % 'random', 'healthy', 'small_macular_scotoma', 'hemianopia' -- as used in getSimulatedObserverDLS()
excludeAroundMidlines       = true;             % the MEDTEG.m object also includes a further option to exclude candidates in-or-around the blindspot (and will also optional ignore PMF data from points within the physiologic blindspot region). However, we've left all this logic out of this present script for brevity/simplicity.
excludeOutsideConvexHull    = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  !!!!!!!!!!!!! No need to edit anything below this point !!!!!!!!!!!!!  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% finish initialising the workspace: clear/close everything
close all
clc
warning('on');
set(0,'defaulttextinterpreter','tex');
%#ok<*UNRCH>    % disable all code analyser warnings for unreachable lines
%#ok<*SAGROW>   % disable all code analyser warnings for growing arrays

% delete any pre-existing likelihoods cache file (in case params have changed
% since last run simulation)
likelihoodsFn = 'likelihoods.mat';
if exist(likelihoodsFn,'file')
    delete(likelihoodsFn);
end

% start new plot
if do_plots
    hFigCandidates = figure('Position', [800 100 512 512], 'Name','Candidates', 'Color','w');
    hFigPMF = figure('Position',[50 50 512 640], 'Name','PMF', 'Color','w');
end

% #1 get all test points (not excluding the blindspot)
[x, y] = getAllTestLocations(grid_type, false); 

% invent some PMFs for each point (normally these would be determined by
% user responses into QUEST+, but for present purposes we'll just invent
% PMF values based on specified normative data)
[DLS_mu, DLS_sd] = getSimulatedObserverSensitivities(x, y, sim_observer_type);
questPlusObjs = createQuestPlusObjs(length(x), [], DLS_mu, DLS_sd);

for newPointN = 1:nPointsToAddToGrid
    fprintf('Adding point %i', newPointN);

    % start plot by clearing everything and plotting the raw points
    if do_plots
        figure(hFigCandidates);
        clf();
        hold on
        plot(x, y, 'o');
        axis square; 
        box on;
        set(gca, 'xlim',[-30 30], 'ylim',[-30 30]);
        xlabel(['X' char(176)], 'FontSize',13, 'FontWeight','bold');
        ylabel(['Y' char(176)], 'FontSize',13, 'FontWeight','bold');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Stage A: Get Candidate Locations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 1.1 Compute convex hull (everything outside of the starting grid
    % shall be ignored) 
    [k,av] = convhull([x y]);
    chverts = [x(k) y(k)]; % [x y] vertices of the convexhull
    chverts = chverts*1.1; % *1.0 to clip strictly to edges
    if do_plots
        plot(chverts(:,1), chverts(:,2), 'r-'); % plot
    end

    % 1.2 Perform delauney triangulation
    dt = delaunayTriangulation(x, y);
    if do_plots
        triplot(dt, 'k:'); % plot
    end

    % 1.3 Compute the Voronoi diagram.
    V = voronoiDiagram(dt);
    if do_plots
        h = voronoi(dt); delete(h(1)); % plot (but only show the lines, and remove the dots)
    end
    %
    % by convention the first vertex is at infinity. Remove this one
    V(1,:) = [];
    % remove any (near) duplicate vertices (rounding to 4 decimal places)
    V = round(V, 4);
    V = unique(V,'row','stable');

    % 1.4 Get list of candidate locations, where the possible candidates
    % are all the verticies of the voronoi cells 
    candidates = V;

    % 1.5 Get candidate weights (0 to 1; inclusive). These weights will be
    % used in Step 3 to give preferential weighting to the various
    % candidate locations. For now, however, we'll just exclude altogether
    % any points with a weight of 0 (i.e., for computational expediency,
    % since these points by definition will never be selected)
    convHullVerts = []; %#ok
    if (excludeOutsideConvexHull)
        convHullVerts = chverts;
    end
    candWeights = getWeightings(candidates, weight_type, convHullVerts, excludeAroundMidlines);
    candidates(candWeights==0,:) = []; % remove
    candWeights(candWeights==0) = [];

    % 1.6 Check that any candidates available. Abort if not
    if isempty(candidates)
        error('No candidates found. Weightings too aggressive?');
    end

    % 1.7 Plot and label candidates
    nCandidates = size(candidates,1);
    vlabels = arrayfun(@(n) {sprintf('%d', n)}, 1:nCandidates);
    if do_plots
        hTxt = text(candidates(:,1), candidates(:,2)+.2, vlabels, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'BackgroundColor', 'white');
        title('(click or press a key to continue)');
    end

    if ~mode_auto
        waitforbuttonpress;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Stage B: Initialise each point and determine best candidate (the
    % one that maximises the reduction in entropy/mm2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % For each candidate location, compute PMF and Area. For PMF we
    % will compute the weighted average of all the PMFs from the new points
    % Natural Neighbours. For Area we will compute what the new voronoi
    % cell would be if we recomputed the voronoi diagram with this
    % additional location included.    
    [~,V0] = polybnd_voronoi([x y], chverts); % (re)compute cells

    % get Probability Mass Function (PMF) for all points in the
    % original/pre-existing voronoi diagram. Only compute this once here,
    % as it will be constant for all candidate points.
    PMFs = tabulatePMFs(questPlusObjs);
    % [DLSs,sigmas] = tabulateDLSestimates(questPlusObjs)
    % H = -nansum(PMFs .* log2(PMFs), 2); %#ok

    % Find best candidate...
    wEDHdeg2Max = -Inf;
    bestCandidatePoint = [NaN NaN];
    bestCandidateQPlus = NaN;
    for i = 1:nCandidates
        fprintf('.');

        % get the next candidate point
        candidatePoint = candidates(i,:);

        % add candidate to existing points
        xy1 = [[x y]; candidatePoint];
        idx = length(xy1);

        % use polybnd_voronoi as we want to bound the peripheral voronoi
        % cells to the convex hull, rather than simply allowing them to
        % extend out to infinity (in C# this is done slightly differently,
        % by calculating the intersection between the shape object and the
        % convex hull, but is the same basic principle)
        [Neighbours,V1] = polybnd_voronoi(xy1, chverts); % neighbours, verticies

        % get the candidate cell
        candidateCellVertices = V1{idx};

        % get the neighbour points and associated voronoi cells
        nidx = Neighbours{idx};
        neighbourPoints = [x(nidx) y(nidx)];
        neighbourCellsVertices = V0(nidx);

        % plot
        if do_plots
            figure(hFigCandidates);
            % hPlot = plot(candidateCellVertices(:,1), candidateCellVertices(:,2), 'b-', 'linewidth',4);
            hPlot = plot(polyshape(candidateCellVertices), 'facecolor','b','linewidth',1);
            for j = 1:length(nidx)
                hPlot(end+1) = plot([candidatePoint(1) neighbourPoints(j,1)], [candidatePoint(2) neighbourPoints(j,2)], 'c:', 'linewidth',3);
            end
            uistack(hTxt(i),'top'); % ensure candidate point visible
            %
            figure(hFigPMF);
            clf();
        end

        % for each Natural Neighbour of the candidate point, get weights
        % using sibson weighting method. Also get the area of the candidate
        % cell at same time (since we get this for free, as Sibson
        % weighting requires us to compute the area of the overlapping
        % region between the new and existing cells)  
        [w, A] = getSibsonWeights(candidateCellVertices, neighbourCellsVertices, candidatePoint, neighbourPoints, do_plots);

        % for each neighbour, get PMF, then compute the weighted sum of these
        neighboursPMF = PMFs(nidx,:);
        candidatePMF = sum(bsxfun(@times, w, neighboursPMF),1);

        % create a new Q+ object with candidatePMF as the prior
        objs = createQuestPlusObjs(1, candidatePMF);
        candQPlus = objs{1};

        % get the current negative Entropy, H. The closer to 0, the more CERTAIN we are (i.e., the narrower the PMF). The bigger the value the more UNCERTAIN we are
        H = candQPlus.getEntropy();

        % get the expected entropy after next trial, EH
        [~,~,EH] = candQPlus.getTargetStim();

        % get absolute change in entropy (bigger is better)
        EDH = abs(H - EH); % abs is not technically needed since entropy should always go down after next trial. But included here anyway for clarity and defensive completeness
        % fprintf('%1.30f\n%1.30f\n%1.30f', H, EH, H-EH)

        % multiply by area to get change in etropy volume (bigger is better)
        EDHdeg2 = EDH * A;
        % i.e., will  be greatest (best) if the area is large and the expected information gain (change in entropy) is large
        % Std of FoS curve      Area    H           EH          EDH         EDHdeg2
        % 1                     10      3.1712      1.8056      1.3656      13.6565
        %                       40                                          54.6260
        % 10                    10      3.1712      2.1834      0.9878      9.8780
        %                       40                                          39.5122

        % apply location weighting coefficient
        wEDHdeg2 = EDHdeg2 * candWeights(i); % bigger wEDHdeg2 is better. Values with 0 weight will equal 0

        % plot
        if do_plots
            plotPMFs([neighboursPMF; candidatePMF]);
        end

        % store the point that will lead to the greatest gain in
        % information (i.e., greatest decrease in entropy)
        if (wEDHdeg2 > wEDHdeg2Max)
            wEDHdeg2Max = wEDHdeg2;
            bestCandidatePoint = candidatePoint;
            bestCandidateQPlus = candQPlus;
        end

        % return to main window. Add labels
        if do_plots
            figure(hFigCandidates);
            if (bestCandidatePoint == candidatePoint)
                color = 'r';
            else
                color = 'b';
            end
            txt = sprintf('Candidate %i:  Area=%1.1f; H=%1.1f; EDH=%1.2f;  w=%1.1f;  wEDHdeg2=%1.2f', i, A, H, EDH, candWeights(i), wEDHdeg2);
            title({txt, '(click or press a key to continue)'}, 'color',color);
        end

        % wait for button press
        if ~mode_auto
            waitforbuttonpress;
        end
        if do_plots
            delete(hPlot);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Done. Store/report the best candidate location
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % if no valid candidate point could be found, just pick a random
    % location
    if any(isnan(bestCandidatePoint))
        error('no valid point found!');
    end

    % plot final point
    if do_plots
        plot(bestCandidatePoint(1), bestCandidatePoint(2), 'gs', 'MarkerSize',20, 'LineWidth',4);
        txt = sprintf('x = %1.1f; y = %1.1f; wEDHdeg2 = %1.2f', bestCandidatePoint(1), bestCandidatePoint(2), wEDHdeg2Max);
        textLoc(txt, 'SouthEast', 'color','g');
    end

    % store
    x(end+1) = bestCandidatePoint(1);
    y(end+1) = bestCandidatePoint(2);
    questPlusObjs{end+1} = bestCandidateQPlus;

    if ~mode_auto
        waitforbuttonpress;
    end

    fprintf('\n');
end


%% plot final results -----------------------------------------------------
hFigResults = figure('Position',[475 350 560 420], 'Name','Results', 'Color','w');

% plot grid
subplot(1,2,1); hold on
hDat = plot(x, y, 'o');
[xOriginal, yOriginal] = getAllTestLocations(grid_type, false);
hDat(end+1) = plot(xOriginal, yOriginal, 'o', 'markerfacecolor','b');  % overplot original points
axis equal; axis square; 
plot(xlim(),[0 0],'k-');
plot([0 0], ylim(), 'k-')
hLeg = legend(fliplr(hDat), 'Original', 'New Points', 'Location','North','Orientation','horizontal');
hLeg.Position = hLeg.Position.*[1 1.15 1 1]; % nudge upwards

% plot sensitivity estimates
subplot(1,2,2); hold on
DLS = tabulateDLSestimates(questPlusObjs);
txt = strsplit(sprintf('%i,',DLS), ','); txt(end) = [];
%
idx = ismember(x,xOriginal) & ismember(y,yOriginal);
txt(~idx) = {'?'};
%
hDat = plot(x, y, '.w');
axis equal; axis square; 
plot(xlim(),[0 0],'k-');
plot([0 0], ylim(), 'k-')
text(x, y, txt, 'horizontalalignment','center', 'fontweight','bold', 'BackgroundColor','w');
hLeg = legend(fliplr(hDat), 'DLS, dB', 'Location','North','Orientation','horizontal');
hLeg.Position = hLeg.Position.*[1 1.15 1 1]; % nudge upwards

% label axes
hTxt = suplabel(['X' char(176)],'x',[0.25, 0.29, .5 .5]);
hTxt(end+1) = suplabel(['Y' char(176)],'y');
set(hTxt, 'FontSize',13, 'FontWeight','bold');


%% Internal Helper Functions ----------------------------------------------
function [x, y] = getAllTestLocations(grid_type, excludeBlindspot)
    switch grid_type
        case '24-2'
            id = [ % 24-2, right eye format (OD)
                NaN, NaN, NaN,   1,   2,   3,   4, NaN, NaN, NaN % +21
                NaN, NaN,   5,   6,   7,   8,   9,  10, NaN, NaN % +15
                NaN,  11,  12,  13,  14,  15,  16,  17,  18, NaN % +9
                19,   20,  21,  22,  23,  24,  25,  26,  27, NaN % +3
                28,   29,  30,  31,  32,  33,  34,  35,  36, NaN % -3
                NaN,  37,  38,  39,  40,  41,  42,  43,  44, NaN % -9
                NaN, NaN,  45,  46,  47,  48,  49,  50, NaN, NaN % -15
                NaN, NaN, NaN,  51,  52,  53,  54, NaN, NaN, NaN % -21
            ];  %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
            x_deg = ones(8,1) * (-27:6:27);
            y_deg = (21:-6:-21)' * ones(1,10);
    
            % exclude the blindspot
            if (excludeBlindspot)
                idx = ~ismember(id, [26 35]) & ~isnan(id);
            else
                idx = ~isnan(id); % exclude null points
            end
    
            % return vals
            x = x_deg(idx);
            y = y_deg(idx);
            return;
        case 'microperimetry'
            theta = pi * linspace(0,2,9); % [0 pi/4 pi/2 pi]
            theta(end) = [];
            rho = [2 4 6 8 10];
    
            x = []; y = [];
            for t = theta
                for r = rho
                    [x(end+1),y(end+1)] = pol2cart(t,r); %#ok
                end
            end
    
            % add centre
            x(end+1) = 0;
            y(end+1) = 0;
    
            % return
            x = x(:);
            y = y(:);
            return;
        case 'seeds'
            id = [ % right eye format (OD)
                NaN, NaN, NaN,  1 , NaN, NaN,  1 , NaN, NaN, NaN
                NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN
                NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN,  1 , NaN
                 1 , NaN, NaN, NaN,  1 , NaN, NaN, NaN, NaN, NaN
                 1 , NaN, NaN, NaN, NaN,  1 , NaN, NaN, NaN, NaN
                NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN,  1 , NaN
                NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN
                NaN, NaN, NaN,  1 , NaN, NaN,  1 , NaN, NaN, NaN];      
            x_deg = ones(8,1) * (-27:6:27);
            y_deg = (21:-6:-21)' * ones(1,10);
    
            % exclude null points
            idx = ~isnan(id);
    
            % return vals
            x = x_deg(idx);
            y = y_deg(idx);
            return;
        otherwise
            error('getAllTestLocations(): Unknown input: %s ?????', grid_type);
    end
end

function [DLS_mu, DLS_sd] = getSimulatedObserverSensitivities(x_degs, y_degs, sim_observer_type)

    % #1 define "true" sensitivity at various predefined locations -- shall
    % assume Right Eye (OD)
    switch sim_observer_type
        case 'random'
            DLS_mu = randi([0 35]); % randomly pick a value between 0 and 35
            DLS_sd = 3;
            return; % in this special case can just return a random number and move on
        case 'healthy' % these values are from Brenton & Phelps (1986), for young (18 years) normally sighted adults
            % typical DLS values for each location (i.e., the mu parameter
            % of the psychometric function):
            norm_DLS = [ % right eye format (OD)
                    NaN  NaN  NaN  26.8 26.3 25.5 25.1 NaN  NaN  NaN    % +27
                    NaN  NaN  29.0 29.2 28.3 27.6 27.8 28.0 NaN  NaN    % +21
                    NaN  28.9 29.1 29.8 30.2 29.8 30.9 30.4 28.8 NaN    % +15
                    28.3 29.6 31.2 32.6 32.3 31.7 31.2 30.6 30.4 30.3   % +9
                    28.9 30.4 32.1 32.8 33.7 33.5 32.2 NaN  31.3 31.9   % +3
                    29.4 31.0 32.6 33.3 33.8 33.8 32.7 NaN  31.5 31.3   % -3
                    29.6 30.4 31.5 33.4 32.9 32.8 32.8 31.1 31.6 31.1   % -9
                    NaN  30.5 30.6 31.2 32.1 31.3 30.9 32.4 31.2 NaN    % -15
                    NaN  NaN  30.1 29.8 30.9 31.4 31.8 31.2 NaN  NaN    % -21
                    NaN  NaN  NaN  29.0 30.1 30.8 30.7 NaN  NaN  NaN    % -27
                ]; %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
            % typical measurement variability for each location:
            norm_SD = [
                    NaN  NaN  NaN  5.8  3.5  4.7  3.9  NaN  NaN  NaN    % +27
                    NaN  NaN  3.6  2.6  4.1  2.6  2.6  2.0  NaN  NaN    % +21
                    NaN  2.9  3.2  2.7  2.2  2.4  2.2  2.3  2.7  NaN    % +15
                    2.7  2.4  1.6  1.7  1.6  2.1  2.1  2.7  2.8  3.2    % +9
                    2.5  2.2  1.7  1.5  2.0  1.7  1.9  NaN  2.9  2.4    % +3
                    2.3  1.7  1.4  2.2  1.4  1.4  2.1  NaN  2.5  2.0    % -3
                    2.2  2.1  2.5  1.7  1.9  1.8  2.2  3.2  2.0  2.6    % -9
                    NaN  1.8  1.8  2.8  2.1  2.0  2.3  1.9  2.4  NaN    % -15
                    NaN  NaN  2.5  2.4  2.1  1.9  2.1  2.9  NaN  NaN    % -21
                    NaN  NaN  NaN  3.0  2.2  2.6  2.4  NaN  NaN  NaN    % -27
                ]; %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
            norm_x_degs = ones(10,1) * (-27:6:27);
            norm_y_degs = (27:-6:-27)' * ones(1,10);
        case 'small_macular_scotoma'
            error("implement me");
        case 'hemianopia'
            norm_DLS = [ % right eye format (OD)
                     0    0    0    0    0   25.5 25.1 NaN  NaN  NaN    % +27
                     0    0    0    0    0   27.6 27.8 28.0 NaN  NaN    % +21
                     0    0    0    0    0   29.8 30.9 30.4 28.8 NaN    % +15
                     0    0    0    0    0   31.7 31.2 30.6 30.4 30.3   % +9
                     0    0    0    0    0   33.5 32.2 NaN  31.3 31.9   % +3
                     0    0    0    0    0   33.8 32.7 NaN  31.5 31.3   % -3
                     0    0    0    0    0   32.8 32.8 31.1 31.6 31.1   % -9
                     0    0    0    0    0   31.3 30.9 32.4 31.2 NaN    % -15
                     0    0    0    0    0   31.4 31.8 31.2 NaN  NaN    % -21
                     0    0    0    0    0   30.8 30.7 NaN  NaN  NaN    % -27
                ]; %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
            % let's assume that we're pretty certain about the level of
            % loss in the damaged hemifield:
            norm_SD = [
                    0.1  0.1  0.1  0.1  0.1  4.7  3.9  NaN  NaN  NaN    % +27
                    0.1  0.1  0.1  0.1  0.1  2.6  2.6  2.0  NaN  NaN    % +21
                    0.1  0.1  0.1  0.1  0.1  2.4  2.2  2.3  2.7  NaN    % +15
                    0.1  0.1  0.1  0.1  0.1  2.1  2.1  2.7  2.8  3.2    % +9
                    0.1  0.1  0.1  0.1  0.1  1.7  1.9  NaN  2.9  2.4    % +3
                    0.1  0.1  0.1  0.1  0.1  1.4  2.1  NaN  2.5  2.0    % -3
                    0.1  0.1  0.1  0.1  0.1  1.8  2.2  3.2  2.0  2.6    % -9
                    0.1  0.1  0.1  0.1  0.1  2.0  2.3  1.9  2.4  NaN    % -15
                    0.1  0.1  0.1  0.1  0.1  1.9  2.1  2.9  NaN  NaN    % -21
                    0.1  0.1  0.1  0.1  0.1  2.6  2.4  NaN  NaN  NaN    % -27
                ]; %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
            % alternatively, if we want to maintain the "default" level of
            % measurement certainty:
            % norm_SD = [
            %         NaN  NaN  NaN  5.8  3.5  4.7  3.9  NaN  NaN  NaN    % +27
            %         NaN  NaN  3.6  2.6  4.1  2.6  2.6  2.0  NaN  NaN    % +21
            %         NaN  2.9  3.2  2.7  2.2  2.4  2.2  2.3  2.7  NaN    % +15
            %         2.7  2.4  1.6  1.7  1.6  2.1  2.1  2.7  2.8  3.2    % +9
            %         2.5  2.2  1.7  1.5  2.0  1.7  1.9  NaN  2.9  2.4    % +3
            %         2.3  1.7  1.4  2.2  1.4  1.4  2.1  NaN  2.5  2.0    % -3
            %         2.2  2.1  2.5  1.7  1.9  1.8  2.2  3.2  2.0  2.6    % -9
            %         NaN  1.8  1.8  2.8  2.1  2.0  2.3  1.9  2.4  NaN    % -15
            %         NaN  NaN  2.5  2.4  2.1  1.9  2.1  2.9  NaN  NaN    % -21
            %         NaN  NaN  NaN  3.0  2.2  2.6  2.4  NaN  NaN  NaN    % -27
            %     ]; %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27

            norm_x_degs = ones(10,1) * (-27:6:27);
            norm_y_degs = (27:-6:-27)' * ones(1,10);
        otherwise
            error('getSimulatedObserverDLS(): Unknown input: %s ?????', sim_observer_type);
    end

    % #2 use interpolatation to get sensitivity at requested location (NB:
    % cleverer forms of interpolation exist, but for present purposes we
    % just need something simple)
    idx = ~isnan(norm_DLS);
    mapFunc_mu = scatteredInterpolant(norm_x_degs(idx), norm_y_degs(idx), norm_DLS(idx)); % in a real setting this would likely only be computed once, at object initialisation. But written this way to minimize lines of code.
    mapFunc_sd = scatteredInterpolant(norm_x_degs(idx), norm_y_degs(idx), norm_SD(idx));

    % compute
    DLS_mu = mapFunc_mu(x_degs, y_degs);
    DLS_sd = mapFunc_sd(x_degs, y_degs);
    
    % validate
    if (any(isnan(DLS_mu)) || any(isnan(DLS_sd)))
        fprintf('%1.2f %1.2f %1.2f\n', [norm_x_degs(:) norm_y_degs(:) DLS_mu(:) DLS_sd(:)]')
        error('one or more XY values produced NaN (see above). Provide normative data for this location, or restrict your test grid')
    end
end

% Natural neighbor (or "Sibson") interpolation is a method of spatial interpolation, developed by Robin Sibson. The method is based on Voronoi tessellation of a discrete set of spatial points. This has advantages over simpler methods of interpolation, such as nearest-neighbor interpolation, in that it provides a smoother approximation to the underlying "true" function.
function [weights, candidateCellArea] = getSibsonWeights(candidateCellVertices, neighbourCellVertices, candidateXY, neighbourXYs, do_plots)
    % get the voronoi cell of the proposed new point
    P0 = polyshape(candidateCellVertices);
    candidateCellArea = area(P0);

    % plot
    if do_plots
        subplot(8,2,1:8);
        cla();
        axis square;
        hold on
        plot(P0);
        plot(candidateXY(1), candidateXY(2),'ko', 'MarkerFaceColor','k');
        set(gca, 'XTick',[], 'YTick',[]);
        hTxt = nan(1, length(neighbourCellVertices));
    end

    % compute how "close" this new cell is to each of its neighbours, where closeness is determined by how much this new cell overlaps with each of its neighbours in the previous/original voronoi graph
    n = length(neighbourCellVertices);
    weights = nan(n,1);
    for i = 1:n
        % get polygon defining the neighbouring cell (from the pre-existing
        % point)
        V = round(neighbourCellVertices{i}, 4);
        V = unique(V,'row','stable'); % remove any (near) duplicate vertices (rounding to 4 decimal places)
        P1 = polyshape(V);

        % compute the area of the intersecting polygon (i.e., the degree of overlap between the Voronoi cell of the new point and this neighbour)
        I = intersect(P0, P1);
    
        % compute Sibson weight (proportion of the new cell taken up by the intersection with this neighbour)
        weights(i) = area(I) / candidateCellArea;

        % plot
        if do_plots
            box on;
            plot(P1);
            plot(I); % shade intersection
            txt = sprintf('w_%i = %1.3f', i, weights(i));
            hTxt(i) = text(neighbourXYs(i,1), neighbourXYs(i,2), txt, 'interpreter','tex');
        end
    end

    % finally, manually override any neighbours that lie on the opposing
    % side of the horizontal meridian (i.e., arrising from different
    % retinal nerve fibres), and renormalize the weights to sum to 1
    idx = sign(neighbourXYs(:,2)) ~= sign(candidateXY(2)); % find any points not in same vertical hemifield
    weights(idx) = 0;
    weights = weights .* 1./sum(weights);
    % update reported weights (if required)
    if any(idx) && do_plots
        for i = find(idx)'
            set(hTxt(i), 'String', sprintf('w_%i = 0', i));
        end
    end
end

function sigma = getPsySlope(DLS)
    % NB: By making the slope of the likelihood function
    % larger for some points, the algorithm will be disinclined to test
    % these points. In practice, the slope might be set to vary with DLS or
    % retinal location, and/or be estimated by QUEST+, or for simplicity be
    % fixed at a constant value.
    %
    % Here we shall set it to be proportional to DLS (i.e., the mu of the
    % psychometric function AKA frequency of seeing curve AKA likelihood
    % function). NB: this function is entirely made-up, but similar real
    % functions can be found in the literature.
    sigma = 1 + 1.1^-DLS * 5; % e.g., DLS=0, slope=6.0000    DLS=10, slope=2.9277    DLS=20, slope=1.7432    DLS=30, slope=1.2865    DLS=40, slope=1.1105

    % ALT1: constant:
    % sigma = 2;

    % ALT2: vary with eccentricity:
    % (not implemented)
end

function objs = createQuestPlusObjs(n, manualPrior, true_DLS_mu, true_DLS_sd)
    % set model
    lambda = 0.05; % lapse rate. Hardcoded. Could (should!) be specified as part of the QUEST model (in which case could easily be made a free parameter -- see QuestPlus.m for examples), but have done it this way for simplicity
    gamma = 0.02; % guess rate. See above.
    %F = @(x,mu)(gamma + (1-gamma-lambda) * (1-normcdf(x,mu,2)  )); % simple, fixed slope
    F = @(x,mu)(gamma + (1-gamma-lambda) * (1-normcdf(x,mu,getPsySlope(mu))  )); % only specified probability of responding correctly, p(c). QuestPlus will automatically assume that the probability of response correct is 1 - p(c).
    %x = 0:40; y = F(0:40, 20); close all; plot(x,y,'-o'); ylim([0 1])

    % set common QUEST+ object params
    stimDomain      = linspace(0, 40, 41);
    paramDomain     = linspace(0, 40, 41); % the possible DLS_mu values
    respDomain    	= [0 1];
    stopRule       	= 'entropy';
    stopCriterion  	= 3;
    minNTrials    	= 2;
    maxNTrials   	= 5;

    % try to load in likelihoods to prevent having to recompute
    likelihoodsDat = [];
    likelihoodsFn = 'likelihoods.mat';
    if exist(likelihoodsFn,'file')
        likelihoodsDat = load(likelihoodsFn);
    end

    % create each
    objs = cell(1,n);
    for i = 1:n
        % create obj
        objs{i} = QuestPlus(F, stimDomain, paramDomain, respDomain, stopRule, stopCriterion, minNTrials, maxNTrials);

        % construct new (random) prior(s)
        if nargin < 2 || isempty(manualPrior)
            priors = normpdf(paramDomain, true_DLS_mu(i), true_DLS_sd(i));
        else
            priors = manualPrior;
        end
        priors = priors./sum(priors);

        % initialise priors/likelihoods
        objs{i}.initialise(priors, likelihoodsDat);

        % if not done so already, save/cache these values to speed up subsequent initialisations
        if isempty(likelihoodsDat)
            objs{i}.saveLikelihoods('likelihoods.mat');
            likelihoodsDat = load(likelihoodsFn);
        end
    end
end
                    
function PMFs = tabulatePMFs(questPlusObjs)
    tmp = cellfun(@(obj)obj.posterior, questPlusObjs, 'UniformOutput', false);
    PMFs = cat(1, tmp{:}); % n x m matrix, where n is the number of questplus objects, and m is the length of the posterior distribution vector
end

function DLS = tabulateDLSestimates(questPlusObjs)
    tmp = cellfun(@(obj)obj.getParamEsts('mode'), questPlusObjs, 'UniformOutput', false);
    DLS = [tmp{:}]'; % n x 1 vector, where n is the number of questplus objects
end

function plotPMFs(PMFs)
    n = size(PMFs,1);
    x = 0:40;

    % get handles to polygons so that we can match the color
    hPolygons = findobj('Type', 'polygon');

    % defensive
    if n > 8
        warning('modify n subplots if wanting to plot more than 7 neighbours');
        return;
    end

    % plot each PMF
    for i = 1:n
        if i==n
            %label = [char(9679) '  candidate'];
            label = {'', '$\circ$ \textbf{candidate}', '$\sum_i w_i p_i(x)$'};
            interpreter = 'latex';
            color = 'k'; % candidate color
        else
            label = sprintf('i = %i',i);
            interpreter = 'tex';
            color = hPolygons(end-i*2).FaceColor;
        end

        % plot
        subplot(8, 2, 8+i); 
        plot(x, PMFs(i,:), '-', 'color',color, 'linewidth',3);

        % format
        textLoc(label, 'North', 'color',color, 'interpreter',interpreter);
        set(gca, 'YLim',[0 1], 'YTick',[0 1]);
        if i < (n-1)
            set(gca,'XTickLabel',[]);
        end
    end

    % label axes
    hTxt = suplabel('DLS, dB');
    hTxt(end+1) = suplabel('PMF','y',[.13 .075 .85 .6]);
    set(hTxt, 'FontSize',13, 'FontWeight','bold');

    % done
    refresh();
end

function w = getWeightings(xy, weight_type, convHullVerts, excludeAroundMidlines)
    % 0. init column vector
    w = ones(size(xy,1),1);

    % 1. Optionally exclude points that fall outside of the convex hull
    % (this may have been done already by voronoiDiagram(), but we'll
    % include this step anyway, just as an explicit sanity check)
    if ~isempty(convHullVerts)
        idx = ~inhull(xy, convHullVerts, [], 0.1);
        w(idx) = 0;
    end

    % 2. Optionally prevent any points within 1 degree of the horizontal meridian:
    if (excludeAroundMidlines)
        %idx = any(abs(xy)<1,2); % to also exclude vertical meridian
        idx = abs(xy(:,2)) < 1; % exclude horizontal meridian
        w(idx) = 0;
    end

    % 3. then also enforce any other rules that the user specifies:
    switch weight_type
        case 'flat'
            w = w*1;
        case 'macula_only'
            w = w .* (sqrt(sum(xy.^2,2)) < 10); % only allow points within central 10 degrees
        case 'prefer_central'
            maxEccentricity_deg = 40; % arbitrary
            w1 = 1 - sqrt(sum(xy.^2,2))/maxEccentricity_deg;
            w1 = min(max(w1, 0), 1);
            w1 = w1 .^ 2; % make even stronger
            w = w .* w1; % give more weight to more central locations
        otherwise
            error('getWeighting(): Unknown input: %s ?????', weight_type);
    end

    % 4. Defensive: clamp to ensure w is between 0 and 1 (inclusive)
    w = min(max(w, 0), 1);
end