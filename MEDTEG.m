classdef MEDTEG < handle
    % Matlab implementation of MEDTEG (Minimum Entropy Dynamic TEst Grid)
    %
    %   >>> How to use:
    %   For examples of use, see MEDTEG.runExample()
    %
    %   >>> Background Info:
    %   TBC
    %
    %   >>> For further info, and to cite:
    %   TBC
    %
    % MEDTEG Methods:
    %   * MEDTEG                - MEDTEG Constructor.
    %   * GetNewTestLocation    - Determine the optimal next/new test location.
    %   * disp                  - Print to console info regarding the internal state of the current MEDTEG object.
    %
    % Public Static Methods:
    %   * runExample	        - Minimal-working-example(s) of usage
    %
    % Examples of use:
    %   MEDTEG.runExample(1)
    %
    % Author:
    %   Pete R Jones <petejonze@gmail.com>
    %
    % Verinfo:
    %   0.0.1	PRJ	06/09/2024 : first_build
    %   0.0.2	PRJ	09/09/2024 : added example usage demos
    %
    % Copyright 2024 : P R Jones <petejonze@gmail.com>
    % *********************************************************************
    %
    % @TODO finalise and test
    % @TODO improve definition of physiologic blindspot

    %% ====================================================================
    %  -----PROPERTIES-----
    %$ ====================================================================      

    properties (GetAccess = public, SetAccess = private)

        % XXXXX
        candWeightingFunc;
        weightFuncParams;

        % object holding QUEST+ parameters (see getExampleQuestPlusParams()
        % for details)
        QPlusParams

        % XXXXXX
        excludeOutsideConvexHull = true; % TODO currently no way of setting in constructor
    end


    properties (GetAccess = private, SetAccess = private)
        chverts; % convex hull vertices. Stored for expedience, as re-used at various points during the code
    end

    
 	%% ====================================================================
    %  -----PUBLIC METHODS-----
    %$ ====================================================================
    
    methods (Access = public)
        
        %% == CONSTRUCTOR =================================================
        
        function obj = MEDTEG(candWeightingFunc, weightFuncParams, QPlusParams)
            % MEDTEG Constructor.
            %
            %   Create a new MEDTEG object.
            %       
            % @param    XXXXXX          XXXX
            %                             E.g.: XXXX
            %                             Default: XXXX 
            % @return   QuestPlus       MEDTEG object handle
            %
            % @date     07/09/24
            % @author   PRJ
            %
            
            % store user specified vals
            %obj.excludeOutsideConvexHull = excludeOutsideConvexHull;
            obj.candWeightingFunc = candWeightingFunc;
            obj.weightFuncParams = weightFuncParams;
            obj.QPlusParams = QPlusParams;
        end
        
        %% == METHODS =================================================
        
        function [bestCandidateXY, bestCandidateQPlus, wEDHdeg2Max, candidates, cand_weights, cand_wEDHdeg2] = GetNewTestLocation(obj, existingTestLocations_XY, existingTestLocations_QuestPlus)
            % XXXXX.
            %
            %   XXXXX.
            %       
            % @param    existingTestLocations_XY        % Nx2 matrix of xy locations for the N existing test-locations
            %                             E.g.: XXXX
            %                             Default: XXXX 
            % @param    existingTestLocations_QuestPlus % Nx1 vector of Quest+ objects for the N existing test-locations (which we shall query to extract their latest posterior density PMF values)
            %                             E.g.: XXXX
            %                             Default: XXXX 
            % @return   XXXXX           XXXXX
            %
            % @date     07/09/24
            % @author   PRJ
            %

            [candidates, cand_weights] = obj.getListOfCandidates(existingTestLocations_XY);
            [bestCandidateXY, bestCandidateQPlus, wEDHdeg2Max, cand_wEDHdeg2, ~] = obj.selectBestCandidate(existingTestLocations_XY, existingTestLocations_QuestPlus, candidates, cand_weights);
        end

        % function [] = disp(obj)
        %     % Print to console info regarding the internal state of the
        %     % current MEDTEG object.
        %     %
        %     %   NB: Uses fprintf to write to current print feed, so the
        %     %   output could in principle be rerouted to an offline log
        %     %   file.
        %     %
        %     % @date     07/09/24
        %     % @author   PRJ
        %     %
        % 
        %     fprintf('   a MEDTEG object. Further info not yet available\n\n');
        % end


    end

 	%% ====================================================================
    %  -----PRIVATE METHODS-----
    %$ ====================================================================
    
    methods (Access = private)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Stage A: Get Candidate Locations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [candidates, candWeights] = getListOfCandidates(obj, existingTestLocations_XY)
            % 1.1 Compute convex hull (everything outside of the starting grid
            % shall be ignored)
            k = convhull(existingTestLocations_XY);
            if (isempty(obj.chverts)) % only do this on first run. After which just use the cached value
                obj.chverts = existingTestLocations_XY(k,:); % [x y] vertices of the convexhull
                obj.chverts = obj.chverts*1.1; % *1.0 to clip strictly to edges
            end

            % 1.2 Perform delauney triangulation
            dt = delaunayTriangulation(existingTestLocations_XY(:,1), existingTestLocations_XY(:,2));

            % 1.3 Compute the Voronoi diagram.
            V = voronoiDiagram(dt);
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
            convHullVerts = [];
            if (obj.excludeOutsideConvexHull)
                convHullVerts = obj.chverts;
            end
            candWeights = obj.candWeightingFunc(candidates, convHullVerts, obj.weightFuncParams);
            candidates(candWeights==0,:) = []; % remove
            candWeights(candWeights==0) = [];

            % 1.6 Check that any candidates available. Throw error if not
            if isempty(candidates)
                error('No candidates found. Weightings too aggressive?');
            end
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Stage B: Initialise each point and determine best candidate (the
        % one that maximises the reduction in entropy/mm2)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [bestCandidateXY, bestCandidateQPlusObj, wEDHdeg2Max, all_wEDHdeg2, all_candQPlusObjs] = selectBestCandidate(obj, existingTestLocations_XY, existingTestLocations_QuestPlus, candidates, candWeights)
            % defensive user-input validation
            if size(existingTestLocations_XY,1) ~= length(existingTestLocations_QuestPlus)
                error('mismatch ????? n=%i xy locations, but %i QUEST+ objects ', size(existingTestLocations_XY,1), length(existingTestLocations_QuestPlus));
            elseif length(candidates) ~= length(candWeights)
                error('mismatch ????? n=%i candidate xy locations, but %i weights', length(candidates), length(candWeights));
            end

            % compute params
            nCandidates = size(candidates, 1);

            % For each candidate location, compute PMF and Area. For PMF we
            % will compute the weighted average of all the PMFs from the new points
            % Natural Neighbours. For Area we will compute what the new voronoi
            % cell would be if we recomputed the voronoi diagram with this
            % additional location included.
            [~,V0] = polybnd_voronoi(existingTestLocations_XY, obj.chverts); % (re)compute cells
    
            % get Probability Mass Function (PMF) for all points in the
            % original/pre-existing voronoi diagram. Only compute this once here,
            % as it will be constant for all candidate points.
            PMFs = MEDTEG.tabulatePMFs(existingTestLocations_QuestPlus);
    
            % Compute H, EDH, EDHdeg2, & wEDHdeg2 for each candidate
            all_wEDHdeg2 = nan(nCandidates, 1);
            all_candQPlusObjs = cell(nCandidates, 1);
            for i = 1:nCandidates
                % get the next candidate point
                candidatePoint = candidates(i,:);
    
                % add candidate to existing points
                xy1 = [existingTestLocations_XY; candidatePoint];
                idx = length(xy1);
    
                % use polybnd_voronoi as we want to bound the peripheral voronoi
                % cells to the convex hull, rather than simply allowing them to
                % extend out to infinity (in C# this is done slightly differently,
                % by calculating the intersection between the shape object and the
                % convex hull, but is the same basic principle)
                % xy1
                [Neighbours,V1] = polybnd_voronoi(xy1, obj.chverts); % neighbours, verticies
    
                % get the candidate cell
                candidateCellVertices = V1{idx};
    
                % get the neighbour points and associated voronoi cells
                nidx = Neighbours{idx};
                neighbourPoints = existingTestLocations_XY(nidx,:);
                neighbourCellsVertices = V0(nidx);

                % for each Natural Neighbour of the candidate point, get weights
                % using sibson weighting method. Also get the area of the candidate
                % cell at same time (since we get this for free, as Sibson
                % weighting requires us to compute the area of the overlapping
                % region between the new and existing cells)
                [w, A] = MEDTEG.getSibsonWeights(candidateCellVertices, neighbourCellsVertices, candidatePoint, neighbourPoints);
                
                % if no valid neighbors best practice is undefined. For now we
                % shall skip this point
                if any(isnan(w))
                    all_wEDHdeg2(i) = -inf;
                    all_candQPlusObjs{i} = MEDTEG.createQuestPlusObj(obj.QPlusParams, ones(1, size(PMFs,2))); % no/flat prior
                    continue
                end

                % for each neighbour, get PMF, then compute the weighted sum of these
                neighboursPMF = PMFs(nidx,:);
                candidatePMF = sum(bsxfun(@times, w, neighboursPMF),1);
    
                % create a new Q+ object with candidatePMF as the prior
                candQPlusObj = MEDTEG.createQuestPlusObj(obj.QPlusParams, candidatePMF);
    
                % get the current negative Entropy, H. The closer to 0, the more CERTAIN we are (i.e., the narrower the PMF). The bigger the value the more UNCERTAIN we are
                H = candQPlusObj.getEntropy();
    
                % get the expected entropy after next trial, EH
                [~,~,EH] = candQPlusObj.getTargetStim();
    
                % get absolute change in entropy (bigger is better)
                EDH = abs(H - EH); % abs is not technically needed since entropy should always go down after next trial. But included here anyway for clarity and defensive completeness
                % fprintf('%1.30f\n%1.30f\n%1.30f', H, EH, H-EH)
    
                % multiply by area to get change in etropy volume (bigger is better)
                EDHdeg2 = EDH * A;
    
                % apply location weighting coefficient
                wEDHdeg2 = EDHdeg2 * candWeights(i); % bigger totalEHvol is better. Values with 0 weight will equal 0
    
                % store
                all_wEDHdeg2(i) = wEDHdeg2;
                all_candQPlusObjs{i} = candQPlusObj;

                % debugging (defensive)
                if isnan(wEDHdeg2)
                    % fprintf('candidate_xy: %.1f, %1.1f\n', candidatePoint);
                    % fprintf('neighbours_xy:\n'); neighbourPoints %#ok
                    % fprintf('neighbours_w:\n'); w %#ok
                    % fprintf('neighbours_PMF:\n'); neighboursPMF %#ok
                    % 
                    % for j = 1:length(existingTestLocations_QuestPlus)
                    %     xy = existingTestLocations_XY(j,:);
                    %     dB = existingTestLocations_QuestPlus{j}.getParamEsts('mode');
                    %     fprintf('<x=%i, y=%i> %1.1f dB\n', xy, dB);
                    % end
                    error('wEDHdeg2 is NaN????\n  PMF = [%s]\n   H = %1.1f\n  EH = %1.1f\n  EDH = %1.1f\n  EDHdeg2 = %1.1f\n  wEDHdeg2 = %1.1f\n', num2str(candidatePMF), H, EH, EDH, EDHdeg2, wEDHdeg2);
                end
            end
    
            % Select best candidate (i.e., the point that will lead to the
            % greatest gain in information (i.e., greatest decrease in
            % entropy). In future this could be updated to allow for different
            % selection strategies (e..g, top N, stochastic, etc.). For now
            % we'll also return the raw vals anyway, so user is free to use
            % them in whatever way they chose.
            [wEDHdeg2Max, idx] = max(all_wEDHdeg2);
            bestCandidateXY = candidates(idx,:);
            bestCandidateQPlusObj = all_candQPlusObjs{i};
        end
    
    end



   	%% ====================================================================
    %  -----STATIC METHODS (private) [Internal helper functions] -----
    %$ ====================================================================
      
    methods (Static, Access = private)

        function [weights, candidateCellArea] = getSibsonWeights(candidateCellVertices, neighbourCellVertices, candidateXY, neighbourXYs)
            % get the voronoi cell of the proposed new point
            P0 = polyshape(candidateCellVertices);
            candidateCellArea = area(P0);
    
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
            end
    
            % manually override (ignore) any neighbours that lie on the opposing
            % side of the horizontal meridian (i.e., arrising from different
            % retinal nerve fibres)
            idx = sign(neighbourXYs(:,2)) ~= sign(candidateXY(2)); % find any points not in same vertical hemifield
            weights(idx) = 0;

            % manually override (ignore) any neighbours that fall within
            % the physiologic blindspot (which for simplicity will define as
            % any point within 3 deg of <15,3> OR <15,-3>)
            crit_deg = 3;
            d_upper = sqrt( (neighbourXYs(:,1) - 15).^2 + (neighbourXYs(:,2) - 3).^2 );
            d_lower = sqrt( (neighbourXYs(:,1) - 15).^2 + (neighbourXYs(:,2) - -3).^2 );
            idx = d_upper<crit_deg | d_lower<crit_deg;
            weights(idx) = 0;

            % renormalize the weights to sum to 1 (will return NaN if all
            % neighbors == 0)
            weights = weights .* 1./sum(weights);
        end

        function QP = createQuestPlusObj(QPlusParams, prior)
            % create obj
            QP = QuestPlus(QPlusParams.F, QPlusParams.stimDomain, QPlusParams.paramDomain, QPlusParams.respDomain, QPlusParams.stopRule, QPlusParams.stopCriterion, QPlusParams.minNTrials, QPlusParams.maxNTrials);

            % ensure that the prior sums to 1 (defensive)
            prior = prior./sum(prior);

            % initialise priors/likelihoods
            QP.initialise(prior, QPlusParams.likelihoodsDat);
        end

        function PMFs = tabulatePMFs(questPlusObjs)
            tmp = cellfun(@(QP)QP.posterior, questPlusObjs, 'UniformOutput', false);
            PMFs = cat(1, tmp{:}); % n x m matrix, where n is the number of questplus objects, and m is the length of the posterior distribution vector. Note that for simplicity we are largely assuming that the PMF is 1D, but the code could be easily updated/tweaked to work with N-D PMFs.
        end

    end

    

   	%% ====================================================================
    %  -----STATIC METHODS (public)-----
    %$ ====================================================================
      
    methods (Static, Access = public)

        function weightFuncParams = getDefaultCandWeightingParameters()
            weightFuncParams = struct();
            weightFuncParams.weight_type = 'flat'; % 'macula_only', 'prefer_central', 'flat'                   -- as used in getWeighting()
            weightFuncParams.excludeCandidatesAroundMidlines = true;
            weightFuncParams.excludeCandidatesInOrAroundBlindspot = true;
        end

        function QPlusParams = getExampleQuestPlusParams()            
            % set model
            lambda = 0.05; % lapse rate. Hardcoded. Could (should!) be specified as part of the QUEST model (in which case could easily be made a free parameter -- see QuestPlus.m for examples), but have done it this way for simplicity
            gamma = 0.02; % guess rate. See above.
            QPlusParams.F = @(x,mu)(gamma + (1-gamma-lambda) * (1-normcdf(x,mu,2)  )); % simple, fixed slope
            % %x = 0:40; y = F(0:40, 20); close all; plot(x,y,'-o'); ylim([0 1])

            % set common QUEST+ object params
            QPlusParams.stimDomain      = linspace(0, 40, 41);
            QPlusParams.paramDomain     = linspace(0, 40, 41); % the possible DLS_mu values
            QPlusParams.respDomain    	= [0 1];
            QPlusParams.stopRule       	= 'entropy';
            QPlusParams.stopCriterion  	= 2.6; % once entropy reaches this will stop...
            QPlusParams.minNTrials    	= 2;
            QPlusParams.maxNTrials   	= 6; % .. or once we've reached this N trials (whichever comes first)

            % precompute likelihoods (see QuestPlus.m)
            L = nan(length(QPlusParams.stimDomain), length(QPlusParams.paramDomain), length(QPlusParams.respDomain));
            for j = 1:size(QPlusParams.paramDomain,2)
                L(:,j,2) = QPlusParams.F(QPlusParams.stimDomain, QPlusParams.paramDomain(j)); % P('success')
            end
            L(:,:,1) = 1 - L(:,:,2); % P('failure'): complement of P('success')
            % store (since in this example will be same for all candidate locations. Otherwise leave this variable empty to force QuestPlus to manually recompute every time)
            QPlusParams.likelihoodsDat = struct();
            QPlusParams.likelihoodsDat.stimDomain  = QPlusParams.stimDomain; % included for validation purposes
            QPlusParams.likelihoodsDat.paramDomain = QPlusParams.paramDomain;
            QPlusParams.likelihoodsDat.respDomain  = QPlusParams.respDomain;
            QPlusParams.likelihoodsDat.likelihoods = L;
        end

        function PMF = getExamplePrior(possibleDLSvals, mu, sigma)       
            % define PMF as a mixture of normal and abnormal PDF. Similar
            % in concept to: Vingrys & Pianta, Optom Vis Sci, 1999: "A new look at..."?
            % but here we'll use made up values for convenience & simplicty

            % validation
            if (sigma <= 0)
                error('sigma (%1.2f) must be a positive scalar value', sigma);
            end

            % define glaucomatous distribution
            PMF_glau = normpdf(possibleDLSvals, 10, 5);
            PMF_glau = PMF_glau ./ sum(PMF_glau);

            % define normally sighted distribution (centered on specified
            % mean value, which will vary with eccentricity)
            PMF_ctrl = normpdf(possibleDLSvals, mu, sigma);
            PMF_ctrl = PMF_ctrl ./ sum(PMF_ctrl);

            % define flat distribution (as we want to ensure no prior = 0)
            PMF_flat = ones(1, length(possibleDLSvals));
            PMF_flat = PMF_flat ./ sum(PMF_flat);

            % weights (should sum to 1)
            w_glau = 0.2;
            w_ctrl = 0.6;
            w_flat = 0.1;

            % compute linear weighted sum
            PMF = w_glau*PMF_glau + w_ctrl*PMF_ctrl + w_flat*PMF_flat;

            % defensive (ensure still sum to 1)
            PMF = PMF ./ sum(PMF);

            % plot (for debugging)
            % figure();
            % plot(possibleDLSvals, PMF);

            % defensive validation
            if any(isnan(PMF))
                error('PMF is NaN????\n  mu = %1.1f\n  sigma = %1.1f\n  PMF_glau = [%s]\n  PMF_ctrl = [%s]\n  PMF_flat = [%s]', mu, sigma, num2str(PMF_glau), num2str(PMF_ctrl), num2str(PMF_flat));
            end
        end

        function w = exampleCandWeightingFunc(xy, convHullVerts, params)
            % -1. unpack input parameters structure
            weight_type = params.weight_type;
            excludeCandidatesAroundMidlines = params.excludeCandidatesAroundMidlines;
            excludeCandidatesInOrAroundBlindspot = params.excludeCandidatesInOrAroundBlindspot;

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
            if (excludeCandidatesAroundMidlines)
                idx = any(abs(xy)<1,2); % to also exclude vertical meridian
                % idx = abs(xy(:,2)) < 1; % exclude horizontal meridian
                w(idx) = 0;
            end

            % 3. Optionally prevent any points within 3 degree of the
            % physiologic blindspot (which for simplicity will define as
            % any point within N deg of <15,3> OR <15,-3>)
            if (excludeCandidatesInOrAroundBlindspot)
                crit_deg = 3;
                d_upper = sqrt( (xy(:,1) - 15).^2 + (xy(:,2) - 3).^2 );
                d_lower = sqrt( (xy(:,1) - 15).^2 + (xy(:,2) - -3).^2 );
                idx = d_upper<crit_deg | d_lower<crit_deg;
                w(idx) = 0;
            end

            % 4. then also enforce any other rules that the user specifies:
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

        function [w, wStructural] = customStructureGuidedCandWeightingFunc(xy, convHullVerts, params)

            % -1. unpack input parameters structure
            % weight_type = params.weight_type;
            % excludeCandidatesAroundMidlines = params.excludeCandidatesAroundMidlines;
            % excludeCandidatesInOrAroundBlindspot = params.excludeCandidatesInOrAroundBlindspot;
            noiseLevel = params.noiseLevel; % 'high', 'low', 'none'
            isVeridical = params.isVeridical; % true / false

            % start with the default/vanilla weights
            w0 = MEDTEG.exampleCandWeightingFunc(xy, convHullVerts, params);

            % get structural info (weights 0 - 1)
            switch noiseLevel
                case 'none' % sigma = 0
                    wStructural = [
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % 27
                        0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.2, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % 25
                        0.0, 0.0, 0.0, 0.0, 0.1, 0.2, 0.3, 0.2, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % 23
                        0.0, 0.0, 0.0, 0.0, 0.2, 0.3, 0.5, 0.3, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % 21
                        0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.6, 0.4, 0.3, 0.2, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % 19
                        0.0, 0.2, 0.4, 0.7, 0.7, 0.8, 0.8, 0.8, 0.7, 0.7, 0.4, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % 17
                        0.0, 0.3, 0.7, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.7, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % 15
                        0.0, 0.3, 0.7, 1.0, 1.0, 1.0, 1.0, 0.9, 0.8, 0.7, 0.4, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % 13
                        0.0, 0.3, 0.7, 1.0, 1.0, 1.0, 1.0, 0.8, 0.6, 0.3, 0.2, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % 11
                        0.0, 0.3, 0.7, 1.0, 1.0, 1.0, 1.0, 0.7, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % 9
                        0.2, 0.4, 0.6, 0.8, 0.8, 0.7, 0.7, 0.4, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % 7
                        0.3, 0.4, 0.6, 0.7, 0.6, 0.4, 0.3, 0.2, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % 5
                        0.5, 0.5, 0.5, 0.5, 0.3, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % 3
                        0.3, 0.3, 0.3, 0.3, 0.2, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % 1
                        0.2, 0.2, 0.2, 0.2, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % -1
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % -3
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % -5
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % -7
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % -9
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % -11
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % -13
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % -15
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % -17
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % -19
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % -21
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % -23
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % -25
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 % -27
                        %-27, -25, -23, -21, -19, -17, -15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27,
                        ];
                case 'low' % sigma = 1/3
                    wStructural = [
                        0.07, 0.04, 0.00, 0.00, 0.00, 0.04, 0.06, 0.13, 0.33, 0.46, 0.30, 0.00, 0.36, 0.26, 0.21, 0.00, 0.40, 0.00, 0.00, 0.00, 0.17, 0.06, 0.00, 0.09, 0.16, 0.00, 0.26, 0.00 % 27
                        0.05, 0.00, 0.23, 0.11, 0.00, 0.60, 0.92, 0.30, 0.21, 0.00, 0.13, 0.00, 0.29, 0.00, 0.36, 0.31, 0.11, 0.00, 0.00, 0.00, 0.37, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 % 25
                        0.00, 0.62, 0.25, 0.00, 0.00, 0.04, 0.41, 0.00, 0.45, 0.37, 0.00, 0.06, 0.00, 0.79, 0.00, 0.00, 0.14, 0.38, 0.00, 0.00, 0.00, 0.79, 0.00, 0.01, 0.13, 0.00, 0.00, 0.00 % 23
                        0.54, 0.00, 0.36, 0.04, 0.00, 0.24, 0.23, 1.00, 0.10, 0.00, 0.00, 0.00, 0.00, 0.00, 0.23, 0.00, 0.12, 0.59, 0.08, 0.00, 0.00, 0.27, 0.00, 0.00, 0.56, 0.15, 0.19, 0.03 % 21
                        0.00, 0.86, 0.49, 0.16, 0.94, 0.13, 0.13, 0.37, 0.46, 0.57, 0.13, 0.45, 0.11, 0.00, 0.02, 0.03, 0.56, 0.59, 0.00, 0.00, 0.08, 0.00, 0.00, 0.00, 0.00, 0.00, 0.52, 0.00 % 19
                        0.00, 0.40, 1.00, 0.44, 1.00, 0.00, 0.89, 0.93, 0.65, 0.76, 0.00, 0.36, 0.64, 0.00, 0.00, 0.39, 0.30, 0.00, 0.00, 0.35, 0.00, 0.18, 0.27, 0.00, 0.06, 0.04, 0.01, 0.41 % 17
                        0.00, 0.66, 0.34, 1.00, 0.93, 0.79, 0.70, 0.84, 0.59, 1.00, 1.00, 0.26, 0.00, 0.22, 0.55, 0.01, 0.00, 0.62, 0.00, 0.00, 0.52, 0.00, 0.00, 0.07, 0.00, 0.00, 0.09, 0.00 % 15
                        0.14, 0.10, 0.56, 0.51, 0.96, 1.00, 1.00, 0.43, 0.59, 0.60, 0.27, 0.33, 0.00, 0.00, 0.11, 0.27, 0.04, 0.00, 0.13, 0.00, 0.00, 0.00, 0.19, 0.00, 0.00, 0.00, 0.00, 1.00 % 13
                        0.00, 0.12, 1.00, 1.00, 1.00, 0.87, 0.46, 0.77, 0.55, 0.41, 0.01, 0.42, 0.00, 0.26, 0.00, 0.00, 0.86, 0.00, 0.00, 0.00, 0.45, 0.00, 0.46, 0.00, 0.00, 0.00, 0.12, 0.00 % 11
                        0.07, 0.65, 1.00, 0.88, 0.25, 1.00, 0.73, 0.46, 0.30, 0.00, 0.21, 0.00, 0.00, 0.00, 0.00, 0.00, 0.16, 0.71, 0.08, 0.00, 0.00, 0.17, 0.01, 0.58, 0.06, 0.00, 0.00, 0.00 % 9
                        0.39, 0.39, 0.83, 1.00, 0.47, 0.69, 0.51, 0.18, 0.00, 0.29, 0.00, 0.00, 0.24, 0.00, 0.09, 0.28, 0.00, 0.14, 0.00, 0.04, 0.50, 0.00, 0.00, 0.00, 0.01, 0.00, 0.00, 0.00 % 7
                        0.76, 0.00, 0.62, 0.35, 0.70, 0.73, 0.70, 0.00, 0.00, 0.00, 0.08, 0.07, 0.00, 0.00, 0.00, 0.00, 0.00, 0.13, 0.00, 0.11, 0.00, 0.36, 0.00, 0.04, 0.00, 0.00, 0.00, 0.25 % 5
                        0.16, 0.40, 0.24, 0.60, 0.88, 0.64, 0.00, 0.00, 0.23, 0.00, 0.44, 0.00, 0.00, 0.00, 0.00, 0.11, 0.00, 0.00, 0.00, 0.15, 0.00, 0.30, 0.00, 0.20, 0.00, 0.27, 0.14, 0.29 % 3
                        0.00, 0.59, 0.52, 0.39, 0.22, 0.08, 0.40, 0.00, 0.00, 0.00, 0.00, 0.00, 0.46, 0.14, 1.00, 0.00, 0.00, 0.30, 0.45, 0.32, 0.00, 0.41, 0.00, 0.12, 0.00, 0.00, 0.44, 0.16 % 1
                        0.13, 0.24, 0.00, 0.45, 0.00, 0.24, 0.37, 0.25, 0.26, 0.33, 0.34, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.34, 0.00, 0.00, 0.00, 0.01, 0.00, 0.00, 0.48, 0.00, 0.00, 0.06 % -1
                        0.00, 0.00, 0.19, 0.00, 0.00, 0.32, 0.00, 0.01, 0.00, 0.84, 0.00, 0.00, 0.00, 0.18, 0.00, 0.16, 0.44, 0.00, 0.48, 0.18, 0.42, 0.33, 0.10, 0.15, 0.00, 0.33, 0.34, 0.02 % -3
                        0.56, 0.00, 0.21, 0.00, 0.15, 0.24, 0.00, 0.00, 0.18, 0.18, 0.00, 0.44, 0.05, 0.39, 0.39, 0.00, 0.19, 0.00, 0.00, 0.00, 0.00, 0.35, 0.12, 0.06, 0.54, 0.00, 0.00, 0.39 % -5
                        0.00, 0.27, 0.00, 0.58, 0.40, 0.08, 0.00, 0.00, 0.31, 0.51, 0.10, 0.00, 0.56, 0.46, 0.32, 0.25, 0.00, 0.21, 0.00, 0.24, 0.00, 0.00, 0.00, 0.45, 0.45, 0.02, 0.00, 0.00 % -7
                        0.00, 0.01, 0.00, 0.59, 0.35, 0.00, 0.00, 0.00, 0.00, 0.04, 0.00, 0.00, 0.00, 0.40, 0.00, 0.00, 0.50, 0.06, 0.00, 0.00, 0.12, 0.00, 0.37, 0.64, 0.00, 0.00, 0.00, 0.04 % -9
                        0.42, 0.00, 0.00, 0.42, 0.00, 0.00, 0.00, 0.24, 1.00, 0.35, 0.42, 0.00, 0.45, 0.12, 0.25, 0.00, 0.12, 0.29, 0.00, 0.08, 0.51, 0.26, 0.00, 0.00, 0.00, 0.10, 0.00, 0.00 % -11
                        0.00, 0.43, 0.37, 0.23, 0.00, 0.44, 0.00, 0.42, 0.00, 0.00, 0.00, 0.00, 0.00, 0.33, 0.33, 0.18, 0.52, 0.00, 0.00, 0.00, 0.27, 0.00, 0.31, 0.08, 0.30, 0.00, 0.00, 0.00 % -13
                        0.09, 0.00, 0.05, 0.19, 0.25, 0.03, 0.00, 0.00, 0.14, 0.40, 0.00, 0.52, 0.12, 0.00, 0.00, 0.03, 0.00, 0.06, 0.04, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.01, 0.33, 0.21 % -15
                        0.00, 0.00, 0.00, 0.19, 0.00, 0.00, 0.63, 0.01, 0.14, 0.00, 0.34, 0.00, 0.00, 0.29, 0.47, 0.29, 0.00, 0.64, 0.37, 0.48, 0.27, 0.00, 0.14, 0.08, 0.00, 0.00, 0.00, 0.56 % -17
                        0.67, 0.00, 0.06, 0.00, 0.00, 0.00, 0.00, 0.00, 0.08, 0.00, 0.18, 0.10, 0.00, 0.26, 0.00, 0.41, 0.41, 0.06, 0.00, 0.20, 0.38, 0.01, 0.23, 0.58, 0.00, 0.02, 0.40, 0.35 % -19
                        0.00, 0.00, 0.00, 0.37, 0.31, 0.15, 0.00, 0.31, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.36, 0.00, 0.00, 0.43, 0.00, 0.03, 0.00, 0.30, 0.24, 0.60, 0.70, 0.00, 0.01 % -21
                        0.00, 0.00, 0.00, 0.00, 0.02, 0.27, 0.00, 0.07, 0.11, 0.00, 0.00, 0.00, 0.00, 0.06, 0.00, 0.11, 0.00, 0.01, 0.43, 0.00, 0.11, 0.00, 0.00, 0.00, 0.00, 0.00, 0.14, 0.14 % -23
                        0.21, 0.00, 0.00, 0.43, 0.13, 0.00, 0.00, 0.26, 0.00, 0.00, 0.10, 0.20, 0.40, 0.53, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.16, 0.00, 0.35, 0.33, 0.19, 0.41 % -25
                        0.00, 0.00, 0.02, 0.00, 0.13, 0.07, 0.05, 0.00, 0.21, 0.08, 0.00, 0.00, 0.00, 0.00, 0.00, 0.07, 0.00, 0.00, 0.00, 0.20, 0.00, 0.00, 0.00, 0.22, 0.48, 0.00, 0.41, 0.00 % -27
                        ]; %-27, -25, -23, -21, -19, -17, -15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27,
                case 'medium' % sigma = 2/3
                    wStructural = [
                        0.00, 0.00, 0.00, 0.00, 0.88, 0.00, 0.00, 0.49, 0.00, 0.00, 0.02, 0.94, 0.16, 0.18, 1.00, 0.25, 0.35, 0.00, 0.30, 1.00, 0.00, 1.00, 0.18, 0.59, 0.12, 0.15, 0.00, 0.00 % 27
                        0.01, 0.00, 0.00, 0.00, 0.00, 0.78, 0.30, 0.00, 0.01, 0.00, 0.00, 0.00, 0.00, 0.43, 0.39, 0.18, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.37, 0.00, 0.00, 0.00, 0.26 % 25
                        0.00, 0.00, 0.00, 0.00, 0.56, 0.03, 1.00, 0.13, 0.65, 0.00, 0.06, 0.00, 0.00, 0.00, 0.00, 1.00, 0.13, 0.37, 0.15, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.58, 0.41, 0.68 % 23
                        0.33, 0.00, 0.00, 0.00, 0.69, 0.47, 0.44, 0.97, 0.55, 0.00, 0.00, 0.00, 0.26, 0.01, 0.00, 0.85, 0.26, 1.00, 1.00, 0.00, 0.29, 0.00, 0.09, 0.00, 0.24, 0.31, 0.19, 0.07 % 21
                        1.00, 0.00, 1.00, 0.36, 0.95, 0.00, 0.15, 0.00, 1.00, 0.33, 0.00, 0.20, 0.61, 0.41, 1.00, 0.00, 1.00, 0.00, 0.63, 0.33, 0.45, 0.00, 0.00, 0.14, 0.78, 0.10, 0.00, 0.68 % 19
                        0.40, 1.00, 0.00, 0.23, 1.00, 0.00, 0.27, 0.00, 0.58, 0.65, 0.95, 1.00, 0.00, 1.00, 0.57, 0.00, 0.00, 0.10, 0.00, 0.00, 0.57, 0.09, 0.53, 0.00, 1.00, 0.87, 0.54, 0.38 % 17
                        0.02, 0.29, 0.67, 1.00, 1.00, 1.00, 0.00, 1.00, 1.00, 1.00, 0.00, 0.07, 0.00, 0.00, 0.00, 0.00, 0.24, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.43, 0.58, 0.71 % 15
                        0.22, 1.00, 0.74, 1.00, 0.94, 0.60, 1.00, 0.97, 0.35, 0.91, 0.10, 0.72, 0.28, 0.00, 0.17, 0.72, 1.00, 0.00, 0.00, 0.00, 0.00, 0.64, 0.00, 0.00, 0.52, 0.36, 0.06, 0.79 % 13
                        1.00, 0.36, 1.00, 1.00, 1.00, 0.88, 0.99, 0.71, 0.00, 0.08, 0.05, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.06, 0.00, 0.00, 0.00, 0.26, 0.34 % 11
                        0.00, 1.00, 1.00, 0.96, 0.71, 0.99, 0.29, 0.93, 0.62, 0.00, 0.61, 0.00, 0.78, 0.00, 0.00, 0.77, 1.00, 0.48, 0.92, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.32 % 9
                        0.00, 0.00, 0.00, 0.27, 1.00, 0.36, 0.64, 0.03, 1.00, 0.57, 0.00, 0.00, 0.00, 0.23, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 1.00, 0.74, 0.57, 0.79, 0.00, 0.81, 1.00, 0.00 % 7
                        0.00, 1.00, 0.57, 0.94, 0.54, 0.00, 0.71, 0.54, 0.00, 0.17, 0.00, 1.00, 1.00, 0.00, 0.00, 0.03, 0.00, 0.19, 0.00, 0.00, 0.03, 0.18, 0.00, 0.66, 0.00, 0.46, 0.00, 0.00 % 5
                        1.00, 0.92, 0.00, 1.00, 0.73, 0.00, 0.00, 1.00, 1.00, 0.94, 0.00, 0.84, 0.87, 0.29, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.77, 0.17, 0.13, 0.00, 0.00, 0.00, 0.66, 1.00 % 3
                        0.17, 0.00, 0.00, 0.84, 0.47, 0.00, 0.00, 0.00, 0.33, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.19, 0.00, 0.59, 0.90, 0.66, 0.00, 0.21, 0.00, 0.00, 0.58, 0.28, 0.02, 0.01 % 1
                        0.99, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.27, 0.00, 0.00, 0.01, 0.92, 0.22, 0.52, 0.00, 0.00, 0.46, 0.03, 1.00, 0.00, 0.00, 0.42, 0.00, 0.00, 0.04, 1.00, 0.00 % -1
                        0.00, 0.00, 0.02, 0.40, 0.17, 0.50, 0.00, 0.96, 0.00, 0.00, 0.16, 1.00, 0.79, 0.00, 0.00, 0.85, 0.00, 0.00, 0.10, 0.00, 0.00, 0.03, 0.00, 0.00, 0.79, 0.71, 1.00, 0.52 % -3
                        0.69, 0.00, 0.00, 1.00, 1.00, 1.00, 0.00, 0.00, 0.99, 0.09, 0.00, 0.00, 0.44, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.21, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.02 % -5
                        0.00, 0.11, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.91, 1.00, 0.00, 0.00, 0.00, 0.73, 0.00, 0.54, 0.00, 0.00, 0.00, 0.89, 0.94, 0.00, 1.00, 0.26, 0.21, 0.14, 0.00 % -7
                        0.00, 1.00, 0.00, 0.71, 0.11, 0.00, 0.18, 0.00, 0.05, 0.00, 0.00, 0.84, 0.14, 0.71, 0.00, 0.00, 0.00, 0.62, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.84, 1.00, 0.30, 0.41 % -9
                        1.00, 0.15, 0.32, 0.00, 0.00, 0.44, 0.09, 0.71, 0.12, 0.30, 0.00, 0.15, 0.40, 0.00, 1.00, 0.00, 0.29, 0.00, 0.77, 0.00, 0.00, 0.00, 0.76, 0.20, 0.14, 0.70, 0.05, 0.61 % -11
                        0.77, 0.20, 0.00, 0.00, 0.10, 0.68, 0.21, 0.00, 0.40, 0.29, 0.00, 0.04, 0.32, 0.26, 0.00, 0.00, 0.00, 1.00, 1.00, 0.00, 0.00, 0.00, 0.15, 0.31, 0.00, 0.01, 0.65, 0.00 % -13
                        0.00, 0.56, 0.00, 0.72, 0.00, 0.00, 0.86, 0.48, 0.00, 0.67, 1.00, 0.18, 0.00, 0.70, 0.00, 0.38, 0.00, 0.31, 0.00, 0.11, 0.00, 0.00, 0.00, 0.00, 0.11, 0.03, 0.00, 0.00 % -15
                        0.47, 0.90, 0.17, 0.96, 0.11, 0.71, 0.27, 0.77, 0.00, 0.40, 0.00, 0.00, 0.37, 0.89, 0.00, 0.19, 1.00, 0.00, 0.00, 0.53, 0.01, 0.00, 0.15, 0.39, 0.00, 0.89, 0.32, 0.00 % -17
                        0.00, 0.00, 0.47, 0.00, 0.00, 1.00, 0.00, 0.65, 0.00, 0.53, 0.47, 0.00, 1.00, 0.00, 0.29, 0.01, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.36, 0.00, 0.38, 0.00, 0.10, 0.25 % -19
                        0.00, 0.00, 0.00, 0.00, 0.36, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.38, 0.11, 0.39, 0.00, 0.00, 0.00, 0.00, 0.34, 0.81, 0.34, 0.25, 0.40, 0.00, 0.00, 0.62, 0.70, 0.33 % -21
                        0.31, 0.11, 0.50, 0.00, 1.00, 0.25, 0.89, 0.34, 0.61, 0.00, 0.03, 0.51, 0.42, 0.00, 0.22, 0.53, 0.20, 0.00, 0.24, 0.24, 0.00, 0.00, 0.00, 0.83, 1.00, 0.47, 0.75, 0.00 % -23
                        0.08, 0.11, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.13, 0.00, 0.04, 0.00, 0.00, 0.00, 0.00, 0.00, 0.04, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.18, 0.73, 0.73, 0.15 % -25
                        0.02, 0.74, 0.56, 0.00, 0.25, 0.00, 1.00, 0.00, 0.99, 0.00, 0.98, 0.00, 0.69, 0.00, 0.19, 0.00, 0.00, 1.00, 0.29, 0.07, 0.00, 0.81, 0.53, 0.00, 0.00, 0.00, 1.00, 0.00 % -27
                        ]; %-27, -25, -23, -21, -19, -17, -15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27,
                case 'high' % sigma = 1
                    wStructural = [
                        0.00, 0.62, 0.63, 0.00, 0.00, 0.00, 0.00, 0.74, 1.00, 0.41, 0.00, 0.14, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 1.00, 0.57, 0.00, 1.00, 0.00, 0.18, 0.00, 0.20, 0.00, 0.34 % 27
                        0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.37, 0.61, 0.00, 0.00, 0.29, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.04, 0.00, 0.24, 0.65, 0.53, 1.00, 0.30, 0.54 % 25
                        0.83, 1.00, 0.00, 0.00, 0.16, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.78, 0.30, 0.50, 0.00, 0.00, 0.51, 0.96, 0.06, 0.72, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00 % 23
                        0.00, 0.13, 0.00, 0.11, 0.00, 1.00, 1.00, 1.00, 0.00, 0.80, 1.00, 0.00, 0.00, 0.75, 0.90, 0.82, 0.30, 0.65, 1.00, 0.00, 0.44, 0.00, 0.47, 0.61, 0.92, 0.49, 0.00, 0.42 % 21
                        1.00, 0.00, 0.00, 0.73, 0.74, 0.38, 1.00, 1.00, 1.00, 0.00, 0.00, 0.71, 0.00, 1.00, 0.00, 0.84, 0.18, 0.00, 0.56, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 1.00, 0.00 % 19
                        0.00, 0.58, 1.00, 0.00, 0.67, 0.00, 0.83, 0.23, 1.00, 0.98, 0.00, 0.00, 1.00, 0.00, 0.03, 0.17, 0.36, 0.71, 1.00, 0.00, 0.00, 0.93, 0.00, 0.00, 0.19, 0.25, 0.92, 0.00 % 17
                        0.39, 0.41, 0.60, 1.00, 0.45, 0.00, 1.00, 0.61, 1.00, 1.00, 0.70, 0.41, 0.00, 0.39, 1.00, 0.00, 0.00, 0.00, 0.00, 0.67, 0.00, 0.80, 0.00, 0.40, 0.40, 0.78, 1.00, 0.00 % 15
                        0.92, 0.35, 0.00, 0.00, 0.62, 1.00, 0.99, 1.00, 1.00, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.42, 0.00, 0.31, 1.00, 0.21, 0.08, 0.00, 1.00 % 13
                        0.00, 1.00, 1.00, 1.00, 0.00, 0.46, 0.00, 1.00, 1.00, 0.91, 0.00, 0.64, 0.65, 0.79, 0.37, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.09, 0.83, 0.00, 1.00, 0.09, 0.36, 0.00 % 11
                        0.09, 0.56, 0.00, 0.00, 1.00, 0.00, 0.00, 0.12, 0.00, 0.00, 1.00, 0.00, 0.31, 0.00, 1.00, 0.00, 0.22, 0.93, 0.00, 1.00, 0.00, 0.91, 0.00, 0.15, 0.93, 1.00, 1.00, 0.00 % 9
                        0.00, 0.40, 1.00, 1.00, 0.78, 0.90, 0.69, 1.00, 0.00, 0.00, 0.53, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.90, 0.00, 0.00, 0.38, 0.00, 1.00, 0.00, 0.61, 0.47, 1.00 % 7
                        0.82, 1.00, 1.00, 0.62, 0.00, 0.08, 1.00, 0.00, 0.83, 0.00, 0.00, 0.25, 1.00, 0.21, 0.00, 0.85, 0.16, 0.88, 0.48, 0.00, 0.00, 0.00, 0.00, 0.89, 0.00, 0.69, 0.00, 0.40 % 5
                        1.00, 0.17, 1.00, 1.00, 0.47, 1.00, 0.00, 0.68, 0.00, 0.22, 0.00, 0.00, 1.00, 0.00, 1.00, 0.64, 0.00, 0.00, 0.45, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.95, 0.00, 0.00 % 3
                        0.09, 0.00, 0.00, 0.59, 0.29, 0.27, 0.00, 0.28, 0.82, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.86, 0.00, 0.00, 1.00, 0.00, 0.00, 0.20, 0.00, 0.00 % 1
                        1.00, 0.00, 1.00, 0.00, 1.00, 0.00, 0.00, 0.50, 0.12, 0.00, 0.00, 1.00, 0.24, 0.00, 0.59, 0.00, 0.00, 0.00, 0.00, 0.00, 0.87, 0.00, 0.31, 0.00, 0.00, 0.00, 1.00, 0.00 % -1
                        0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 1.00, 1.00, 0.66, 1.00, 0.78, 0.70, 0.00, 0.50, 1.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.30, 0.35, 1.00, 0.00, 0.00, 0.00, 1.00, 0.97 % -3
                        0.00, 1.00, 0.00, 1.00, 0.00, 0.00, 0.28, 0.00, 0.00, 1.00, 0.49, 0.00, 0.00, 1.00, 0.41, 0.00, 1.00, 0.04, 0.43, 0.00, 0.00, 0.00, 0.00, 0.43, 0.00, 0.46, 1.00, 0.89 % -5
                        0.50, 0.69, 0.99, 0.00, 1.00, 0.00, 0.72, 0.00, 0.08, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.51, 0.17, 0.00, 0.12, 0.00, 0.19, 0.00, 0.00, 0.32, 0.00 % -7
                        0.00, 0.79, 0.68, 0.00, 0.15, 1.00, 1.00, 0.12, 0.00, 0.54, 0.00, 0.00, 1.00, 1.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.45, 0.00, 0.25, 0.00, 0.12, 0.00, 1.00, 0.29 % -9
                        0.00, 0.71, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.10, 0.00, 0.00, 0.91, 0.07, 1.00, 0.00, 0.00, 0.00, 0.51, 0.31, 0.04, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 % -11
                        0.44, 0.48, 1.00, 0.00, 1.00, 0.44, 0.50, 1.00, 0.32, 0.99, 0.00, 0.00, 0.70, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.91, 1.00, 1.00, 1.00, 0.00, 0.22, 0.00, 0.45 % -13
                        0.44, 0.35, 0.00, 0.00, 1.00, 0.00, 0.73, 0.00, 1.00, 0.57, 0.80, 0.12, 0.43, 0.92, 1.00, 1.00, 0.00, 0.00, 0.87, 0.00, 0.51, 0.59, 0.00, 0.66, 0.00, 0.02, 1.00, 0.00 % -15
                        0.00, 0.00, 0.45, 0.30, 1.00, 0.00, 0.00, 0.00, 0.00, 0.21, 0.02, 0.00, 0.00, 0.00, 0.00, 0.00, 0.27, 0.80, 0.31, 0.36, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00 % -17
                        0.20, 0.51, 0.00, 0.06, 0.00, 1.00, 0.00, 0.00, 0.27, 0.42, 0.00, 0.00, 0.00, 1.00, 1.00, 0.09, 0.62, 1.00, 1.00, 0.82, 0.87, 0.00, 1.00, 0.00, 0.53, 0.61, 0.00, 0.00 % -19
                        0.46, 1.00, 0.00, 0.12, 1.00, 0.00, 0.07, 0.00, 0.00, 0.00, 0.19, 0.00, 0.81, 0.00, 0.17, 0.51, 0.73, 0.00, 0.00, 0.04, 0.93, 0.37, 0.00, 0.56, 1.00, 0.26, 0.53, 0.00 % -21
                        0.00, 0.85, 1.00, 0.12, 1.00, 0.00, 0.96, 0.00, 0.00, 0.00, 0.00, 0.53, 0.00, 0.00, 0.00, 0.11, 0.18, 0.00, 0.14, 0.00, 0.00, 0.00, 1.00, 0.00, 1.00, 0.00, 0.00, 0.00 % -23
                        0.00, 0.09, 0.00, 1.00, 0.00, 0.74, 1.00, 0.09, 0.00, 0.00, 0.00, 0.00, 1.00, 1.00, 0.87, 0.00, 1.00, 0.41, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.63, 1.00 % -25
                        0.00, 1.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.09, 0.42, 0.00, 0.64, 0.24, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.18, 0.00, 1.00, 1.00, 0.39, 0.08, 0.00, 0.00, 1.00 % -27
                        ]; %-27, -25, -23, -21, -19, -17, -15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27,
                otherwise
                    error('Specified param ("noiseLevel") value not recognised: %s', noiseLevel);
            end
            % wStructural = fliplr(wStructural); % try flipping to confirm not a fluke!
            [x_deg, y_deg] = meshgrid(linspace(-27, 27, 28), linspace(27, -27, 28));
            
            % add a modest baseline to all values (otherwise may end up
            % with no valid locations, particularly early on when N
            % candidates may be low)
            wStructural = wStructural + .01;

            % simulate technician flipping the eye(!!)
            if ~isVeridical
                wStructural = fliplr(wStructural);
                % wStructural = flipud(wStructural);
            end

            w1 = interp2(x_deg, y_deg, wStructural, xy(:,1), xy(:,2));


            % 0 any points outside of measured region
            w1(isnan(w1)) = 0;

            % combine structure with default
            w = w0 .* w1.^0.5;
        end

        function [x_deg, y_deg, norm_DLS, norm_sigma] = getExampleNormativeValues()
            % example normative values. In this case extracted from Brenton
            % & Phelps (1986) The normal visual field on the Humphrey field
            % analyzer, Ophthalmologica, 193(1-2):56-74. doi: 10.1159/000309679.
            norm_DLS = [ % right eye format (OD)
                NaN  NaN  NaN  26.8 26.3 25.5 25.1 NaN  NaN  NaN    % +27
                NaN  NaN  29.0 29.2 28.3 27.6 27.8 28.0 NaN  NaN    % +21
                NaN  28.9 29.1 29.8 30.2 29.8 30.9 30.4 28.8 NaN    % +15
                28.3 29.6 31.2 32.6 32.3 31.7 31.2 30.6 30.4 30.3   % +9
                28.9 30.4 32.1 32.8 33.7 33.5 32.2   0  31.3 31.9   % +3
                29.4 31.0 32.6 33.3 33.8 33.8 32.7   0  31.5 31.3   % -3
                29.6 30.4 31.5 33.4 32.9 32.8 32.8 31.1 31.6 31.1   % -9
                NaN  30.5 30.6 31.2 32.1 31.3 30.9 32.4 31.2 NaN    % -15
                NaN  NaN  30.1 29.8 30.9 31.4 31.8 31.2 NaN  NaN    % -21
                NaN  NaN  NaN  29.0 30.1 30.8 30.7 NaN  NaN  NaN    % -27
                ]; %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
            norm_sigma = [
                NaN  NaN  NaN  5.8  3.5  4.7  3.9  NaN  NaN  NaN    % +27
                NaN  NaN  3.6  2.6  4.1  2.6  2.6  2.0  NaN  NaN    % +21
                NaN  2.9  3.2  2.7  2.2  2.4  2.2  2.3  2.7  NaN    % +15
                2.7  2.4  1.6  1.7  1.6  2.1  2.1  2.7  2.8  3.2    % +9
                2.5  2.2  1.7  1.5  2.0  1.7  1.9   3  2.9  2.4     % +3
                2.3  1.7  1.4  2.2  1.4  1.4  2.1   3  2.5  2.0     % -3
                2.2  2.1  2.5  1.7  1.9  1.8  2.2  3.2  2.0  2.6    % -9
                NaN  1.8  1.8  2.8  2.1  2.0  2.3  1.9  2.4  NaN    % -15
                NaN  NaN  2.5  2.4  2.1  1.9  2.1  2.9  NaN  NaN    % -21
                NaN  NaN  NaN  3.0  2.2  2.6  2.4  NaN  NaN  NaN    % -27
                ]; %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
            x_deg = ones(10,1) * (-27:6:27);
            y_deg = (27:-6:-27)' * ones(1,10);
        end

        function [x_deg, y_deg, simObs_true_DLS, F_simObs_DLS, dat] = runExample(exampleN, optionalParams)
            % Examples of use, demonstrating/testing key functionalities.
            % Examples include:
            %   1. Simple, "enhanced 24-2" case, where we test the 24-2 and then test an additional N points using MEDTEG
            %
            % @param    exampleN	Example to run [1|2|3|4|5|6|7|8|9]. Defaults to 1.
            %
            % @date     07/09/24
            % @author   PRJ
            %

            % parse inputs: if no example specified, run example 1
            if nargin<1 || isempty(exampleN)
                fprintf('Defaulting to example 0..\n');
                exampleN = 0;
            end

            % run selected example
            switch exampleN
                case 0 % Simple, "enhanced 24-2" case with made-up data, where we generate a 24-2 of pre-tested points and then request an additional 1 location using MEDTEG
                    
                    % define simulated observer (try setting one item to
                    % an extreme value, e.g., 0, to see how the algorithm
                    % will automatically attempt to test around there)
                    simObs_true_DLS = [ % right eye format (OD)
                        NaN  NaN  29.0 29.2 28.3 27.6 27.8 28.0 NaN  NaN    % +21
                        NaN  28.9 29.1 29.8 30.2 29.8 30.9 30.4 28.8 NaN    % +15
                        28.3 29.6 31.2 32.6 32.3 31.7 31.2 30.6 30.4 30.3   % +9
                        28.9 30.4 32.1 32.8 33.7 33.5 32.2   0  31.3 31.9   % +3
                        29.4 31.0 32.6 33.3 33.8 33.8 32.7   0  31.5 31.3   % -3
                        29.6 30.4 31.5 33.4 32.9 32.8 32.8 31.1 31.6 31.1   % -9
                        NaN  30.5 30.6 31.2 32.1 31.3 30.9 32.4 31.2 NaN    % -15
                        NaN  NaN  30.1 29.8 30.9 31.4 31.8 31.2 NaN  NaN    % -21
                    ]; %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27

                    % get standard Q+ params
                    QPlusParams = MEDTEG.getExampleQuestPlusParams();

                    % define fixed, starting test grid
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
                    idx = ~isnan(id); % exclude null points
                    x = x_deg(idx);
                    % y = y_deg(idx);

                    % make up fake data (create QUEST+ objects)
                    existingTestLocations_XY = nan(length(x), 2);
                    existingTestLocations_QPlus = cell(1, length(x));
                    for i = 1:length(x)
                        existingTestLocations_XY(i,:) = [x_deg(id==i) y_deg(id==i)];
                        existingTestLocations_QPlus{i} = QuestPlus(QPlusParams.F, QPlusParams.stimDomain, QPlusParams.paramDomain, QPlusParams.respDomain, QPlusParams.stopRule, QPlusParams.stopCriterion, QPlusParams.minNTrials, QPlusParams.maxNTrials);
                        PMF = normpdf(QPlusParams.paramDomain, simObs_true_DLS(id == i), 5); % 5 <-- assume constant measurement error across the VF (not true!)
                        PMF = PMF./sum(PMF);
                        existingTestLocations_QPlus{i}.initialise(PMF, QPlusParams.likelihoodsDat);
                    end

                    % create MEDTEG object
                    M = MEDTEG(@MEDTEG.exampleCandWeightingFunc, MEDTEG.getDefaultCandWeightingParameters(), QPlusParams);

                    % now get an additional 1 point using MEDTEG
                    [bestCandidateXY, bestCandidateQPlus, ~, candidates, ~, cand_wEDHdeg2] = M.GetNewTestLocation(existingTestLocations_XY, existingTestLocations_QPlus);

                    % store
                    newTestLocations_XY = bestCandidateXY; 
                    newTestLocations_QPlus = {bestCandidateQPlus};
                case 1 % Monte Carlo "enhanced 24-2" case with simulated trial-by-trial responses, in which we first test the standard 24-2 grid and then test an additional N=5 points using MEDTEG (NB: "Simulation 1" in paper)
                    % define normative values (used for setting priors)
                    [x_deg, y_deg, norm_DLS, norm_sigma] = MEDTEG.getExampleNormativeValues();

                    % init simulation params
                    QPlusParams = MEDTEG.getExampleQuestPlusParams();%  get standard Q+ params
                    totalNTrialsCounter = 0;
                    n = 5; % n additional test locations to add

                    % define simulated observer
                    simObs_true_DLS = [ % right eye format (OD) - with a small, deep, central scotoma
                        NaN  NaN  NaN  26.8 26.3 25.5 25.1 NaN  NaN  NaN    % +27
                        NaN  NaN  29.0 29.2 28.3 27.6 27.8 28.0 NaN  NaN    % +21
                        NaN  28.9 29.1 29.8 30.2 29.8 30.9 30.4 28.8 NaN    % +15
                        28.3 29.6 31.2  2.6  2.3 31.7 31.2 30.6 30.4 30.3   % +9
                        28.9 30.4 32.1  2.8 33.7 33.5 32.2  0   31.3 31.9   % +3
                        29.4 31.0 32.6 33.3 33.8 33.8 32.7  0   31.5 31.3   % -3
                        29.6 30.4 31.5 33.4 32.9 32.8 32.8 31.1 31.6 31.1   % -9
                        NaN  30.5 30.6 31.2 32.1 31.3 30.9 32.4 31.2 NaN    % -15
                        NaN  NaN  30.1 29.8 30.9 31.4 31.8 31.2 NaN  NaN    % -21
                        NaN  NaN  NaN  29.0 30.1 30.8 30.7 NaN  NaN  NaN    % -27
                    ]; %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27                    
                    simObs_true_DLS = simObs_true_DLS + randn(size(simObs_true_DLS)).*norm_sigma; % also add some random jitter to all points, just for some added complexity
                    simObs_true_sigma = [ % slope of psychometric function
                        NaN  NaN  NaN  5.8  3.5  4.7  3.9  NaN  NaN  NaN    % +27
                        NaN  NaN  3.6  2.6  4.1  2.6  2.6  2.0  NaN  NaN    % +21
                        NaN  2.9  3.2  2.7  2.2  2.4  2.2  2.3  2.7  NaN    % +15
                        2.7  2.4  1.6  1.7  1.6  2.1  2.1  2.7  2.8  3.2    % +9
                        2.5  2.2  1.7  1.5  2.0  1.7  1.9   3   2.9  2.4     % +3
                        2.3  1.7  1.4  2.2  1.4  1.4  2.1   3   2.5  2.0     % -3
                        2.2  2.1  2.5  1.7  1.9  1.8  2.2  3.2  2.0  2.6    % -9
                        NaN  1.8  1.8  2.8  2.1  2.0  2.3  1.9  2.4  NaN    % -15
                        NaN  NaN  2.5  2.4  2.1  1.9  2.1  2.9  NaN  NaN    % -21
                        NaN  NaN  NaN  3.0  2.2  2.6  2.4  NaN  NaN  NaN    % -27
                    ]; %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
                    simObs_true_lambda = 0.05; % lapse rate
                    simObs_true_gamma = 0.02; % guess rate
                    simObs_true_F = @(x, mu,sigma, lambda, gamma)(gamma + (1-gamma-lambda) * (1-normcdf(x,mu,sigma)  )); % "probability of anscorrect" psychometric function (v. similar in form to the one assumed by QUEST+, but here for ease of manipulation we shall make all the variables free, including the slope of the psychometric function

                    % define fixed, starting test grid
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
                    idx = ~isnan(ids); % exclude null points
                    x = x_deg(idx);
                    % y = y_deg(idx);

                    % run simulations
                    existingTestLocations_XY = nan(length(x), 2);
                    existingTestLocations_QPlus = cell(1, length(x));
                    k = 0;
                    for id = unique(ids(idx))'
                        k = k + 1;

                        % get/set xy location
                        existingTestLocations_XY(k,:) = [x_deg(ids == id) y_deg(ids == id)];
 
                        % normative values (used to initialise the PMF
                        % prior, and we shall also use a jittered version
                        % of these values for the simulated observer)
                        prior_DLS = norm_DLS(ids == id);
                        prior_sigma = norm_sigma(ids == id);

                        % init psychophysical algorithm
                        QP = QuestPlus(QPlusParams.F, QPlusParams.stimDomain, QPlusParams.paramDomain, QPlusParams.respDomain, QPlusParams.stopRule, QPlusParams.stopCriterion, QPlusParams.minNTrials, QPlusParams.maxNTrials);
                        prior = MEDTEG.getExamplePrior(QPlusParams.paramDomain, prior_DLS, prior_sigma*1.96);
                        QP.initialise(prior, QPlusParams.likelihoodsDat);

                        % get sim observer true sensitivity (jittered versions of normative value)
                        obs_DLS = simObs_true_DLS(ids == id);
                        obs_sigma = simObs_true_sigma(ids == id);
                        obs_lambda = simObs_true_lambda;
                        obs_gamma = simObs_true_gamma;

                        % run trials
                        startGuess = QP.getParamEsts('mode'); %#ok starting guess (for debugging)
                        while ~QP.isFinished() % for trial = 1:6
                            dB = QP.getTargetStim();
                            pC = simObs_true_F(dB, obs_DLS, obs_sigma, obs_lambda, obs_gamma); % QPlusParams.F(dB, obs_true_DLS); % get probability of answering correctly
                            anscorrect = rand() < pC;
                            QP.update(dB, anscorrect);
                            totalNTrialsCounter = totalNTrialsCounter + 1;
                        end
                        endGuess = QP.getParamEsts('mode'); %#ok final estimate (for debugging)
                        % fprintf('<x=%i, y=%i>  True Sensitivity: %1.1f dB;  Start Guess = %1.1f dB;  Final Estimate: %1.1f\n', x_deg(ids==id), y_deg(ids==id), obs_true_DLS, startGuess, endGuess)% <- for debugging. Uncomment to see that Q+ really  does get closer to the true answer over time

                        % store result
                        existingTestLocations_QPlus{k} = QP;
                    end

                    % create MEDTEG object
                    weightParams = MEDTEG.getDefaultCandWeightingParameters();
                    weightParams.weight_type = 'flat'; % 'macula_only', 'prefer_central', 'flat'
                    M = MEDTEG(@MEDTEG.exampleCandWeightingFunc, weightParams, QPlusParams);

                    % create interpolants for simulated-observer data (used below)
                    % idx = ~isnan(norm_DLS) & norm_DLS > 1; % exclude untested points and points in physiologic blindspots
                    idx = ~isnan(norm_DLS); % exclude untested points
                    x = x_deg(idx);
                    y = y_deg(idx);
                    v = simObs_true_DLS(idx);
                    F_simObs_DLS = scatteredInterpolant(x,y,v);
                    v = simObs_true_sigma(idx);
                    F_simObs_sigma = scatteredInterpolant(x,y,v);

                    % similarly, create interpolants for normative data (used below)
                    v = norm_DLS(idx);
                    F_norm_DLS = scatteredInterpolant(x,y,v);

                    % now get and test an additional N points using MEDTEG
                    for i = 1:n
                        [bestCandidateXY, bestCandidateQPlus, ~, candidates, ~, cand_wEDHdeg2] = M.GetNewTestLocation(existingTestLocations_XY, existingTestLocations_QPlus);

                        % init psychophysical algorithm
                        QP = QuestPlus(QPlusParams.F, QPlusParams.stimDomain, QPlusParams.paramDomain, QPlusParams.respDomain, QPlusParams.stopRule, QPlusParams.stopCriterion, QPlusParams.minNTrials, QPlusParams.maxNTrials);
                        prior_DLS = bestCandidateQPlus.getParamEsts('mode');
                        prior_sigma = 5;
                        prior = MEDTEG.getExamplePrior(QPlusParams.paramDomain, prior_DLS, prior_sigma); % preferable to the raw prior computed by MEDTEG
                        QP.initialise(prior, QPlusParams.likelihoodsDat);

                        % get sim observer true sensitivity by interpolation
                        x = bestCandidateXY(1);
                        y = bestCandidateXY(2);
                        obs_DLS = F_simObs_DLS(x, y);
                        obs_sigma = F_simObs_sigma(x, y);
                        obs_lambda = 0.05; % lapse rate
                        obs_gamma = 0.02; % guess rate

                        % run trials
                        startGuess = QP.getParamEsts('mode'); % starting guess (for debugging)
                        while ~QP.isFinished() % for trial = 1:6
                            dB = QP.getTargetStim();
                            pC = simObs_true_F(dB, obs_DLS, obs_sigma, obs_lambda, obs_gamma); % QPlusParams.F(dB, obs_true_DLS); % get probability of answering correctly
                            anscorrect = rand() < pC;
                            QP.update(dB, anscorrect);
                            totalNTrialsCounter = totalNTrialsCounter + 1;
                        end
                        endGuess = QP.getParamEsts('mode'); % final estimate (for debugging)
                        fprintf('<x=%1.0f, y=%1.0f>  True Sensitivity: %1.1f dB;  Start Guess = %1.1f dB;  Final Estimate: %1.1f\n', x, y, obs_DLS, startGuess, endGuess)% <- for debugging. Uncomment to see that Q+ really  does get closer to the true answer over time

                        % store result
                        existingTestLocations_XY(end+1, :) = bestCandidateXY; %#ok
                        existingTestLocations_QPlus{end+1} = QP; %#ok
                    end

                    % for plotting purposes we'll explicitly store the original and new test locations separately
                    newTestLocations_XY                             = existingTestLocations_XY((end-(n-1)):end, :);
                    newTestLocations_QPlus                          = existingTestLocations_QPlus((end-(n-1)):end);
                    existingTestLocations_XY((end-(n-1)):end, :)    = [];
                    existingTestLocations_QPlus((end-(n-1)):end)    = [];

                case 2 % Monte Carlo "seed" case with simulated trial-by-trial responses, in which we first test just a handful of seed locations, and then build a dynamic grid up (NB: "Simulation 2" in paper)
                    % define normative values (used for setting priors)
                    [x_deg, y_deg, norm_DLS, norm_sigma] = MEDTEG.getExampleNormativeValues();

                    % init simulation params
                    QPlusParams = MEDTEG.getExampleQuestPlusParams();%  get standard Q+ params
                    totalNTrialsCounter = 0;
                    n = 14; % n additional test locations to add (44 + 10 seeds, to bring up to parity with 24-2)

                    % define simulated observer
                    simObs_true_DLS = [ % right eye format (OD) - with an extensive superior hemifield loss
                        NaN   NaN  NaN  0.8  0.3  0.5  0.1  NaN  NaN  NaN    % +27
                        NaN   NaN  0.0  0.2  0.3  0.6  0.8  0.0  NaN  NaN    % +21
                        NaN   0.9  0.1  0.8  0.2  2.8  1.9  0.4  0.8  NaN    % +15
                         1.3  1.6  1.2  2.6  2.3  1.7  1.2  0.6  0.4  0.3   % +9
                         2.9  2.4  2.1  2.8  2.7  2.5  2.2  0    1.3  1.9   % +3
                        29.4 31.0 32.6 33.3 33.8 33.8 32.7  0   31.5 31.3   % -3
                        29.6 30.4 31.5 33.4 32.9 32.8 32.8 31.1 31.6 31.1   % -9
                        NaN  30.5 30.6 31.2 32.1 31.3 30.9 32.4 31.2 NaN    % -15
                        NaN  NaN  30.1 29.8 30.9 31.4 31.8 31.2 NaN  NaN    % -21
                        NaN  NaN  NaN  29.0 30.1 30.8 30.7 NaN  NaN  NaN    % -27
                    ]; %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27                    
                    simObs_true_DLS = simObs_true_DLS + randn(size(simObs_true_DLS)).*norm_sigma; % also add some random jitter to all points, just for some added complexity
                    if rand() > 0.5
                        simObs_true_DLS = flipud(simObs_true_DLS); % half the time flip the field so that the loss is in the inferior field, to show that it isn't a fluke
                    end
                    simObs_true_sigma = [ % slope of psychometric function
                        NaN  NaN  NaN  1.2  1.2  1.3  1.2  NaN  NaN  NaN    % +27
                        NaN  NaN  1.3  1.4  1.3  1.4  1.2  1.0  NaN  NaN    % +21
                        NaN  1.6  1.2  1.5  1.2  1.4  1.2  1.3  1.7  NaN    % +15
                        1.7  1.4  1.6  1.7  1.6  1.1  1.1  1.7  1.8  1.2    % +9
                        1.5  1.2  1.7  1.5  1.0  1.7  1.9   3   1.9  1.4     % +3
                        2.3  1.7  1.4  2.2  1.4  1.4  2.1   3   2.5  2.0     % -3
                        2.2  2.1  2.5  1.7  1.9  1.8  2.2  3.2  2.0  2.6    % -9
                        NaN  1.8  1.8  2.8  2.1  2.0  2.3  1.9  2.4  NaN    % -15
                        NaN  NaN  2.5  2.4  2.1  1.9  2.1  2.9  NaN  NaN    % -21
                        NaN  NaN  NaN  3.0  2.2  2.6  2.4  NaN  NaN  NaN    % -27
                    ]; %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
                    simObs_true_lambda = 0.05; % lapse rate
                    simObs_true_gamma = 0.02; % guess rate
                    simObs_true_F = @(x, mu,sigma, lambda, gamma)(gamma + (1-gamma-lambda) * (1-normcdf(x,mu,sigma)  )); % "probability of anscorrect" psychometric function (v. similar in form to the one assumed by QUEST+, but here for ease of manipulation we shall make all the variables free, including the slope of the psychometric function

                    % define starting test grid: 10 seed points only
                    ids = [ % right eye format (OD)
                        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN   % +27
                        NaN, NaN, NaN,  3 , NaN, NaN,  5 , NaN, NaN, NaN   % +21
                        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN   % +15
                        NaN, NaN, NaN,  1 , NaN, NaN, NaN, NaN,  10, NaN   % +9
                         7 , NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN   % +3
                         9 , NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN   % -3
                        NaN, NaN, NaN, NaN, NaN, NaN,  2 , NaN,  8 , NaN   % -9
                        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN   % -15
                        NaN, NaN, NaN,  4 , NaN, NaN,  6 , NaN, NaN, NaN   % -21
                        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN   % -27
                    ];  %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
                    idx = ~isnan(ids); % exclude null points
                    x = x_deg(idx);

                    % run simulations (note that from here on code
                    % unchanged from example 2)
                    existingTestLocations_XY = nan(length(x), 2);
                    existingTestLocations_QPlus = cell(1, length(x));
                    k = 0;
                    for id = unique(ids(idx))'
                        k = k + 1;

                        % get/set xy location
                        existingTestLocations_XY(k,:) = [x_deg(ids == id) y_deg(ids == id)];
 
                        % normative values (used to initialise the PMF
                        % prior, and we shall also use a jittered version
                        % of these values for the simulated observer)
                        prior_DLS = norm_DLS(ids == id);
                        prior_sigma = norm_sigma(ids == id);

                        % init psychophysical algorithm
                        QP = QuestPlus(QPlusParams.F, QPlusParams.stimDomain, QPlusParams.paramDomain, QPlusParams.respDomain, QPlusParams.stopRule, QPlusParams.stopCriterion, QPlusParams.minNTrials, QPlusParams.maxNTrials);
                        prior = MEDTEG.getExamplePrior(QPlusParams.paramDomain, prior_DLS, prior_sigma*1.96);
                        QP.initialise(prior, QPlusParams.likelihoodsDat);

                        % get sim observer true sensitivity (jittered versions of normative value)
                        obs_DLS = simObs_true_DLS(ids == id);
                        obs_sigma = simObs_true_sigma(ids == id);
                        obs_lambda = simObs_true_lambda;
                        obs_gamma = simObs_true_gamma;

                        % run trials
                        startGuess = QP.getParamEsts('mode'); %#ok starting guess (for debugging)
                        while ~QP.isFinished() % for trial = 1:6
                            dB = QP.getTargetStim();
                            pC = simObs_true_F(dB, obs_DLS, obs_sigma, obs_lambda, obs_gamma); % QPlusParams.F(dB, obs_true_DLS); % get probability of answering correctly
                            anscorrect = rand() < pC;
                            QP.update(dB, anscorrect);
                            totalNTrialsCounter = totalNTrialsCounter + 1;
                        end
                        endGuess = QP.getParamEsts('mode'); %#ok final estimate (for debugging)
                        % fprintf('<x=%i, y=%i>  True Sensitivity: %1.1f dB;  Start Guess = %1.1f dB;  Final Estimate: %1.1f\n', x_deg(ids==id), y_deg(ids==id), obs_DLS, startGuess, endGuess)% <- for debugging. Uncomment to see that Q+ really  does get closer to the true answer over time

                        % store result
                        existingTestLocations_QPlus{k} = QP;
                    end

                    % create MEDTEG object
                    weightParams = MEDTEG.getDefaultCandWeightingParameters();
                    weightParams.weight_type = 'flat'; % 'macula_only', 'prefer_central', 'flat'
                    weightParams.excludeCandidatesInOrAroundBlindspot = true;
                    M = MEDTEG(@MEDTEG.exampleCandWeightingFunc, weightParams, QPlusParams);

                    % create interpolants for simulated-observer data (used below)
                    idx = ~isnan(norm_DLS) & norm_DLS > 1; % exclude untested points and points in physiologic blindspots
                    x = x_deg(idx);
                    y = y_deg(idx);
                    v = simObs_true_DLS(idx);
                    F_simObs_DLS = scatteredInterpolant(x,y,v);
                    v = simObs_true_sigma(idx);
                    F_simObs_sigma = scatteredInterpolant(x,y,v);

                    % similarly, create interpolants for normative data (used below)
                    v = norm_DLS(idx);
                    F_norm_DLS = scatteredInterpolant(x,y,v);

                    % now get and test an additional N points using MEDTEG
                    for i = 1:n
                        [bestCandidateXY, bestCandidateQPlus, ~, candidates, ~, cand_wEDHdeg2] = M.GetNewTestLocation(existingTestLocations_XY, existingTestLocations_QPlus);

                        % init psychophysical algorithm
                        QP = QuestPlus(QPlusParams.F, QPlusParams.stimDomain, QPlusParams.paramDomain, QPlusParams.respDomain, QPlusParams.stopRule, QPlusParams.stopCriterion, QPlusParams.minNTrials, QPlusParams.maxNTrials);
                        prior_DLS = bestCandidateQPlus.getParamEsts('mode');
                        prior_sigma = 5;
                        prior = MEDTEG.getExamplePrior(QPlusParams.paramDomain, prior_DLS, prior_sigma); % preferable to the raw prior computed by MEDTEG
                        QP.initialise(prior, QPlusParams.likelihoodsDat);

                        % get sim observer true sensitivity by interpolation
                        x = bestCandidateXY(1);
                        y = bestCandidateXY(2);
                        obs_DLS = F_simObs_DLS(x, y);
                        obs_sigma = F_simObs_sigma(x, y);
                        obs_lambda = 0.05; % lapse rate
                        obs_gamma = 0.02; % guess rate

                        % run trials
                        startGuess = QP.getParamEsts('mode'); % starting guess (for debugging)
                        while ~QP.isFinished() % for trial = 1:6
                            dB = QP.getTargetStim();
                            pC = simObs_true_F(dB, obs_DLS, obs_sigma, obs_lambda, obs_gamma); % QPlusParams.F(dB, obs_true_DLS); % get probability of answering correctly
                            anscorrect = rand() < pC;
                            QP.update(dB, anscorrect);
                            totalNTrialsCounter = totalNTrialsCounter + 1;
                        end
                        endGuess = QP.getParamEsts('mode'); % final estimate (for debugging)
                        fprintf('<x=%1.0f, y=%1.0f>  True Sensitivity: %1.1f dB;  Start Guess = %1.1f dB;  Final Estimate: %1.1f\n', x, y, obs_DLS, startGuess, endGuess)% <- for debugging. Uncomment to see that Q+ really  does get closer to the true answer over time

                        % store result
                        existingTestLocations_XY(end+1, :) = bestCandidateXY; %#ok
                        existingTestLocations_QPlus{end+1} = QP; %#ok
                    end

                    % for plotting purposes we'll explicitly store the original and new test locations separately
                    newTestLocations_XY                             = existingTestLocations_XY((end-(n-1)):end, :);
                    newTestLocations_QPlus                          = existingTestLocations_QPlus((end-(n-1)):end);
                    existingTestLocations_XY((end-(n-1)):end, :)    = [];
                    existingTestLocations_QPlus((end-(n-1)):end)    = [];                 
                
                case 3 % Seeds, using a priori structural data as weights (NB: "Simulation 3" in paper)
                    % e.g., MEDTEG.runExample(3, struct('structNoiseLevel','low', 'structIsVerical',false));
                    
                    % define normative values (used for setting priors)
                    [x_deg, y_deg, norm_DLS, norm_sigma] = MEDTEG.getExampleNormativeValues();

                    % init simulation params
                    QPlusParams = MEDTEG.getExampleQuestPlusParams();%  get standard Q+ params
                    totalNTrialsCounter = 0;
                    n = 14; % n additional test locations to add (44 + 10 seeds, to bring up to parity with 24-2)

                    % define simulated observer
                    % healthy vision
                    % simObs_true_DLS = [ % right eye format (OD)
                    %     NaN  NaN  NaN  26.8 26.3 25.5 25.1 NaN  NaN  NaN    % +27
                    %     NaN  NaN  29.0 29.2 28.3 27.6 27.8 28.0 NaN  NaN    % +21
                    %     NaN  28.9 29.1 29.8 30.2 29.8 30.9 30.4 28.8 NaN    % +15
                    %     28.3 29.6 31.2 32.6 32.3 31.7 31.2 30.6 30.4 30.3   % +9
                    %     28.9 30.4 32.1 32.8 33.7 33.5 32.2   0  31.3 31.9   % +3
                    %     29.4 31.0 32.6 33.3 33.8 33.8 32.7   0  31.5 31.3   % -3
                    %     29.6 30.4 31.5 33.4 32.9 32.8 32.8 31.1 31.6 31.1   % -9
                    %     NaN  30.5 30.6 31.2 32.1 31.3 30.9 32.4 31.2 NaN    % -15
                    %     NaN  NaN  30.1 29.8 30.9 31.4 31.8 31.2 NaN  NaN    % -21
                    %     NaN  NaN  NaN  29.0 30.1 30.8 30.7 NaN  NaN  NaN    % -27
                    %     ]; %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
                    simObs_true_DLS = [ % right eye format (OD)
                        NaN  NaN  NaN  26.8 26.3 25.5 25.1 NaN  NaN  NaN    % +27
                        NaN  NaN  11.0 29.2 28.3 27.6 27.8 28.0 NaN  NaN    % +21
                        NaN   6.9  6.1  5.8 30.2 29.8 30.9 30.4 28.8 NaN    % +15
                        9.3  4.6  5.2 32.6 32.3 31.7 31.2 30.6 30.4 30.3   % +9
                        9.9  10.4 32.1 32.8 33.7 33.5 32.2   0  31.3 31.9   % +3
                        29.4 31.0 32.6 33.3 33.8 33.8 32.7   0  31.5 31.3   % -3
                        29.6 30.4 31.5 33.4 32.9 32.8 32.8 31.1 31.6 31.1   % -9
                        NaN  30.5 30.6 31.2 32.1 31.3 30.9 32.4 31.2 NaN    % -15
                        NaN  NaN  30.1 29.8 30.9 31.4 31.8 31.2 NaN  NaN    % -21
                        NaN  NaN  NaN  29.0 30.1 30.8 30.7 NaN  NaN  NaN    % -27
                    ]; %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
                    simObs_true_DLS = simObs_true_DLS + randn(size(simObs_true_DLS)).*norm_sigma; % also add some random jitter to all points, just for some added complexity
                    simObs_true_sigma = [ % slope of psychometric function
                        NaN  NaN  NaN  5.8  3.5  4.7  3.9  NaN  NaN  NaN    % +27
                        NaN  NaN  3.6  2.6  4.1  2.6  2.6  2.0  NaN  NaN    % +21
                        NaN  2.9  3.2  2.7  2.2  2.4  2.2  2.3  2.7  NaN    % +15
                        2.7  2.4  1.6  1.7  1.6  2.1  2.1  2.7  2.8  3.2    % +9
                        2.5  2.2  1.7  1.5  2.0  1.7  1.9   3   2.9  2.4     % +3
                        2.3  1.7  1.4  2.2  1.4  1.4  2.1   3   2.5  2.0     % -3
                        2.2  2.1  2.5  1.7  1.9  1.8  2.2  3.2  2.0  2.6    % -9
                        NaN  1.8  1.8  2.8  2.1  2.0  2.3  1.9  2.4  NaN    % -15
                        NaN  NaN  2.5  2.4  2.1  1.9  2.1  2.9  NaN  NaN    % -21
                        NaN  NaN  NaN  3.0  2.2  2.6  2.4  NaN  NaN  NaN    % -27
                    ]; %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
                    simObs_true_lambda = 0.05; % lapse rate
                    simObs_true_gamma = 0.02; % guess rate
                    simObs_true_F = @(x, mu,sigma, lambda, gamma)(gamma + (1-gamma-lambda) * (1-normcdf(x,mu,sigma)  )); % "probability of anscorrect" psychometric function (v. similar in form to the one assumed by QUEST+, but here for ease of manipulation we shall make all the variables free, including the slope of the psychometric function

                    % define starting test grid: 10 seed points only
                    ids = [ % right eye format (OD)
                        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN   % +27
                        NaN, NaN, NaN,  3 , NaN, NaN,  5 , NaN, NaN, NaN   % +21
                        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN   % +15
                        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN,  10, NaN   % +9
                         7 , NaN, NaN, NaN, 1, NaN, NaN, NaN, NaN, NaN   % +3
                         9 , NaN, NaN, NaN, NaN, 2, NaN, NaN, NaN, NaN   % -3
                        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN,  8 , NaN   % -9
                        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN   % -15
                        NaN, NaN, NaN,  4 , NaN, NaN,  6 , NaN, NaN, NaN   % -21
                        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN   % -27
                    ];  %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
                    idx = ~isnan(ids); % exclude null points
                    x = x_deg(idx);

                    % run simulations (note that from here on code
                    % unchanged from example 2)
                    existingTestLocations_XY = nan(length(x), 2);
                    existingTestLocations_QPlus = cell(1, length(x));
                    k = 0;
                    for id = unique(ids(idx))'
                        k = k + 1;

                        % get/set xy location
                        existingTestLocations_XY(k,:) = [x_deg(ids == id) y_deg(ids == id)];
 
                        % normative values (used to initialise the PMF
                        % prior, and we shall also use a jittered version
                        % of these values for the simulated observer)
                        prior_DLS = norm_DLS(ids == id);
                        prior_sigma = norm_sigma(ids == id);

                        % init psychophysical algorithm
                        QP = QuestPlus(QPlusParams.F, QPlusParams.stimDomain, QPlusParams.paramDomain, QPlusParams.respDomain, QPlusParams.stopRule, QPlusParams.stopCriterion, QPlusParams.minNTrials, QPlusParams.maxNTrials);
                        prior = MEDTEG.getExamplePrior(QPlusParams.paramDomain, prior_DLS, prior_sigma*1.96);
                        QP.initialise(prior, QPlusParams.likelihoodsDat);

                        % get sim observer true sensitivity (jittered versions of normative value)
                        obs_DLS = simObs_true_DLS(ids == id);
                        obs_sigma = simObs_true_sigma(ids == id);
                        obs_lambda = simObs_true_lambda;
                        obs_gamma = simObs_true_gamma;

                        % run trials
                        % startGuess = QP.getParamEsts('mode'); %#ok starting guess (for debugging)
                        while ~QP.isFinished() % for trial = 1:6
                            dB = QP.getTargetStim();
                            pC = simObs_true_F(dB, obs_DLS, obs_sigma, obs_lambda, obs_gamma); % QPlusParams.F(dB, obs_true_DLS); % get probability of answering correctly
                            anscorrect = rand() < pC;
                            QP.update(dB, anscorrect);
                            totalNTrialsCounter = totalNTrialsCounter + 1;
                        end
                        % endGuess = QP.getParamEsts('mode'); %#ok final estimate (for debugging)
                        % fprintf('<x=%i, y=%i>  True Sensitivity: %1.1f dB;  Start Guess = %1.1f dB;  Final Estimate: %1.1f\n', x_deg(ids==id), y_deg(ids==id), obs_DLS, startGuess, endGuess)% <- for debugging. Uncomment to see that Q+ really  does get closer to the true answer over time

                        % store result
                        existingTestLocations_QPlus{k} = QP;
                    end

                    % create MEDTEG object
                    weightParams = MEDTEG.getDefaultCandWeightingParameters();
                    weightParams.weight_type = 'flat'; % 'macula_only', 'prefer_central', 'flat'
                    weightParams.excludeCandidatesInOrAroundBlindspot = true;
                    weightParams.noiseLevel = 'low'; % 'none', 'low', 'high', 'extreme'
                    weightParams.isVeridical = true; % false;
                    % overwrite if user requested
                    if (exist('optionalParams','var') && ~isempty(optionalParams))
                        weightParams.noiseLevel = optionalParams.structNoiseLevel;
                        weightParams.isVeridical = optionalParams.structIsVerical;
                    end
                    M = MEDTEG(@MEDTEG.customStructureGuidedCandWeightingFunc, weightParams, QPlusParams);
                    

                    % create interpolants for simulated-observer data (used below)
                    idx = ~isnan(norm_DLS) & norm_DLS > 1; % exclude untested points and points in physiologic blindspots
                    x = x_deg(idx);
                    y = y_deg(idx);
                    v = simObs_true_DLS(idx);
                    F_simObs_DLS = scatteredInterpolant(x,y,v);
                    v = simObs_true_sigma(idx);
                    F_simObs_sigma = scatteredInterpolant(x,y,v);

                    % similarly, create interpolants for normative data (used below)
                    v = norm_DLS(idx);
                    F_norm_DLS = scatteredInterpolant(x,y,v);

                    % now get and test an additional N points using MEDTEG
                    for i = 1:n
                        % get new location to test
                        [bestCandidateXY, bestCandidateQPlus, ~, candidates, ~, cand_wEDHdeg2] = M.GetNewTestLocation(existingTestLocations_XY, existingTestLocations_QPlus);

                        % init psychophysical algorithm
                        QP = QuestPlus(QPlusParams.F, QPlusParams.stimDomain, QPlusParams.paramDomain, QPlusParams.respDomain, QPlusParams.stopRule, QPlusParams.stopCriterion, QPlusParams.minNTrials, QPlusParams.maxNTrials);
                        prior_DLS = bestCandidateQPlus.getParamEsts('mode');
                        prior_sigma = 5;
                        prior = MEDTEG.getExamplePrior(QPlusParams.paramDomain, prior_DLS, prior_sigma); % preferable to the raw prior computed by MEDTEG
                        QP.initialise(prior, QPlusParams.likelihoodsDat);

                        % get sim observer true sensitivity by interpolation
                        x = bestCandidateXY(1);
                        y = bestCandidateXY(2);
                        obs_DLS = F_simObs_DLS(x, y);
                        obs_sigma = F_simObs_sigma(x, y);
                        obs_lambda = 0.05; % lapse rate
                        obs_gamma = 0.02; % guess rate

                        % run trials
                        % startGuess = QP.getParamEsts('mode'); % starting guess (for debugging)
                        while ~QP.isFinished() % for trial = 1:6
                            dB = QP.getTargetStim();
                            pC = simObs_true_F(dB, obs_DLS, obs_sigma, obs_lambda, obs_gamma); % QPlusParams.F(dB, obs_true_DLS); % get probability of answering correctly
                            anscorrect = rand() < pC;
                            QP.update(dB, anscorrect);
                            totalNTrialsCounter = totalNTrialsCounter + 1;
                        end
                        % endGuess = QP.getParamEsts('mode'); % final estimate (for debugging)
                        % fprintf('<x=%1.0f, y=%1.0f>  True Sensitivity: %1.1f dB;  Start Guess = %1.1f dB;  Final Estimate: %1.1f\n', x, y, obs_DLS, startGuess, endGuess)% <- for debugging. Uncomment to see that Q+ really  does get closer to the true answer over time

                        % store result
                        existingTestLocations_XY(end+1, :) = bestCandidateXY; %#ok
                        existingTestLocations_QPlus{end+1} = QP; %#ok
                    end

                    % for plotting purposes we'll explicitly store the original and new test locations separately
                    newTestLocations_XY                             = existingTestLocations_XY((end-(n-1)):end, :);
                    newTestLocations_QPlus                          = existingTestLocations_QPlus((end-(n-1)):end);
                    existingTestLocations_XY((end-(n-1)):end, :)    = [];
                    existingTestLocations_QPlus((end-(n-1)):end)    = [];

                case 4 % Small scotoma plus a complete loss of attention after N locations (NB: "Simulation 4" in paper)

                    % define normative values (used for setting priors)
                    [x_deg, y_deg, norm_DLS, norm_sigma] = MEDTEG.getExampleNormativeValues();

                    % init simulation params
                    QPlusParams = MEDTEG.getExampleQuestPlusParams();%  get standard Q+ params
                    maxNLocs = 52;

                    simObs_true_DLS = [ % right eye format (OD)
                        NaN  NaN  NaN  26.8 26.3 25.5 25.1 NaN  NaN  NaN    % +27
                        NaN  NaN  11.0 29.2 28.3 27.6 27.8 28.0 NaN  NaN    % +21
                        NaN   6.9  6.1  5.8 30.2 29.8 30.9 30.4 28.8 NaN    % +15
                        9.3   4.6  5.2 32.6 32.3 31.7 31.2 30.6 30.4 30.3   % +9
                        9.9  10.4 32.1 32.8 33.7 33.5 32.2   0  31.3 31.9   % +3
                        29.4 31.0 32.6 33.3 33.8 33.8 32.7   0  31.5 31.3   % -3
                        29.6 30.4 31.5 33.4 32.9 32.8 32.8 31.1 31.6 31.1   % -9
                        NaN  30.5 30.6 31.2 32.1 31.3 30.9 32.4 31.2 NaN    % -15
                        NaN  NaN  30.1 29.8 30.9 31.4 31.8 31.2 NaN  NaN    % -21
                        NaN  NaN  NaN  29.0 30.1 30.8 30.7 NaN  NaN  NaN    % -27
                    ]; %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
                    simObs_true_DLS = simObs_true_DLS + randn(size(simObs_true_DLS)).*norm_sigma; % also add some random jitter to all points, just for some added complexity
                    simObs_true_sigma = [ % slope of psychometric function
                        NaN  NaN  NaN  5.8  3.5  4.7  3.9  NaN  NaN  NaN    % +27
                        NaN  NaN  3.6  2.6  4.1  2.6  2.6  2.0  NaN  NaN    % +21
                        NaN  2.9  3.2  2.7  2.2  2.4  2.2  2.3  2.7  NaN    % +15
                        2.7  2.4  1.6  1.7  1.6  2.1  2.1  2.7  2.8  3.2    % +9
                        2.5  2.2  1.7  1.5  2.0  1.7  1.9   3   2.9  2.4     % +3
                        2.3  1.7  1.4  2.2  1.4  1.4  2.1   3   2.5  2.0     % -3
                        2.2  2.1  2.5  1.7  1.9  1.8  2.2  3.2  2.0  2.6    % -9
                        NaN  1.8  1.8  2.8  2.1  2.0  2.3  1.9  2.4  NaN    % -15
                        NaN  NaN  2.5  2.4  2.1  1.9  2.1  2.9  NaN  NaN    % -21
                        NaN  NaN  NaN  3.0  2.2  2.6  2.4  NaN  NaN  NaN    % -27
                    ]; %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
                    simObs_true_lambda = 0.05; % lapse rate
                    simObs_true_gamma = 0.02; % guess rate
                    simObs_true_F = @(x, mu,sigma, lambda, gamma)(gamma + (1-gamma-lambda) * (1-normcdf(x,mu,sigma)  )); % "probability of anscorrect" psychometric function (v. similar in form to the one assumed by QUEST+, but here for ease of manipulation we shall make all the variables free, including the slope of the psychometric function

                    % define starting test grid: 10 seed points only
                    ids = [ % right eye format (OD)
                        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN   % +27
                        NaN, NaN, NaN,  3 , NaN, NaN,  5 , NaN, NaN, NaN   % +21
                        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN   % +15
                        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN,  10, NaN   % +9
                         7 , NaN, NaN, NaN, 1,   NaN, NaN, NaN, NaN, NaN   % +3
                         9 , NaN, NaN, NaN, NaN, 2,   NaN, NaN, NaN, NaN   % -3
                        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN,  8 , NaN   % -9
                        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN   % -15
                        NaN, NaN, NaN,  4 , NaN, NaN,  6 , NaN, NaN, NaN   % -21
                        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN   % -27
                    ];  %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
                    idx = ~isnan(ids); % exclude null points
                    seed_ids = ids(idx);
                    x = x_deg(idx);

                    % for 24-2 randomly shuffle the order of the remaining
                    % 42 points
                    phase2_ids = [ % 24-2 remaining points
                        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN % +27
                        NaN, NaN, NaN, NaN,   1,   1, NaN, NaN, NaN, NaN % +21
                        NaN, NaN,   1,   1,   1,   1,   1,   1, NaN, NaN % +15
                        NaN,   1,   1,   1,   1,   1,   1,   1, NaN, NaN % +9
                         1,    1,   1,   1,   1,   1,   1, NaN,   1, NaN % +3
                         1,    1,   1,   1,   1,   1,   1, NaN,   1, NaN % -3
                        NaN,   1,   1,   1,   1,   1,   1,   1,   1, NaN % -9
                        NaN, NaN,   1,   1,   1,   1,   1,   1,   1, NaN % -15
                        NaN, NaN, NaN,   1,   1,   1,   1, NaN, NaN, NaN % -21
                        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN % -27
                    ];  %-27  -21  -15   -9   -3   +3   +9  +15  +21  +27
                    phase2_idx = find(~isnan(phase2_ids) & isnan(ids)); % exclude null points or points already tested in phase 1
                    phase2_idx=phase2_idx(randperm(length(phase2_idx))); % randomly shuffle order


                    % randomly determine how many points we'll test (before
                    % losing attentiveness)
                    minN = sum(idx(:)); % n in seed grid
                    maxN = maxNLocs-minN;
                    nAdditionalTestLocsToTest = ceil(rand()*maxN);
                    nTestLocsToTest = nAdditionalTestLocsToTest + minN; % +minN since we'll assume that all observers can do the initial 10 seed points
                    

                    % 1. run 24-2 location-by-location in random order
                    testLocations_standard_XY = nan(nTestLocsToTest, 2);
                    testLocations_standard_QPlus = cell(1, nTestLocsToTest);
                    k = 0;
                    for idx = phase2_idx'
                        k = k+ 1;

                        % get/set xy location
                        testLocations_standard_XY(k,:) = [x_deg(idx) y_deg(idx)];
 
                        % normative values (used to initialise the PMF
                        % prior, and we shall also use a jittered version
                        % of these values for the simulated observer)
                        prior_DLS = norm_DLS(idx);
                        prior_sigma = norm_sigma(idx);

                        % init psychophysical algorithm
                        QP = QuestPlus(QPlusParams.F, QPlusParams.stimDomain, QPlusParams.paramDomain, QPlusParams.respDomain, QPlusParams.stopRule, QPlusParams.stopCriterion, QPlusParams.minNTrials, QPlusParams.maxNTrials);
                        prior = MEDTEG.getExamplePrior(QPlusParams.paramDomain, prior_DLS, prior_sigma*1.96);
                        QP.initialise(prior, QPlusParams.likelihoodsDat);

                        % get sim observer true sensitivity (jittered versions of normative value)
                        obs_DLS = simObs_true_DLS(idx);
                        obs_sigma = simObs_true_sigma(idx);
                        obs_lambda = simObs_true_lambda;
                        obs_gamma = simObs_true_gamma;

                        % run trials
                        % startGuess = QP.getParamEsts('mode'); %#ok starting guess (for debugging)
                        while ~QP.isFinished() % for trial = 1:6
                            dB = QP.getTargetStim();
                            pC = simObs_true_F(dB, obs_DLS, obs_sigma, obs_lambda, obs_gamma); % QPlusParams.F(dB, obs_true_DLS); % get probability of answering correctly
                            anscorrect = rand() < pC;
                            QP.update(dB, anscorrect);
                        end
                        % endGuess = QP.getParamEsts('mode'); %#ok final estimate (for debugging)
                        % fprintf('<x=%i, y=%i>  True Sensitivity: %1.1f dB;  Start Guess = %1.1f dB;  Final Estimate: %1.1f\n', x_deg(ids==id), y_deg(ids==id), obs_DLS, startGuess, endGuess)% <- for debugging. Uncomment to see that Q+ really  does get closer to the true answer over time

                        % store result
                        testLocations_standard_QPlus{k} = QP;

                        % check if done
                        if (k >= nTestLocsToTest)
                            break
                        end
                    end


                    % 2. run MEDTEG
                    
                    % 2.1 seed points
                    testLocations_MEDTEG_XY = nan(length(x), 2);
                    testLocations_MEDTEG_QPlus = cell(1, length(x));
                    k = 0;
                    for id = unique(seed_ids)'
                        k = k + 1;

                        % get/set xy location
                        testLocations_MEDTEG_XY(k,:) = [x_deg(ids == id) y_deg(ids == id)];
 
                        % normative values (used to initialise the PMF
                        % prior, and we shall also use a jittered version
                        % of these values for the simulated observer)
                        prior_DLS = norm_DLS(ids == id);
                        prior_sigma = norm_sigma(ids == id);

                        % init psychophysical algorithm
                        QP = QuestPlus(QPlusParams.F, QPlusParams.stimDomain, QPlusParams.paramDomain, QPlusParams.respDomain, QPlusParams.stopRule, QPlusParams.stopCriterion, QPlusParams.minNTrials, QPlusParams.maxNTrials);
                        prior = MEDTEG.getExamplePrior(QPlusParams.paramDomain, prior_DLS, prior_sigma*1.96);
                        QP.initialise(prior, QPlusParams.likelihoodsDat);

                        % get sim observer true sensitivity (jittered versions of normative value)
                        obs_DLS = simObs_true_DLS(ids == id);
                        obs_sigma = simObs_true_sigma(ids == id);
                        obs_lambda = simObs_true_lambda;
                        obs_gamma = simObs_true_gamma;

                        % run trials
                        startGuess = QP.getParamEsts('mode'); %#ok starting guess (for debugging)
                        while ~QP.isFinished() % for trial = 1:6
                            dB = QP.getTargetStim();
                            pC = simObs_true_F(dB, obs_DLS, obs_sigma, obs_lambda, obs_gamma); % QPlusParams.F(dB, obs_true_DLS); % get probability of answering correctly
                            anscorrect = rand() < pC;
                            QP.update(dB, anscorrect);
                        end
                        endGuess = QP.getParamEsts('mode'); %#ok final estimate (for debugging)
                        fprintf('<x=%i, y=%i>  True Sensitivity: %1.1f dB;  Start Guess = %1.1f dB;  Final Estimate: %1.1f\n', x_deg(ids==id), y_deg(ids==id), obs_DLS, startGuess, endGuess)% <- for debugging. Uncomment to see that Q+ really  does get closer to the true answer over time

                        % store result
                        testLocations_MEDTEG_QPlus{k} = QP;
                    end

                    % 2.2 additional MEDTEG mediated points

                    % create MEDTEG object
                    weightParams = MEDTEG.getDefaultCandWeightingParameters();
                    weightParams.weight_type = 'flat'; % 'macula_only', 'prefer_central', 'flat'
                    weightParams.excludeCandidatesInOrAroundBlindspot = true;
                    M = MEDTEG(@MEDTEG.exampleCandWeightingFunc, weightParams, QPlusParams);

                    % create interpolants for simulated-observer data (used below)
                    idx = ~isnan(norm_DLS) & norm_DLS > 1; % exclude untested points and points in physiologic blindspots
                    x = x_deg(idx);
                    y = y_deg(idx);
                    v = simObs_true_DLS(idx);
                    F_simObs_DLS = scatteredInterpolant(x,y,v);
                    v = simObs_true_sigma(idx);
                    F_simObs_sigma = scatteredInterpolant(x,y,v);

                    % similarly, create interpolants for normative data (used below)
                    v = norm_DLS(idx);
                    F_norm_DLS = scatteredInterpolant(x,y,v);

                    % now get and test an additional N points using MEDTEG
                    for i = 1:nAdditionalTestLocsToTest
                        % get new location to test
                        [bestCandidateXY, bestCandidateQPlus] = M.GetNewTestLocation(testLocations_MEDTEG_XY, testLocations_MEDTEG_QPlus);

                        % init psychophysical algorithm
                        QP = QuestPlus(QPlusParams.F, QPlusParams.stimDomain, QPlusParams.paramDomain, QPlusParams.respDomain, QPlusParams.stopRule, QPlusParams.stopCriterion, QPlusParams.minNTrials, QPlusParams.maxNTrials);
                        prior_DLS = bestCandidateQPlus.getParamEsts('mode');
                        prior_sigma = 5;
                        prior = MEDTEG.getExamplePrior(QPlusParams.paramDomain, prior_DLS, prior_sigma); % preferable to the raw prior computed by MEDTEG
                        QP.initialise(prior, QPlusParams.likelihoodsDat);

                        % get sim observer true sensitivity by interpolation
                        x = bestCandidateXY(1);
                        y = bestCandidateXY(2);
                        obs_DLS = F_simObs_DLS(x, y);
                        obs_sigma = F_simObs_sigma(x, y);
                        obs_lambda = 0.05; % lapse rate
                        obs_gamma = 0.02; % guess rate

                        % run trials
                        startGuess = QP.getParamEsts('mode'); % starting guess (for debugging)
                        while ~QP.isFinished() % for trial = 1:6
                            dB = QP.getTargetStim();
                            pC = simObs_true_F(dB, obs_DLS, obs_sigma, obs_lambda, obs_gamma); % QPlusParams.F(dB, obs_true_DLS); % get probability of answering correctly
                            anscorrect = rand() < pC;
                            QP.update(dB, anscorrect);
                        end
                        endGuess = QP.getParamEsts('mode'); % final estimate (for debugging)
                        fprintf('<x=%1.0f, y=%1.0f>  True Sensitivity: %1.1f dB;  Start Guess = %1.1f dB;  Final Estimate: %1.1f\n', x, y, obs_DLS, startGuess, endGuess)% <- for debugging. Uncomment to see that Q+ really  does get closer to the true answer over time

                        % store result
                        testLocations_MEDTEG_XY(end+1, :) = bestCandidateXY; %#ok
                        testLocations_MEDTEG_QPlus{end+1} = QP; %#ok
                    end

                    % package up data and return (without plotting)
                    dat = table('Size',[nTestLocsToTest*2 7],'VariableTypes',["double","double","double","double","double","logical","logical"],'VariableNames',["x_deg", "y_deg", "norm_DLS", "observer_true_DLS", "observer_est_DLS", "isNewlyAddedPoint", "isStandard"]);
                    for i = 1:nTestLocsToTest
                        est_DLS_dB = testLocations_standard_QPlus{i}.getParamEsts('mode');
                        norm_DLS_dB = F_norm_DLS(testLocations_standard_XY(i,1), testLocations_standard_XY(i,2));
                        true_DLS_dB = F_simObs_DLS(testLocations_standard_XY(i,1), testLocations_standard_XY(i,2));
                        dat{i,:} = [testLocations_standard_XY(i,:) norm_DLS_dB true_DLS_dB est_DLS_dB 0 true]; % add row to data table
                    end
                    for i = 1:nTestLocsToTest
                        est_DLS_dB = testLocations_MEDTEG_QPlus{i}.getParamEsts('mode');
                        norm_DLS_dB = F_norm_DLS(testLocations_MEDTEG_XY(i,1), testLocations_MEDTEG_XY(i,2));
                        true_DLS_dB = F_simObs_DLS(testLocations_MEDTEG_XY(i,1), testLocations_MEDTEG_XY(i,2));
                        dat{i+nTestLocsToTest,:} = [testLocations_MEDTEG_XY(i,:) norm_DLS_dB true_DLS_dB est_DLS_dB 0 false]; % add row to data table
                    end
                    return;

                otherwise
                    error('Specified example not recognised.\n\nTo run, type:\n   MEDTEG.runExample(n)\nwhere n is an integer %i..%i\n\nE.g., MEDTEG.runExample(6);', 1, 1);
            end
            
         	% display generic status report (@TODO)
            %M.disp(); % show MEDTEG info
            
            % compare estimates to 'ground truth' (known exactly, since simulated)
            % @TODO

            % All done
            fprintf('Simulation complete\n');

            % Plot results
            figure();
            %
            subplot(2,2,1);
            hold on
            plot(existingTestLocations_XY(:,1), existingTestLocations_XY(:,2), 'o')
            plot([15 15],[-3 3],'kx')
            plot(newTestLocations_XY(:,1), newTestLocations_XY(:,2), 'rs');
            axis square
            plot(xlim(), [0 0], 'k--', [0 0], ylim(), 'k--'); % crosshair
            title('Test Locations')
            
            subplot(2,2,2);
            plot(existingTestLocations_XY(:,1), existingTestLocations_XY(:,2), '.')
            hold on
            plot([15 15],[-3 3],'kx')
            nCandidates = size(candidates,1);
            for i = 1:nCandidates
                txt = sprintf('%1.1f', cand_wEDHdeg2(i));
                text(candidates(i,1), candidates(i,2), txt, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 7);
            end
            plot(bestCandidateXY(:,1), bestCandidateXY(:,2), 'rs', 'MarkerSize',15);
            axis square
            plot(xlim(), [0 0], 'k--', [0 0], ylim(), 'k--'); % crosshair
            title('MEDTEG scores for last added point');

            subplot(2,2,3);
            plot(existingTestLocations_XY(:,1), existingTestLocations_XY(:,2), '.w')
            hold on
            for i = 1:size(existingTestLocations_XY,1)
                DLS_dB = existingTestLocations_QPlus{i}.getParamEsts('mode');
                txt = sprintf('%1.0f', DLS_dB);
                text(existingTestLocations_XY(i,1), existingTestLocations_XY(i,2), txt, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 7);
            end
            for i = 1:size(newTestLocations_XY,1)
                DLS_dB = newTestLocations_QPlus{i}.getParamEsts('mode');
                txt = sprintf('%1.0f', DLS_dB);
                text(newTestLocations_XY(i,1), newTestLocations_XY(i,2), txt, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 7, 'Color','r');
            end
            axis square
            plot(xlim(), [0 0], 'k--', [0 0], ylim(), 'k--'); % crosshair
            title('DLS estimates');

            if (exampleN == 0)
                return; % stop here
            end

            % init data storage
            nRows = length(existingTestLocations_QPlus) + length(newTestLocations_QPlus);
            dat = table('Size',[nRows 6],'VariableTypes',["double","double","double","double","double","logical"],'VariableNames',["x_deg", "y_deg", "norm_DLS", "observer_true_DLS", "observer_est_DLS", "isNewlyAddedPoint"]);
            k = 0;
            %
            subplot(2,2,4);
            plot(existingTestLocations_XY(:,1), existingTestLocations_XY(:,2), '.w')
            hold on
            for i = 1:length(existingTestLocations_QPlus)
                est_DLS_dB = existingTestLocations_QPlus{i}.getParamEsts('mode');
                norm_DLS_dB = F_norm_DLS(existingTestLocations_XY(i,1), existingTestLocations_XY(i,2));
                true_DLS_dB = F_simObs_DLS(existingTestLocations_XY(i,1), existingTestLocations_XY(i,2));
                txt = sprintf('%1.0f', est_DLS_dB - true_DLS_dB);
                if ~ismember(existingTestLocations_XY(i,:), [15 3; 15 -3],'rows') % skip blindspot
                    text(existingTestLocations_XY(i,1), existingTestLocations_XY(i,2), txt, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 7);
                end
                k = k + 1;
                dat{k,:} = [existingTestLocations_XY(i,:) norm_DLS_dB true_DLS_dB est_DLS_dB 0]; % add row to data table
            end
            for i = 1:length(newTestLocations_QPlus)
                est_DLS_dB = newTestLocations_QPlus{i}.getParamEsts('mode');
                norm_DLS_dB = F_norm_DLS(newTestLocations_XY(i,1), newTestLocations_XY(i,2));
                true_DLS_dB = F_simObs_DLS(newTestLocations_XY(i,1), newTestLocations_XY(i,2));
                txt = sprintf('%1.0f', est_DLS_dB - true_DLS_dB);
                text(newTestLocations_XY(i,1), newTestLocations_XY(i,2), txt, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize', 7, 'Color','r');
                k = k + 1;
                dat{k,:} = [newTestLocations_XY(i,:) norm_DLS_dB true_DLS_dB est_DLS_dB 1];
            end
            axis square
            plot(xlim(), [0 0], 'k--', [0 0], ylim(), 'k--'); % crosshair
            title('DLS residual error');
        end

    end
  
end