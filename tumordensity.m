% ===============================================================
%  Compute average hitting‐time to 0 from node 10 as L → ∞
%  continuous‐time random walk with state‐dependent rates inside [1..k]
% ===============================================================

clear; clc;

%% 1) USER PARAMETERS
k        = 10;           % tumor region is nodes 1..k
% Choose a base constant for inside‐region rates: lambda_i = c_lambda * i
c_lambda = 0.1;          % example: adjust so that rates are reasonable
% Outside region, equal rate for left and right:
alpha_ext = 1.0;         % rate for left and right when pos > k and pos < L

d        = 1.0;          % (unused now for region 1..k; kept if you want baseline)
N        = 1000;         % number of Monte‐Carlo samples per L
beta    = 0.2

% The list of L‐values to test:
Ls = 1:30;               % e.g., 1 through 30

rng('shuffle');          % different RNG seed each run

%% 2) Preallocate storage for the averaged hitting times
numL = numel(Ls);
avgHittingTime = zeros(size(Ls));

%% 3) Loop over each L and perform N simulations
for idxL = 1:numL
    L = Ls(idxL);
    
    % Sanity check: if L < k, cannot start at node 10 (k=10)
    if L < k
        warning('L = %d is smaller than k = %d. Skipping this L.', L, k);
        avgHittingTime(idxL) = NaN;
        continue;
    end
    
    totalTimeSum = 0;  % accumulate hitting times for this L
    
    % Monte‐Carlo loop
    for trial = 1:N
        pos = k;      % start at node 10 (i.e., pos = k)
        t   = 0;      % current time
        
        % Simulate until pos hits 0
        while pos > 0
            % === COMPUTE RATES λ_left, λ_right AT CURRENT pos ===
            if pos <= k
                % inside tumor region: equal probabilities but state‐dependent rates
                % λ_i = c_lambda * pos^beta
                rate_i = c_lambda * pos^beta;
                rateLeft  = rate_i;
                rateRight = rate_i;
            elseif pos < L
                % outside tumor but not at boundary: equal rate alpha_ext
                rateLeft  = alpha_ext;
                rateRight = alpha_ext;
            else
                % pos == L (reflecting at right boundary): forced left at rate alpha_ext
                rateLeft  = alpha_ext;
                rateRight = 0;
            end
            
            alpha = rateLeft + rateRight;   % total escape‐rate
            
            % === SAMPLE WAITING TIME τ ~ Exp(alpha) ===
            r1  = rand();
            tau = -log(r1) / alpha;
            t   = t + tau;
            
            % === MOVE LEFT or RIGHT ===
            if pos == L
                % forced reflection: step left
                pos = pos - 1;
            else
                r2 = rand();
                if r2 < rateLeft / alpha
                    pos = pos - 1;   % step left
                else
                    pos = pos + 1;   % step right
                end
            end
            % loop until pos == 0
        end
        
        % at this point pos == 0; record hitting time
        totalTimeSum = totalTimeSum + t;
    end
    
    % Compute average over N trials
    avgHittingTime(idxL) = totalTimeSum / N;
    
    fprintf('L = %6d   ⟨T_hit⟩ ≈ %.4f  (over %d runs)\n', ...
            L, avgHittingTime(idxL), N);
end

%% 4) PLOT ⟨Hitting Time⟩ vs L
figure('Color','w');
loglog(Ls, avgHittingTime, 'o-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('L (right boundary)');
ylabel('Average hitting time ⟨T_{hit}⟩ from node 10 to 0');
title(sprintf('⟨T_{hit}⟩ vs L   (k=%d, λ_i=c·i with c=%.3g, α_{ext}=%.3g, N=%d)', ...
              k, c_lambda, alpha_ext, N));

% Optionally annotate each point with its numerical value
for i = 1:numL
    if ~isnan(avgHittingTime(i))
        text( Ls(i)*1.1, avgHittingTime(i), ...
              sprintf('%.1f', avgHittingTime(i)), ...
              'FontSize', 8 );
    end
end

pause;
