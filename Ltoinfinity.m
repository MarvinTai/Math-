% ===============================================================
%  Compute average hitting‐time to 0 from node 10 as L → ∞
%  using a continuous‐time random walk (biased on 1..k, unbiased on k+1..L−1,
%  reflecting at L, absorbing at 0). We fix k = 10 and start at pos = 10.
%
%  For each L ∈ {10, 100, 1000, 10000, 100000}, we simulate N trials
%  and record the time T at which the walk first hits node 0. Then we
%  average those T’s and plot 〈T〉 vs L.
% ===============================================================

clear; clc;

%% 1) USER PARAMETERS
k   = 10;           % biased region is nodes 1..k
P_L = 0.6;          % left‐jump probability in biased region
d   = 1.0;          % base jump‐rate (so α(i)=d for all i< L)
N   = 10000;          % number of Monte‐Carlo samples per L

% The list of L‐values to test:
Ls = 1:30;

rng('shuffle');     % different RNG seed each run

%% 2) Preallocate storage for the averaged hitting times
numL = numel(Ls);
avgHittingTime = zeros(size(Ls));

%% 3) Loop over each L and perform N simulations
for idxL = 1:numL
    L = Ls(idxL);
    
    % Sanity check: if L < k, we cannot start at node = 10, so skip.
    if L < k
        warning('L = %d is smaller than k = %d. Skipping this L.', L, k);
        avgHittingTime(idxL) = NaN;
        continue;
    end
    
    totalTimeSum = 0;  % accumulate hitting times for this L
    
    % Monte‐Carlo loop
    for trial = 1:N
        pos = 10;      % start at node 10
        t   = 0;       % current time
        
        % Simulate until pos hits 0
        while pos > 0
            % === 3.1) COMPUTE RATES λ_left, λ_right AT CURRENT pos ===
            if pos <= k
                % biased region: left‐rate = d*P_L, right‐rate = d*(1−P_L)
                rateLeft  = d * P_L;
                rateRight = d * (1 - P_L);
            elseif pos < L
                % unbiased region: left‐rate = d*0.5, right‐rate = d*0.5
                rateLeft  = d * 0.5;
                rateRight = d * 0.5;
            else
                % pos == L (reflecting at right boundary): forced left
                rateLeft  = d;
                rateRight = 0;
            end
            
            alpha = rateLeft + rateRight;   % total escape‐rate (always = d)
            
            % === 3.2) SAMPLE WAITING TIME τ ~ Exp(alpha) ===
            r1  = rand();
            tau = -log(r1) / alpha;
            t   = t + tau;
            
            % === 3.3) MOVE LEFT or RIGHT ===
            if pos == L
                % forced reflection: step left
                pos = pos - 1;
            else
                r2 = rand();
                if (r2 < rateLeft / alpha)
                    pos = pos - 1;   % step left
                else
                    pos = pos + 1;   % step right
                end
            end
            % loop back until pos == 0
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
title(sprintf('⟨T_{hit}⟩ vs L   (k=%d, P_L=%.2f, N=%d)', k, P_L, N));

% Optionally annotate each point with its numerical value
for i = 1:numL
    text( Ls(i)*1.1, avgHittingTime(i), ...
          sprintf('%.1f', avgHittingTime(i)), ...
          'FontSize', 8 );
end

pause;