clear; clc;

%% 1) USER PARAMETERS
k_orig    = 10;            % initial tumor region is nodes 1..k_orig
L         = 30;            % fixed right boundary (>= k_orig)
c_lambda  = 0.1;           % λ_i = c_lambda * i^beta inside-region
alpha_ext = 1.0;           % external rate outside tumor
N         = 1000;          % Monte Carlo samples per beta

% sweep beta from 0.1 to 5 in 0.1 increments
betas     = 0.1:0.1:5.0;
numB      = numel(betas);

rng('shuffle');

%% 2) Preallocate
avgHittingTime = zeros(1, numB);

%% 3) Loop over each beta
for iB = 1:numB
    beta = betas(iB);
    totalTime = 0;
    
    % Monte‐Carlo trials
    for trial = 1:N
        k_curr = k_orig;    % dynamic tumor boundary
        pos    = k_curr+1;  % start just outside the tumor
        t      = 0;
        
        while pos > 0
            % decide rates at current pos
            if pos > k_curr
                % outside tumor
                rateL = alpha_ext;
                rateR = alpha_ext;
            elseif pos == k_curr
                % right at the front of tumor
                rate_i = c_lambda * pos^beta;
                rateL  = rate_i;
                rateR  = rate_i;
            else
                % strictly inside tumor (pos < k_curr)
                rate_i = c_lambda * pos^beta;
                rateL  = rate_i;
                rateR  = rate_i;
            end
            
            % total jump rate and waiting time
            rate_tot = rateL + rateR;
            tau      = -log(rand) / rate_tot;
            t        = t + tau;
            
            % choose direction
            if pos == L
                new_pos = pos - 1;    % forced reflection
            else
                if rand < rateL/rate_tot
                    new_pos = pos - 1;
                else
                    new_pos = pos + 1;
                end
            end
            
            % **check for tumor‐kill event**:
            % if we stepped from k_curr down into k_curr-1,
            % that means the cell at k_curr was eradicated.
            if pos == k_curr && new_pos == k_curr-1
                k_curr = k_curr - 1;
            end
            
            pos = new_pos;
        end
        
        totalTime = totalTime + t;
    end
    
    avgHittingTime(iB) = totalTime / N;
    fprintf('beta = %.2f   ⟨T_hit⟩ = %.4f\n', beta, avgHittingTime(iB));
end

%% 4) Plot ⟨Hitting Time⟩ vs beta
figure('Color','w');
plot(betas, avgHittingTime, 'o-', 'LineWidth',1.5, 'MarkerSize',6);
grid on;
xlabel('β');
ylabel('Average hitting time ⟨T_{hit}⟩');
title(sprintf('⟨T_{hit}⟩ vs β (k_0=%d, L=%d, c=%.2g, α_{ext}=%.2g, N=%d)', ...
              k_orig, L, c_lambda, alpha_ext, N));
pause;
