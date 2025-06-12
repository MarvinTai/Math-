clear; clc;

%% 1) USER PARAMETERS
k_orig         = 10;            % initial tumor region is nodes 1..k_orig
L              = 30;            % fixed right boundary (>= k_orig)
c_lambda       = 0.1;           % λ_i = c_lambda * i^beta inside-region
alpha_ext      = 1.0;           % external rate outside tumor
N              = 1000;          % Monte Carlo samples per beta

% sweep beta from 0.1 to 5 in 0.1 increments
betas          = 0.1:0.1:5.0;
numB           = numel(betas);

% Regrowth control: number of consecutive jumps outside needed to regrow 1 node
regrowth_jumps = 3;     

rng('shuffle');

%% 2) Preallocate
avgHittingTime = zeros(1, numB);

%% 3) Loop over each beta
for iB = 1:numB
    beta = betas(iB);
    totalTime = 0;
    
    for trial = 1:N
        % initialize dynamic tumor boundary & counters
        k_curr       = k_orig;    
        pos          = k_curr + 1;  
        outside_count= 0;           
        t            = 0;
        
        while pos > 0
            % 3.1) Decide jump rates
            if pos > k_curr
                rateL = alpha_ext;  
                rateR = alpha_ext;
            else
                rate_i = c_lambda * pos^beta;
                rateL  = rate_i;
                rateR  = rate_i;
            end
            
            rate_tot = rateL + rateR;
            tau      = -log(rand) / rate_tot;
            t        = t + tau;
            
            % 3.2) Choose direction
            if pos == L
                new_pos = pos - 1;
            else
                if rand < rateL/rate_tot
                    new_pos = pos - 1;
                else
                    new_pos = pos + 1;
                end
            end
            
            % 3.3) Check for kill event (crossing from k_curr → k_curr-1)
            if pos == k_curr && new_pos == k_curr - 1
                k_curr = k_curr - 1;
                % Reset regrowth counter since we're back at boundary 
                outside_count = 0;
            end
            
            % 3.4) Check for regrowth condition
            if new_pos > k_curr + 1
                % immune cell is “off-guard,” increment counter
                outside_count = outside_count + 1;
                if outside_count >= regrowth_jumps
                    k_curr = k_curr + 1;       % tumor grows outward
                    outside_count = 0;         % reset counter
                end
            else
                % immune back near boundary, reset counter
                outside_count = 0;
            end
            
            pos = new_pos;
        end
        
        totalTime = totalTime + t;
    end
    
    avgHittingTime(iB) = totalTime / N;
    fprintf('beta = %.2f   ⟨T_hit⟩ = %.4f   final_k = %d\n', ...
            beta, avgHittingTime(iB), k_curr);
end

%% 4) Plot ⟨Hitting Time⟩ vs beta
figure('Color','w');
plot(betas, avgHittingTime, 'o-', 'LineWidth',1.5, 'MarkerSize',6);
grid on;
xlabel('\beta');
ylabel('Average hitting time ⟨T_{hit}⟩');
title(sprintf('⟨T_{hit}⟩ vs \\beta (k_0=%d→dynamic, L=%d, c=%.2g, α_{ext}=%.2g, regrowth_jumps=%d, N=%d)', ...
              k_orig, L, c_lambda, alpha_ext, regrowth_jumps, N));
pause;
