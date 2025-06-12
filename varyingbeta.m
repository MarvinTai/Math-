clear; clc;

%% 1) USER PARAMETERS
k         = 10;           % tumor region is nodes 1..k
L         = 30;           % fixed right boundary (>= k)
c_lambda  = 0.1;          % λ_i = c_lambda * i^beta inside-region
alpha_ext = 1.0;          % constant rate for left and right when pos > k and pos < L
N         = 1000;         % Monte Carlo samples per beta

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
    
    for trial = 1:N
        pos = k;
        t   = 0;
        
        while pos > 0
            if pos <= k
                % inside tumor: state‐dependent rate λ_i = c_lambda * i^beta
                rate_i    = c_lambda * pos^beta;
                rateLeft  = rate_i;
                rateRight = rate_i;
                
            elseif pos < L
                % outside tumor: constant α_ext for both directions
                rateLeft  = alpha_ext;
                rateRight = alpha_ext;
                
            else
                % pos == L (right‐boundary reflection)
                rateLeft  = alpha_ext;
                rateRight = 0;
            end
            
            alpha_total = rateLeft + rateRight;
            tau         = -log(rand) / alpha_total;
            t           = t + tau;
            
            % decide step
            if pos == L
                pos = pos - 1;
            else
                if rand < rateLeft/alpha_total
                    pos = pos - 1;
                else
                    pos = pos + 1;
                end
            end
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
xlabel('beta');
ylabel('<T_hit> from node 10 to 0');
title(sprintf('⟨T_{hit}⟩ vs \\beta (k=%d, L=%d, c=%.3g, \\alpha_{ext}=%.2g, N=%d)', ...
              k, L, c_lambda, alpha_ext, N));

pause;