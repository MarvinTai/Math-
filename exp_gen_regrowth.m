% ===============================================================
%  Final tumour-size distribution after T_max = 50    (150 runs)
%  Immune CTMC + regrowth with Exp-distributed *rate* λ_G
% ===============================================================

clear; clc;  rng('shuffle');

%% USER-TUNEABLE PARAMETERS
k0        = 10;        % initial tumour width (cells at sites 1..k0)
c_lambda  = 0.1;       % prefactor inside tumour: λ_i = c*i^β
beta      = 0.2;       % exponent β
alpha_ext = 1.0;       % symmetric rate outside tumour
mu_lambda = 0.30;      % mean of Exp-distribution for regrowth *rate*
T_max     = 50;        % simulation horizon
Nsim      = 150;       % how many runs
L         = 30;        % right wall (reflecting); choose large enough

%% STORAGE
k_final = zeros(Nsim,1);          % tumour size after 50 units in each run

%% ---------------  MONTE-CARLO SIMULATIONS  -------------------
for n = 1:Nsim
    k   = k0;                     % current tumour width
    pos = k + 1;                  % immune cell starts just outside
    t   = 0;                      % master clock
    t_reg_next = inf;             % no regrowth scheduled yet

    while t < T_max
        % -- (re)start or cancel regrowth clock depending on distance --
        far_enough = (pos - k >= 3) && (k < L);
        if far_enough && isinf(t_reg_next)
            lambdaG     = -mu_lambda * log(rand());    %  Exp(mean = mu_lambda)
            tauG        = -log(rand()) / lambdaG;      %  waiting-time
            t_reg_next  = t + tauG;
        elseif ~far_enough && ~isinf(t_reg_next)
            t_reg_next  = inf;                         % cancel pending regrowth
        end

        % -- immune-cell move rates --
        if pos <= k                       % inside tumour
            rate_i = c_lambda * pos^beta;
            rL = rate_i;  rR = rate_i;
        elseif pos < L                   % outside, not at wall
            rL = alpha_ext;  rR = alpha_ext;
        else                             % pos == L  (reflecting wall)
            rL = alpha_ext;  rR = 0;
        end
        R_move = rL + rR;

        % -- draw next move time --
        tau_move     = -log(rand()) / R_move;
        t_move_next  = t + tau_move;

        % -- What happens first? (move, regrowth, or horizon) --
        [t_next, tag] = min([t_move_next, t_reg_next, T_max]);  % tag: 1=move,2=reg,3=horizon
        t = t_next;

        if tag == 3                       % horizon reached: stop run
            break
        elseif tag == 2                   % ----- regrowth fires -----
            k          = k + 1;
            t_reg_next = inf;             % one-shot; will restart if eligible
            continue                      % go to next while-iteration
        end

        % ----- otherwise, a move fires -----
        u = rand()*R_move;
        if u < rL                         % step left
            pos = pos - 1;
            if pos == 0                   % tumour eliminated early
                k = 0;  break
            end
        else                              % step right
            pos = pos + 1;
        end
    end

    k_final(n) = k;                       % record outcome of this run
end

%% ---------------  PLOT: histogram of k_final  -----------------
minK = min(k_final);
maxK = max(k_final);
centers       = minK:maxK;                % integer bar centres
[counts, ~]   = hist(k_final, centers);   % legacy-compatible histogram

figure('Color','w');
bar(centers, counts, 0.8);
grid on;
xlabel('Final tumour size  k_{final}');
ylabel('Number of simulations');
title(sprintf(['Distribution after T_{max}=%.0f  (N=%d runs)\n', ...
               'k_0=%d,  μ_{λ}=%.2g,  inside λ_i=c·i^{β}, c=%.2g, β=%.2g, ', ...
               'α_{ext}=%.2g' ], ...
               T_max, Nsim, k0, mu_lambda, c_lambda, beta, alpha_ext));

pause;