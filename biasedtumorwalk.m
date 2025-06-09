% =============================================================
%  One‐particle continuous‐time random walk (biased on 1..k,
%   unbiased on k+1..L−1, reflecting at L, absorbing at 0)
%  -------------------------------------------------------------
clear; clc;

% ==== 1) USER PARAMETERS =====
L   = 30;         % rightmost node (reflecting)
k   = 10;         % last node in the biased region
P_L = 0.4;        % probability of stepping left in biased region
d   = 1.0;        % base jump‐rate scale (so total escape‐rate is always d)
                  % (if you want a different total rate at each i, modify below)
rng('shuffle');   % or comment out for reproducibility
% =============================

%  Initial condition: one particle at node k+1 at time t=0
pos = k + 1;
t   = 0;

%  Storage arrays for plotting
times     = t;       % record jump times:  T0, T1, …
positions = pos;     % record positions:   X0, X1, …

while (pos > 0)
    % =========================================================
    %  2) DETERMINE THE TWO “REACTION RATES” AT CURRENT pos
    % =========================================================
    if pos <= k
        % Biased region (1..k):  left‐hop rate = d*P_L;  right‐hop rate = d*(1-P_L)
        rateLeft  = d * P_L;
        rateRight = d * (1 - P_L);
    elseif pos < L
        % Unbiased region (k+1..L-1):  left‐hop = d*0.5;  right‐hop = d*0.5
        rateLeft  = d * 0.5;
        rateRight = d * 0.5;
    else
        % pos == L  → reflecting right edge ⇒ forced left jump, but STILL wait first
        rateLeft  = d;
        rateRight = 0;
    end

    % total “escape‐rate” from this node:
    alpha = rateLeft + rateRight;
    % (in your setup α == d at every position, but shown here for clarity)

    % =========================================================
    %  3) SAMPLE WAITING TIME τ ~ Exp(alpha)  (                                )
    %     (inverse transform: τ = –(1/alpha)*log(r1) )  <— FOR EVERY NODE, even L
    % =========================================================
    r1  = rand;
    tau = -log(r1) / alpha;
    t   = t + tau;  

    % =========================================================
    %  4) NOW DECIDE: LEFT vs RIGHT  (since α>0, we always have a tau)
    % =========================================================
    if pos == L
        % At the right boundary, reflect immediately: forced left
        pos = pos - 1;
    else
        % Either in biased (1..k) or unbiased (k+1..L-1) region:
        r2 = rand;
        if (r2 < rateLeft / alpha)
            pos = pos - 1;   % step left
        else
            pos = pos + 1;   % step right
        end
    end

    % =========================================================
    %  5) STORE RESULTS
    % =========================================================
    times     = [times;     t];
    positions = [positions; pos];
end

% =============================================================
%  6) PLOT THE STOCHASTIC TRAJECTORY AS A STAIRCASE
% =============================================================
figure('Color','w', 'Name','CTMC random walk');
stairs(times, positions, 'LineWidth', 2);
xlabel('time  t');
ylabel('node  i');
title('CTMC: 1‐particle random walk with position‐dependent rates');
grid on;

pause;   % keep the figure window open until you press a key
