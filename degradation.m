% ---------- SSA for  A --k--> ∅  (Gillespie direct method) ----------
clear; clc;

% --- user parameters --------------------------------------------------
n0    = 40;      % initial molecule count  A(0)
k     = 0.15;    % rate constant
Tmax  = 100;     % stop if t > Tmax    (set [] to run until extinction)
rng('shuffle')           % seed for reproducibility
% ---------------------------------------------------------------------

% storage (pre-allocate a bit for speed, grow if needed)
cap     = 2*n0 + 100;          % rough upper bound on number of jumps
t       = zeros(cap,1);        % jump times  T0,T1,...
Acount  = zeros(cap,1);        % molecule counts A(Tn)
t(1)      = 0;
Acount(1) = n0;

% main loop ------------------------------------------------------------
idx = 1;                       % current event index
while Acount(idx) > 0
    % current state
    A  = Acount(idx);
    dt = -log(rand)/(k*A);     % exponential waiting time  τ
    
    newTime = t(idx) + dt;
    
    % stop if Tmax set and exceeded
    if ~isempty(Tmax) && newTime > Tmax
        break
    end
    
    % record next jump
    idx = idx + 1;
    if idx > cap                     % enlarge arrays if necessary
        cap      = cap*2;
        t        = [t;       zeros(cap-idx+1,1)];
        Acount   = [Acount;  zeros(cap-idx+1,1)];
    end
    
    t(idx)      = newTime;
    Acount(idx) = A - 1;       % reaction consumes 1 molecule
end

% trim unused tail
t       = t(1:idx);
Acount  = Acount(1:idx);

% plotting -------------------------------------------------------------
figure('Color','w');
stairs(t, Acount, 'LineWidth', 2);
xlabel('Time  t');
ylabel('Molecule count  A(t)');
title('SSA trajectory for  A \rightarrow \emptyset');
grid on;

pause;  % waits for user to press Enter in terminal

