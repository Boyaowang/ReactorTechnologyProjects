function yPellet = solvePellet(w_guess,T_guess,ptot_guess,uz_guess)
constant;
% %%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%
% wCH40 = FRACin(1) .*ones(mpart,1);
% wCO0 = FRACin(2) .*ones(mpart,1);
% wCO20 = FRACin(3) .*ones(mpart,1);
% wH20 = FRACin(4) .*ones(mpart,1);
% wH2O0 = FRACin(5) .*ones(mpart,1);
% wN2 = FRACin(6) .*ones(mpart,1);
% T0 = Tin*ones(mpart,1);

wCH40 = w_guess(1) .*ones(mpart,1);
wCO0 = w_guess(2) .*ones(mpart,1);
wCO20 = w_guess(3) .*ones(mpart,1);
wH20 = w_guess(4) .*ones(mpart,1);
wH2O0 = w_guess(5) .*ones(mpart,1);
wN2 = w_guess(6) .*ones(mpart,1);
T0 = T_guess*ones(mpart,1);


%initial guess vectors:
init = [wCH40; wCO0; wCO20; wH20; wH2O0; T0];
%parameters vectors:
par = [ptot_guess; uz_guess; T_guess; wCH40(1); wCO0(1); wCO20(1); wH20(1); wH2O0(1)];
%%%%%%%%%%%%%%%%%%%%% non-linear system solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create ODE solver:
options = optimoptions('fsolve','Display','none','PlotFcn',@optimplotfirstorderopt);

func = @calcPellet;
yPellet = fsolve(@(init)calcPellet(init,par), init);