%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Advanced Derivative: PS 9 %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Lionel BRODARD, Antoine MARCHAL & Valentin TISSOT-DAGUETTE %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc; warning("off")
%% Parameters
K       = 0.805;
k       = 0.15;
sig_r = 0.01;
theta   = 0.05;
r0      = 0.042;
T1      = 0.25;
T2      = 5+T1;

%% Analytical Solution

B = (1 - exp(-k*(T2-T1)))/k;

sig_tilde = sig_r * sqrt((1-exp(-2*k*T1))/(2*k)) * B;

h = 1/sig_tilde * log(ZC(T2, r0, k, theta, sig_r)...
    /(ZC(T1, r0, k, theta, sig_r)*K)) + sig_tilde/2;

ZBP = K * ZC(T1,r0,k,theta,sig_r) * normcdf(-h+sig_tilde)...
      - ZC(T2,r0,k,theta,sig_r) * normcdf(-h);

%% Finite difference

% Space domain truncation
r_min = 0;

% Derivation of r_max (please refer to the report for further explanations)

[~,A,B] = ZC(T2-T1,r0,k,theta,sig_r);

% Short rate at which the ZCB price equals the strike price at T_1
r_star = log(A/K)/B 

% Choose a tolerance level
epsilon = 1e-12;

% Standard deviation and mean of the short rate in the Vasicek model
sig = @(tau) sqrt(sig_r^2/(2*k) * (1 - exp(- 2 * k *tau )));

mu = @(tau) r0 * exp(-k*tau) + theta * (1 - exp(-k*tau));

% Function fo maximize
g = @(tau) theta + (r_star - theta - norminv(epsilon)*sig(tau))* exp(k*tau);

figure 
fplot(g,[0,T1],'Linewidth',1.5)
xlabel("Time to maturity \tau")
ylabel("r_{max}(\tau,\epsilon)")

% This suggests to take r_max roughly equal to 0.08
r_max= 0.08;

%Initial boundary condition
initial_cond = @(r) max(0,K-ZC(T2-T1,r,k,theta,sig_r));

%Left boundary condition
bc_left = @(tau) 0;

%Right boundary condition
bc_right = @(tau) K * ZC(tau,r_max,k,theta,sig_r)...
         - ZC(tau + T2-T1,r_max,k,theta,sig_r);

% Rate of convergence with respect to the space grid

% Array of number of steps in the space grid
N_rs = 10:5:150; 

% Smallest space interval
min_dr = (r_max-r_min)./max(N_rs);

% Number of equal length intervals in the time grid
% (chosen s.t. the explicit method is stable)
N_t = ceil((sig_r^2 + r_max*(min_dr)^2)*T1/(min_dr)^2);

fprintf("\nRequired number of time steps for N_r = 150: %d\n",N_t)

ZBP_expl = zeros(length(N_rs),1); ZBP_impl = ZBP_expl; ZBP_CN = ZBP_expl;

for i=1:length(N_rs)
    
    % Implicit method: backward Euler finite time difference (nu=1)
    [V_back,r_grid,t_grid]=FiniteDiff(sig_r,k,theta,bc_left,bc_right,...
                           initial_cond,r_max,r_min,N_rs(i),T1,N_t,1);
    
    % Explicit method: forward Euler finite time difference (nu=0)
    [V_for,~,~]=FiniteDiff(sig_r,k,theta,bc_left,bc_right,...
                          initial_cond,r_max,r_min,N_rs(i),T1,N_t,0);
    
    % Crank-Nicholson method (nu = 1/2)
    [V_CN,~,~]=FiniteDiff(sig_r,k,theta,bc_left,bc_right,...
                         initial_cond,r_max,r_min,N_rs(i),T1,N_t,0.5);
    
    % interpolation of the value corresponding to r0
    idx = find(r_grid<r0, 1, 'last');
    
    ZBP_impl(i) = interp1([r_grid(idx),r_grid(idx+1)],...
                        [V_back(idx,end),V_back(idx+1,end)],r0);
                    
    ZBP_expl(i) = interp1([r_grid(idx),r_grid(idx+1)],...
                       [V_for(idx,end),V_for(idx+1,end)],r0);
                   
    ZBP_CN(i) = interp1([r_grid(idx),r_grid(idx+1)],...
                      [V_CN(idx,end),V_CN(idx+1,end)],r0);
    
end

%% Plot 1 analyse of convergence with respect to the space grid
figure

plot(N_rs, ZBP*ones( length(N_rs),1), 'k','LineWidth',2)
hold on
plot(N_rs, ZBP*ones( length(N_rs),1)*1.01, 'k--','LineWidth',2)
plot(N_rs, ZBP*ones( length(N_rs),1)*0.99, 'k--','LineWidth',2)
plot(N_rs, ZBP_CN, '.-','LineWidth',2,'MarkerSize',10)
plot(N_rs, ZBP_expl, '.-','LineWidth',2,'MarkerSize',10)
plot(N_rs, ZBP_impl, '.--','LineWidth',2,'MarkerSize',10)
legend('analytical', 'band +1%', 'band -1%','Crank-Nicholson',...
       'explicit','implicit' )
xlim([min(N_rs),max(N_rs)])
xlabel('N_r')
ylabel('ZBP(0,T1,T2,K)')


%% Rate of convergence in time

% Number of steps in space and the corresponding step size
N_r = 40; dr = (r_max-r_min)/N_r;

% minimum value for N_t to get stability
N_t_min = ceil((sig_r^2 + r_max*dr^2)*T1/dr^2);

fprintf("\nRequired number of time steps for N_r = 40: %d\n",N_t_min)

% Array of N_t's (start slightly outside the stability region)
N_ts = (N_t_min-3):1:(2*N_t_min);
    
ZBP_expl = zeros(length(N_ts),1); ZBP_impl = ZBP_expl; ZBP_CN = ZBP_expl;

for i=1:length(N_ts)

    % Explicit method: backward finite time difference (\nu=1)
    [V_back,~,~]=FiniteDiff(sig_r,k,theta,bc_left,bc_right,...
                            initial_cond,r_max,r_min,N_r,T1,N_ts(i),1);
    % Implicit method: forward finite time difference (\nu=0)
    [V_for,~,~]=FiniteDiff(sig_r,k,theta,bc_left,bc_right,...
                           initial_cond,r_max,r_min,N_r,T1,N_ts(i),0);
    % Crank-Nicholson method
    [V_CN,r_grid,t_grid]=FiniteDiff(sig_r,k,theta,bc_left,bc_right,...
                                    initial_cond,r_max,r_min,N_r  ,...
                                    T1,N_ts(i),0.5);
    
    % interpolation of the value corresponding to r0
    idx = find(r_grid<r0, 1, 'last');
    
    ZBP_impl(i) = interp1([r_grid(idx),r_grid(idx+1)],...
                        [V_back(idx,end),V_back(idx+1,end)],r0);
    
    ZBP_expl(i) = interp1([r_grid(idx),r_grid(idx+1)],...
                          [V_for(idx,end),V_for(idx+1,end)],r0);

    ZBP_CN(i) = interp1([r_grid(idx),r_grid(idx+1)],...
                        [V_CN(idx,end),V_CN(idx+1,end)],r0);
    
end

%% Plot 2 analyse of convergence with respect to the time grid

figure
hold off
plot(N_ts, ZBP*ones( length(N_ts),1), 'k', 'LineWidth',2)
hold on
plot(N_ts, ZBP*ones( length(N_ts),1)*1.01, 'k--','LineWidth',2)
plot(N_ts, ZBP*ones( length(N_ts),1)*0.99, 'k--','LineWidth',2)
plot(N_ts, ZBP_CN, '.-','LineWidth',1,'MarkerSize',10)
plot(N_ts, ZBP_expl, '.-','LineWidth',2,'MarkerSize',10)
plot(N_ts, ZBP_impl, '.--','LineWidth',2,'MarkerSize',10)
legend('analytical', 'band +1%', 'band -1%','Crank-Nicholson', 'explicit','implicit' )
xlabel('N_t')
ylabel('ZBP(0,T1,T2,K)')


%% Supporting functions

function  [P,A,B] = ZC(tau,r,k,theta,sigma)

%OUTPUT: 
% P:  Zero Coupon Bond price P(t,T,r)
%
%INPUT: 
% tau: time to maturity
% r,k,theta,sigma: Vasicek parameter dr=k(theta-r)dt+sigma dW

B = (1 - exp(-k*tau))/k;
A = exp((theta - sigma^2/(2*k^2))*(B - tau) - sigma^2/(4*k)*B^2);

P = A*exp(-B*r);
end



