clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Advanced Derivative: PS 9 %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Lionel BRODARD, Antoine MARCHAL & Valentin TISSOT-DAGUETTE %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameters
K       = 0.805;
k       = 0.15;
sigma_r = 0.01;
theta   = 0.05;
r0      = 0.042;
T1      = 0.25;
T2      = 5+T1;

%% Analytical Solution

B = (1 - exp(-k*(T2-T1)))/k;
sig_tilde = sigma_r * sqrt((1-exp(-2*k*T1))/(2*k)) * B;
h = 1/sig_tilde * log(ZC(T2, r0, k, theta, sigma_r)...
    /(ZC(T1, r0, k, theta, sigma_r)*K)) + sig_tilde/2;

ZBP = K * ZC(T1,r0,k,theta,sigma_r) * normcdf(-h+sig_tilde)...
      - ZC(T2,r0,k,theta,sigma_r) * normcdf(-h);

%% Finite difference

% Space domain truncation
r_max = 0.8;
r_min = 0;

%Initial boundary condition
initial_cond = @(r) max(0,K-ZC(T2-T1,r,k,theta,sigma_r));
%Left boundary condition
bc_left = 0;
%Right boundary condition
bc_right = K-ZC(T2-T1,r_max,k,theta,sigma_r);


%% Rate of convergence with respect to the space grid

N_rs = [200:10:600]; % N_r \in N_rs will be the Number of equal length interval in the space grid
% Time discretization
N_t = 200; %Number of equal length interval in the time grid

ZBP_expl = zeros (length(N_rs),1);
ZBP_impl = zeros (length(N_rs),1);
ZBP_CN = zeros (length(N_rs),1);

for i=1:length(N_rs)
    N_r=N_rs(i);
    % Explicit method: backward finite time difference (\nu=1)
    [V_back,r_grid,t_grid]=FiniteDiff(sigma_r,k,theta,bc_left,bc_right,initial_cond,r_max,r_min,N_r,T1,N_t,1);
    % Implicit method: forward finite time difference (\nu=0)
    [V_for,r_grid,t_grid]=FiniteDiff(sigma_r,k,theta,bc_left,bc_right,initial_cond,r_max,r_min,N_r,T1,N_t,0);
    % Crank-Nicholson method
    [V_CN,r_grid,t_grid]=FiniteDiff(sigma_r,k,theta,bc_left,bc_right,initial_cond,r_max,r_min,N_r,T1,N_t,0.5);
    
    % interpolation of the value corresponding to r0
    idx = find(r_grid<r0, 1, 'last');
    V0r0_back = interp1([r_grid(idx),r_grid(idx+1)], [V_back(idx,end),V_back(idx+1,end)],r0);
    V0r0_for = interp1([r_grid(idx),r_grid(idx+1)], [V_back(idx,end),V_back(idx+1,end)],r0);
    V0r0_CN = interp1([r_grid(idx),r_grid(idx+1)], [V_CN(idx,end),V_CN(idx+1,end)],r0);
    
    ZBP_expl(i) = V0r0_back;
    ZBP_impl(i) = V0r0_for;
    ZBP_CN(i) = V0r0_CN;
    
    
end



%% Plot 1 analyse of convergence with respect to the space grid
figure(1);

plot(N_rs, ZBP*ones( length(N_rs),1), 'k','LineWidth',2)
hold on
plot(N_rs, ZBP*ones( length(N_rs),1)*1.01, 'k--','LineWidth',2)
plot(N_rs, ZBP*ones( length(N_rs),1)*0.99, 'k--','LineWidth',2)
plot(N_rs, ZBP_CN, 'b.-','LineWidth',2,'MarkerSize',10)
plot(N_rs, ZBP_expl, 'r.-','LineWidth',2,'MarkerSize',10)
plot(N_rs, ZBP_impl, 'g.--','LineWidth',2,'MarkerSize',10)
legend('analytical', 'band +1%', 'band -1%','Crank-Nicholson', 'explicit','implicit' )
xlabel('N_r')
ylabel('ZBP(0,T1,T2,K)')


%% Rate of convergence in time
N_r = 400;
N_ts = [5:5:150];
ZBP_expl = zeros (length(N_ts),1);
ZBP_impl = zeros (length(N_ts),1);
ZBP_CN = zeros (length(N_ts),1);

for i=1:length(N_ts)
    N_t=N_ts(i);
    % Explicit method: backward finite time difference (\nu=1)
    [V_back,r_grid,t_grid]=FiniteDiff(sigma_r,k,theta,bc_left,bc_right,initial_cond,r_max,r_min,N_r,T1,N_t,1);
    % Implicit method: forward finite time difference (\nu=0)
    [V_for,r_grid,t_grid]=FiniteDiff(sigma_r,k,theta,bc_left,bc_right,initial_cond,r_max,r_min,N_r,T1,N_t,0);
    % Crank-Nicholson method
    [V_CN,r_grid,t_grid]=FiniteDiff(sigma_r,k,theta,bc_left,bc_right,initial_cond,r_max,r_min,N_r,T1,N_t,0.5);
    
    % interpolation of the value corresponding to r0
    idx = find(r_grid<r0, 1, 'last');
    V0r0_back = interp1([r_grid(idx),r_grid(idx+1)], [V_back(idx,end),V_back(idx+1,end)],r0);
    V0r0_for = interp1([r_grid(idx),r_grid(idx+1)], [V_back(idx,end),V_back(idx+1,end)],r0);
    V0r0_CN = interp1([r_grid(idx),r_grid(idx+1)], [V_CN(idx,end),V_CN(idx+1,end)],r0);
    
    ZBP_expl(i) = V0r0_back;
    ZBP_impl(i) = V0r0_for;
    ZBP_CN(i) = V0r0_CN;
    
end

%% Plot 2 analyse of convergence with respect to the time grid
figure(2);
hold off
plot(N_ts, ZBP*ones( length(N_ts),1), 'k', 'LineWidth',2)
hold on
plot(N_ts, ZBP*ones( length(N_ts),1)*1.01, 'k--','LineWidth',2)
plot(N_ts, ZBP*ones( length(N_ts),1)*0.99, 'k--','LineWidth',2)
plot(N_ts, ZBP_CN, 'b.-','LineWidth',1,'MarkerSize',10)
plot(N_ts, ZBP_expl, 'r.-','LineWidth',2,'MarkerSize',10)
plot(N_ts, ZBP_impl, 'g.--','LineWidth',2,'MarkerSize',10)
legend('analytical', 'band +1%', 'band -1%','Crank-Nicholson', 'explicit','implicit' )
xlabel('N_t')
ylabel('ZBP(0,T1,T2,K)')



%% Supporting functions

function  P = ZC(tau,r,k,theta,sigma)

%OUTPUT: 
% P:  Zero Coupon Bond price P(t,T,r)
%INPUT: 
% tau: time to maturity
% r,k,theta,sigma: Vasicek parameter dr=k(theta-r)dt+sigma dW

B = (1 - exp(-k*tau))/k;
A = exp((theta - sigma^2/(2*k^2))*(B - tau) - sigma^2/(4*k)*B^2);

P = A*exp(-B*r);
end



