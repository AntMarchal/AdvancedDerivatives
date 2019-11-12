%=========================================================================%
%========================== Advanced Derivatives =========================% 
%============================== Problem Set 7 ============================%
%======== BRODARD Lionel, MARCHAL Antoine, TISSOT-DAGUETTE Valentin ======%
%=========================================================================%

close all; clear; clc; format short; warning('off')

% Freeze the random seed
rng(1)

%% 0. Setup

S_0 = 100; K = 98; dt = 1/4; t = dt:dt:1; N_t = length(t);

r = 0; q = 0.02; sig = 0.23; N_MC = 10.^(2:5); N_sim = 1e2;

% Approximated Prices and Standard Deviations
P = zeros(length(N_MC),2); Std = zeros(length(N_MC),2);

%% I. Least-Square Monte Carlo with constant volatility

for j = 1:length(N_MC)
    
% Boolean to show a histogram for the distribution of S_T and A_T
show = (j == length(N_MC));

int = "quarterly intervals";

% Please refer to the function LSMC_Price for further details
[P(j,1),Std(j,1)] = LSMC(S_0,K,r,q,sig,t,t,dt,N_MC(j),N_sim,show,int);
end

%% II. Least-Square Monte Carlo with Local Volatilities

% Get the adequate local volatilities from the previous assignment 
data = readtable('Volatility_Surface.csv');

% Extract the local volatilities for the corresponding moneyness
Local_Vol = data{data{:,'K_over_S'} == K/S_0,2:end}';

N_Vol = length(Local_Vol);

% Time to maturities from the volatility surface
TTM = zeros(N_Vol ,1);

for j = 1:N_Vol
    
    var_name = data.Properties.VariableNames{j+1};
    
    TTM(j) = str2double(var_name(3:end))/1e3;
end

% Use monthly intervals as time grid
dt_M = 1/12;  t_M = round(dt_M:dt_M:1,12);

% Interpolate the local volatilities on the time grid
sig_M = interp1(TTM,Local_Vol,t_M);

for j = 1:length(N_MC)
    
show = 2 * (j == length(N_MC)); int = "monthly intervals";

[P(j,2),Std(j,2)] = LSMC(S_0,K,r,q,sig,t_M,t,dt_M,N_MC(j),N_sim,show,int);

end

%% III. Results

vol = {'constant volatility','local volatilities'};

for m = 1:2
    
    fprintf('\nLeast-Square Monte Carlo with %s\n',vol{m});
    fprintf('\nPrice with %d MC simulations: %2.4f\n',N_MC(end),P(end,m))
    fprintf('\nCorresponding standard deviation: %2.4f\n',Std(end,m))
    fprintf('\n___________________________________________________\n');
    
end

% Comparison with the corresponding European call option price
call = BS_Price(S_0,K,r,1,sig,q);

fprintf('\nCorresponding European Call Option Price: %2.4f\n',call)

figure 
loglog(N_MC,Std,'o-','Linewidth',1.5); grid on; hold on
loglog(N_MC,exp(2)./sqrt(N_MC),'k--','Linewidth',1.5)

legend('LSMC with constant volatility',...
       'LSMC with local volatilies'   ,...
       " O(N_{MC}^{-1/2})" )
   
xlabel('Number of Monte Carlo Simulations')
ylabel('Standard deviation')
title('Standard deviation of the option price versus N_{MC} ')