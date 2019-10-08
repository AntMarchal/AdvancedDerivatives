%=========================================================================%
%========================== Advanced Derivatives =========================% 
%============================== Problem Set 4 ============================%
%======== BRODARD Lionel, MARCHAL Antoine, TISSOT-DAGUETTE Valentin ======%
%=========================================================================%

close all; clear; clc; format short; warning('off')

%% 0. Setup

% Load the data
Table  = readtable("Impvols_SPX_AMZN.xlsx",'Range','A2:F280');

% Remove empty columns
Table  = removevars(Table ,{'Var3','Var4'});

% Remove the first row (empty)
Table = Table(2:end,:);

% Rename the columns (K: strikes, IV: ImpliedVols)
Table.Properties.VariableNames = {'SPX_K' 'SPX_IV' 'AMZN_K' 'AMZN_IV'};

% Parameters
T = 0.296; r = 0.024; 

% Closing prices and dividend rate (amzn and spx respectively)
S_0 = [1971,2921]; delta = [0.019,0.018];

id = ["AMZN","SPX"];

%% I. Derive the marginal implied distributions
%figure 

for i = 1:2
    
    % Rows where strikes and implied vols are available
    J = find(Table{: , id(i) + "_K"} > 0);
    
    % Temporary array of call prices for the different strikes 
    C_tmp = BS_price(S_0(i),Table{J , id(i) + "_K" },...
                     r , T, Table{J , id(i) + "_IV"},delta(i))
     
    % Strikes increments K_i - K_{i-1} for the Finite Difference Scheme
    D_K = diff(Table{J , id(i) + "_K" });
    
    % Use the first-order forward difference to derive the implied cdf
    Table{J(1:end-1),id(i)+"_cdf"} = min(1 + (C_tmp(2:end)-C_tmp(1:end-1))...
                                  ./ D_K, 1 );
                             
    % Extrapolate the cdf by filling the end of the column with "1"'s
    Table{J(end):end,id(i)+"_cdf"} = 1; 
    
    subplot(2,1,i) 
    plot(Table{:,id(i)+"_K"},Table{:,id(i)+"_cdf"},'Linewidth',2)
    xlabel(sprintf("S_T^{%s}",id(i)))
    ylabel("\Phi^{impv}")
    title(sprintf('Implied distribution for %s',id(i)))
    xlim([min(Table{:,id(i)+"_K"}),max(Table{:,id(i)+"_K"})])
    
end


Table

%% Monte Carlo and Gaussian Copula method

% Number of simulations 
N_sim = 1e4;

% Inverse of the implied marginal distribution
Impl_cdf_inv = @(x,id) Table{find(Table{:,id + "_cdf"} >= x,1), id + "_K"};

% Mean and covariance matrix for the multivariate Gaussian distribution
Mu = [0,0]; Sigma = [1,0.5;0.5,1];

% N_sim x 2 array for the inverse transform sampling
x = normcdf(mvnrnd(Mu,Sigma,N_sim));

% Array of terminal values for the stock prices
S_T = zeros(N_sim,2);
figure
for i = 1:2
    
    for j = 1:N_sim
 
        S_T(j,i) = Impl_cdf_inv(x(j,i),id(i));

    end
    subplot(2,1,i)
    histogram(S_T(:,i),40)    
    title(sprintf('Empirical distribution of S_T for %s',id(i)))
end

% Simulated payout values
H = max(S_T(:,2)/S_0(2) - S_T(:,1)/S_0(1),0);

% Price of the exchange option using the MC estimator
P = exp(-r*T) * mean(H)

%% Extra: comparison with the formula of assignment 1

% It cannot be accurate since the BS modeal assumes a constant volatility
% Hence the "best" BS vol is given as follows (/!\ nan):
sigma_BS = (sum(Table.AMZN_IV) + sum(Table.SPX_IV(J)))...
      ./(length(Table.AMZN_IV) + length(Table.SPX_IV(J)))

% Thresholds in the outperformance option at t = 0
d = @(i) ((delta(1)-delta(2))+ (-1)^i * sigma_BS^2/2)* sqrt(T)/sigma_BS;

P_BS = exp(-delta(2)*T) * normcdf(d(2))...
     - exp(-delta(1)*T) * normcdf(d(1))
 
