%=========================================================================%
%========================== Advanced Derivatives =========================% 
%============================== Problem Set 4 ============================%
%======== BRODARD Lionel, MARCHAL Antoine, TISSOT-DAGUETTE Valentin ======%
%=========================================================================%

close all; clear; clc; format short; warning('off')

rng(1331)

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

% Closing prices and dividend rates (amzn and spx respectively)
S_0 = [1971,2921]; delta = [0.019,0.018];

id = ["AMZN","SPX"];

%% I. Derive the marginal implied distributions

for i = 1:2
    
    % Rows where strikes and implied vols are available
    I = find(~ isnan(Table{: , id(i) + "_K"}));
    
    % Array of call prices for the different strikes (using BS_Price.m)
    C = exp(r*T) * BS_Price(S_0(i),Table{I , id(i) + "_K" },...
                            r , T, Table{I , id(i) + "_IV"},delta(i));
     
    % Strikes increments K_i - K_{i-1} for the Finite Difference Scheme
    D_K = diff(Table{I , id(i) + "_K" });
    
    % Use the first-order forward difference to derive the implied cdf
    Table{I(1:end-1),id(i) + "_cdf"} = 1 + diff(C)./ D_K;
                                                                                                     
    % Extrapolate the cdf by filling the end of the column with 1's
    Table{I(end):end,id(i)+"_cdf"} = 1; 
    
    % Plot of the implied distribution
    subplot(2,1,i) 
    plot(Table{:,id(i)+"_K"},Table{:,id(i)+"_cdf"},'Linewidth',2)
    xlabel(sprintf("S_T^{%s}",id(i)))
    ylabel("\Phi^{imp}")
    title(sprintf('Implied distribution for %s',id(i)))
    xlim([min(Table{:,id(i)+"_K"}),max(Table{:,id(i)+"_K"})])
    
end

% Display the head of the table
Table(1:20,:)


%% Monte Carlo and Gaussian Copula method

% Number of simulations 
N_sim = 1e4;

% Inverse of the implied marginal distribution: 
% (first index in the Table whose cdf value is above a given level x)
Impl_cdf_inv = @(x,id) Table{find(Table{:,id + "_cdf"} >= x,1), id + "_K"};

% Mean and covariance matrix for the multivariate Gaussian distribution
Mu = [0,0]; Sigma = [1,0.5;0.5,1];

% N_sim x 2 array for the inverse transform sampling
x = normcdf(mvnrnd(Mu,Sigma,N_sim));

% Array of terminal values for the stock prices
S_T = zeros(N_sim,2);

for i = 1:2
    for j = 1:N_sim
        
        % Use the inverse function defined above  
        S_T(j,i) = Impl_cdf_inv(x(j,i),id(i));
    end    
end

% Simulated payout values
H = max(S_T(:,2)/S_0(2) - S_T(:,1)/S_0(1),0);

fprintf('\nPrice of the exchange option %2.5f\n',exp(-r*T) * mean(H))
