%=========================================================================%
%========================== Advanced Derivatives =========================% 
%============================== Problem Set 2 ============================%
%======== BRODARD Lionel, MARCHAL Antoine, TISSOT-DAGUETTE Valentin ======%
%=========================================================================%

close all
clear

%=========================================================================%
%================================ Exercise 1 =============================%
%=========================================================================%

fprintf('\n=================== Exercise 1 =====================\n')

%% 0. Setup
S = 100; D = 50; r = 0.05; sig = 0.4; T_2 = 10; q = 0;

% Computation of the current value of the firm's equity
E = blsprice(S,D,r,T_2,sig);

fprintf('\nCurrent value of the equity: %2.2f\n',E)

% Array of strikes and maturities
K_1 = [0.6,0.8,1,1.2].* E; T_1 = [2,5,7,9];
                                        
%% I: Find the compound price and extract the implied volatility 
IV = zeros(length(K_1),length(T_1)); 

% Store the errors when computing the implied volatilities
fval = IV;
          
for i=1:length(K_1)
    
    for j=1:length(T_1)
        
        % Find S^* st C_BS(S^*,D,T_2 - T_1) - K_1 = 0
        S_star = fzero(@(S) blsprice(S,D,r,T_2 - T_1(j),sig,q)- K_1(i),S);
        
        % Use the function Compound_Call
        Compound_Price = Compound_Call(S,S_star,T_1(j),T_2,K_1(i),D,r,sig);
        
        % Deduce the implied volatility
        [IV(i,j),fval(i,j)] = fzero(@(sig) blsprice(E,K_1(i),r,T_1(j),sig)...
                            - Compound_Price,0.5);
             
    end
end

fprintf('\nMax error when finding the implied vol: %e\n',max(fval,[],'all'))

fprintf('\nImplied Volatilities:\n'); IV

%% II. Plot the implied volatility surface

figure
surf(T_1,K_1,IV, 'FaceAlpha',0.5); colorbar
zlabel('Implied volatility'); ylabel('Strike'); xlabel('Time to maturity')
set(gca, 'YDir','reverse'); view(-60,15)
zlim([0.98 * min(IV,[],'all'),1.02 * max(IV,[],'all')])
suptitle('Exercise 1: Implied Volatility Surface')

%=========================================================================%
%================================ Exercise 2 =============================%
%=========================================================================%

fprintf('\n=================== Exercise 2 =====================\n')

%% 0. Setup

gamma = -0.08; lambda_Q = 0.2; sigma = 0.2; r = sigma^2; 

T = [0.02,0.08,0.25,0.5]; S = 100; K = S * (0.8:0.1:1.1);

%% I. Determination of the maximum number of jumps

% Our aim is to find a suitable maximum number of jumps J
% to truncate the infinite sum appearing in the Merton formula

% Vector of candidates
J_candidates = 1:20;

% Vector that will contain the approximated call prices
Call_Merton = zeros(length(J_candidates),1);

for j = J_candidates
    Call_Merton(j) = Merton_Price(S,S,T(1),r,sigma,lambda_Q,gamma,j);
end

% Return the first integer J s.t. the prices truncated at the level
% J and J+1 are identical up to machine precision (eps)
J = find((Call_Merton(2:end) - Call_Merton(1:end-1)) < eps == 1,1);

fprintf('\nTruncation level required: J = %d\n',J)

%% II. Plots displaying the implied volatility as a function of strike

maturity = {"1 week","1 month","3 months","6 months"};

figure

Call_Merton = zeros(length(K),length(T));

for t = 1:length(T)
    
    % Store the implied volatilities in a vector
    imp_vol = zeros(length(K),1);
    
    for k = 1:length(K)
        
        % Use the function Merton_price:
        Call_Merton(k,t) = Merton_Price(S,K(k),r,T(t),sigma,lambda_Q,gamma,J);
        
        % Compute the implied volatility thanks to the function blsimpv
        imp_vol(k) = blsimpv(S,K(k),r,T(t),Call_Merton(k,t));
    
    end
    
    % Plot
    subplot(2,2,t)
    plot(K,imp_vol,'o-','linewidth',1.5); 
    xlabel('K'); ylabel('Implied Volatility'); ylim([0.19 0.33]);
    title(sprintf('T = %s',maturity{t}))
end

fprintf('\nCall Option Prices:\n'); Call_Merton

suptitle('Exercise 2: Implied Volatility for different maturities')