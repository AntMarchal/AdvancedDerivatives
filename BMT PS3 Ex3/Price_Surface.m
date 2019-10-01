function P = Price_Surface(S_0,K,r,T,sigma,lambda_Q,gamma,J)

%=============================================================%
%====== Merton Jump-Diffusion Option pricing formula =========%
%=============================================================%

% OUTPUT:
% P: The price(s) of a european call option

% INPUT:
% S_0:      Initial price of the option 
% K:        Strike (can be a vector)
% T:        Maturity (can be a vector)
% r:        Risk-free rate
% sigma :   Volatility 
% lambda_Q: Jump intensity parameter
% gamma:    Jump size
% J:        Maximum number of jumps 

% Array of prices
P = zeros(length(K),length(T));

for t = 1:length(T)
    
    % Modified initial asset price to take account of the jumps
    S_new = S_0 .* (1 - gamma).^(0:J) * exp(lambda_Q * gamma * T(t));
    
    for k = 1:length(K)
        
        % Approximated price where the following functions are used:
        % blsprice: Black-Scholes price
        % poisspdf: pdf of a Poisson distribution
        P(k,t) = sum(blsprice(S_new,K(k),r,T(t),sigma,0)...
              .* poisspdf(0:J , lambda_Q * T(t)));

    end
end

  