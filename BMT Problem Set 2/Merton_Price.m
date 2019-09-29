function P = Merton_Price(S,K,r,T,sigma,lambda_Q,gamma,J)

%=============================================================%
%====== Merton Jump-Diffusion Option pricing formula =========%
%=============================================================%

% OUTPUT:
% P: The price of a european call option

% INPUT:
% S:        Initial price of the option 
% K:        Payoff parameter
% T:        Maturity
% r:        Risk-free rate
% sigma :   Volatility 
% lambda_Q: Jump intensity parameter
% gamma:    Jump size
% J:        Maximum number of jumps considered in the Merton formula


% Modified initial asset price to take account of the jumps
S_new = S .* (1 + gamma).^(0:J) * exp(-(lambda_Q * gamma) * T);

% Approximated price where the following functions are used:
% blsprice: Black-Scholes price
% poisspdf: pdf of a Poisson distribution

P = sum(blsprice(S_new,K,r,T,sigma,0).* poisspdf(0:J,lambda_Q * T));



  