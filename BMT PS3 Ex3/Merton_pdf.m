function PDF = Merton_pdf(S_grid,S_0,r,T,sig,lambda_Q,gamma,J)

%=================================================================%
%=== Probability distribution of S_T in the Merton model =========%
%=================================================================%

% OUTPUT:
% PDF: Probability distribution function of S_T

% INPUT:
% S_grid:   Grid of stock prices to compute the pdf
% S_0:      Initial price of the option 
% K:        Strike (can be a vector)
% T:        Maturity (can be a vector)
% r:        Risk-free rate
% sigma:    Volatility 
% lambda_Q: Jump intensity parameter
% gamma:    Jump size
% J:        Maximum number of jumps

PDF = zeros(length(S_grid),length(T));

for t = 1:length(T)
    
    % Threshold
    d = - (log(S_0 * (1 - gamma).^(0:J)' ./ S_grid')       ...
      + (r + lambda_Q * gamma) * T(t) )/ (sig * sqrt(T(t)))...
      + 1/2 * sig * sqrt(T(t));
  
    % Deduce the pdf of S_T
    PDF(:,t) = sum(normpdf(d)./(S_grid' * sig * sqrt(T(t)))...
            .* poisspdf(0:J,lambda_Q * T(t))');

end