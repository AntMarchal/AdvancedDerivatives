function [C,P]= BS_price(S_0,K,r,T,sigma,delta)

%=========================================================================%
%============== Black-Scholes Formula for vanilla options ================%
%=========================================================================%

% Thresholds
d = @(i) (log(S_0./K) + (r - delta - (-1)^i .* sigma.^2/2) .* T)...
          ./ (sqrt(T).* sigma);

% Call Option Price(s)
C = S_0.* exp(- delta * T) .* normcdf(d(1))...
  - K  .* exp(- r * T)     .* normcdf(d(2));

% Put Option Price(s)
P = K  .* exp(-r * T)      .* normcdf(-d(2))...
  - S_0.* exp(- delta * T) .* normcdf(-d(1));
