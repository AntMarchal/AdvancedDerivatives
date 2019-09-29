function CompoundPrice = Compound_Call(S,S_star,tau1,tau2,K_1,K_2,r,sig)

%=====================================================%
%====== Compound Call Option pricing formula =========%
%=====================================================%

% OUTPUT:
% P: Price of a compound call 

% Thresholds
a = @(i)(log(S/S_star)+ (r -(-1)^i * sig^2/2)*tau1)/(sig*sqrt(tau1)); 
b = @(i)(log(S/K_2)   + (r -(-1)^i * sig^2/2)*tau2)/(sig*sqrt(tau2)); 

% Covariance matrix
SIGMA= [1,sqrt(tau1/tau2);sqrt(tau1/tau2),1];

% Pricing formula
CompoundPrice = S *  mvncdf([a(1),b(1)],[0,0],SIGMA)...
              - K_2 * exp(-r*tau2) * mvncdf([a(2) b(2)],[0,0],SIGMA)...
              - exp(-r*tau1) * K_1 * normcdf(a(2));

end