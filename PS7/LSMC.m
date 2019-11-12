function [Price,Std] = LSMC(S_0,K,r,q,sig,t,t_ex,dt,N_MC,N_sim,show,int)

%=========================================================================%
% Least-Square Monte Carlo Method                                         %
%=========================================================================%
% INPUT:                                                                  %
% S_0   : Initial Stock Price                                             %
% K     : Strike                                                          %
% r     : Interest Rate                                                   %
% q     : Dividend Yield                                                  %
% sig   : Array of local volatilities                                     %
% t     : Simulation dates                                                %
% t_ex  : Exercise dates                                                  %
% dt    : Time step for the simulations                                   %
% N_MC  : Number of Monte Carlo simulations                               %
% N_sim : Number of simulations to compute the standard deviation         %
%=========================================================================%

% Number of simulation dates
N_t = length(t); 

% Number of exercise dates
N_ex = length(t_ex);

% Array of option prices
Prices = zeros(N_sim,1);

for s = 1:N_sim 
    
% Trajectories (Euler-Maruyama Scheme)
S = S_0 * cumprod(1 + (r-q)*dt + sqrt(dt)*sig.*randn(N_MC,N_t),2);

% Moving average of the stock prices
A = movmean(S,[N_t 0],2);

if show >= 1
    
    % Display path trajectories
    figure; subplot(2,1,1);
    plot([0,t],[repmat(S_0,30,1),S(1:30,:)],'Linewidth',1.2) 
    title('S: Stock price')
   
    subplot(2,1,2); 
    plot([0,t],[repmat(S_0,30,1),A(1:30,:)],'Linewidth',1.2)
    title('A: Running mean of the stock price');
    suptitle("Trajectories of S and A, " + int)
   
end
    
if show == 2 
    
    % histogram for the distribution of A_T and S_T
    figure 
    histogram(A(:,end),floor(N_MC^(1/3)),'Normalization','pdf'); hold on
    histogram(S(:,end),floor(N_MC^(1/3)),'Normalization','pdf')
    xlabel('Terminal Value')
    title('Distribution of S_T and A_T')
    legend('A_T','S_T')      
end 

show = 0; 

% Indices of t where one can exercise the option
id_ex = find(ismember(t,t_ex));

% Keep only the (average) stock prices at the exercise dates
S = S(:,id_ex); A = A(:,id_ex);

% Matrix of intrinsic values
Int_Values = max(A - K,0);

% Cash Flow Matrix 
CF = zeros(N_MC,N_ex);

% Store the payouts at maturity
CF(:,end) = Int_Values(:,end);

for k = N_ex-1:-1:1
    
    % Row indices where the intrisic value is positive at time k
    id = find(Int_Values(:,k) > 0);
    
    % Design Matrix
    X = [ones(length(id),1),S(id,k),S(id,k).^2,S(id,k).^3];
                
    if k > 1
        X = [X,A(id,k),A(id,k).^2,A(id,k).^3];
    end
    
    % Response Variable
    Y = sum(exp(-r *(t_ex(k+1:end) - t_ex(k))) .* CF(id,k+1:end),2);

    % Continuation values obtained using OLS regression
    Cont_Values =  X * ((X'*X)\(X'*Y));
    
    % Get the indices where it is better to exercise
    id_ex = find(Int_Values(id,k) > Cont_Values);
    
    % Add the corresponding intrisic values in the cash flow matrix
    CF(id(id_ex),k) = Int_Values(id(id_ex),k);
    
    % Set all future cash flows to 0 when the option is exercised at time k
    CF(id(id_ex),k+1:end) = 0;
    
end

% Average over all paths to get the price of the option  
Prices(s) = mean(sum(CF .* exp(- r * t_ex),2));

end

% Compute the average price and its standard deviation
Price = mean(Prices);   Std = std(Prices) ;