%=========================================================================%
%========================== Advanced Derivatives =========================% 
%============================== Problem Set 3 ============================%
%======== BRODARD Lionel, MARCHAL Antoine, TISSOT-DAGUETTE Valentin ======%
%=========================================================================%

close all
clear

%=========================================================================%
%================================ Exercise 3 =============================%
%=========================================================================%

fprintf('\n=================== Exercise 3 =====================\n')

%% 0. Setup

gamma = 0.1; lambda_Q = 1; sig = 0.2; r = sig^2; 

T = [0.02,0.08]; S_0 = 100;

%% I. Truncation Level

% As in the previous assignment, we look for a suitable maximum number
% of jumps J to truncate the infinite sum in the Merton formula

% Vector of candidates
J_candidates = 1:10;

% Vector that will contain ATM call prices
Call = zeros(length(J_candidates),1);

for j = J_candidates

    Call(j) = Price_Surface(S_0,S_0,T(1),r,sig,lambda_Q,gamma,j);
end

J = find((Call(2:end) - Call(1:end-1)) < eps == 1,1);

fprintf('\nTruncation level required: J = %d\n',J)

%% II. Computation of the implied probability distribution function (pdf)

% Choose lower and upper bound for the strike grid
K_min = 75; K_max = 125;

% Choose a sufficiently small step size
Delta_K = 1e-1;

% Inner grid of strikes
K_inner_grid = (K_min:Delta_K:K_max)';

fprintf('\nNumber of points on the inner grid: %d\n',length(K_inner_grid))

% Consider a sligthly larger grid to get the implied pdf
% also at the boundaries (as we will use a centered approximation)
K_grid = [K_min-Delta_K ; K_inner_grid ; K_max + Delta_K];

% Surface of undiscounted call option prices C(K,T)
C_K_T = exp(r*T).* Price_Surface(S_0,K_grid,r,T,sig,lambda_Q,gamma,J);

% Use a second-order central difference to derive the implied pdf
Implied_pdf = (C_K_T(3:end,:) - 2 * C_K_T(2:end-1,:)...
            +  C_K_T(1:end-2,:))  * Delta_K^(-2)      ; 

%% Exact probability distribution of S_T in the Merton model

% Call the function Merton_pdf (see the corresponding file attached)
Exact_pdf = Merton_pdf(K_inner_grid,S_0,r,T,sig,lambda_Q,gamma,J);

% Maximum approximation error using the implied pdf
Max_error = max(abs(Exact_pdf - Implied_pdf));

%% Visual comparison between the exact and implied pdf

maturity = {'1 week','1 month'};
 
figure

for t = 1:length(T)
    
    subplot(2,1,t)
    
    plot(K_inner_grid,Exact_pdf(:,t),'color',[0.4,0.7,0],'linewidth',2)
    hold on
    plot(K_inner_grid',Implied_pdf(:,t),'k:','linewidth',2)
    
    xlabel('S'); ylabel('\phi(S,T)');
    
    set(get(gca,'ylabel'),'rotation',0)
    ylim([-1e-4,1.1 * max(Implied_pdf,[],'all')])
    
    legend('Exact pdf','Implied pdf')
    title(sprintf('T = %s',maturity{t}))
    
    fprintf('\nMaximum error (T = %s): %e\n', maturity{t},Max_error(t))

end

suptitle('Exercise 3: Exact and Implied probability distribution for different maturities')
