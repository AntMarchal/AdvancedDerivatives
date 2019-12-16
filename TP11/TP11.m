clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Advanced Derivative: PS 11 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Lionel BRODARD, Antoine MARCHAL & Valentin TISSOT-DAGUETTE %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
K   = 0.05;
T_s = 1;
T_x = 2;
dt  = 0.25;
T_N = 4;
sigma = [ones(1,7)*0.2, ones(1,4)*0.22, ones(1,5)*0.24];

%% Paths simulation and intrinsic value of swap
Nsim = 10000;
N_x = length(dt:dt:T_x);
N = length(dt:dt:T_N);

[F_paths,V] = PathSimulation(Nsim,N_x,N,dt,sigma,T_s,K);




%% Optimal exercise rule
%Keep only intrinsic value at the exercise dates
V = V(:,5:end);
%max(V(V>0))\simeq 0.0498
Hs = 0.0001:0.0001:0.05;

exercise_idx = zeros(size(V));

%Matrix of optimal H
OptimalH = zeros(4,1);
%Two periods scheme
%Initialisation, at the bigining the second period correspond to the last
%value, this will be updated at each iteration
V_2 = V(:,end);

for i=1:4
    %Swaption price for 2 periods
    SP = zeros(length(Hs),1);
    V_1 = V(:,end-i);
    DF = zeros(Nsim,1);
    DF(:,1) = 1./(1+dt*F_paths(:,9-i,9-i));


    for h=1:length(Hs) 
        H = Hs(h);
        CF_2period = TwoPeriodSwationCF(V_1,V_2,H,DF);
        %Swaption price for 2 periods
        SP(h) = mean(CF_2period);  
    end
    %Optimal H
    [argvalue, argmax] = max(SP);
    OptimalH(i) = Hs(argmax);
    %Assign continuation value for next iteration knowing the optimal H
    V_2 = TwoPeriodSwationCF(V_1,V_2,OptimalH(i),DF);
end


%% Out-of-sample price computation using optimal exercise rule

[F_paths,V] = PathSimulation(Nsim,N_x,N,dt,sigma,T_s,K);
%Initialisation
V_2 = V(:,end);


for i=1:4
    V_1 = V(:,end-i);
    DF = zeros(Nsim,1);
    DF(:,1) = 1./(1+dt*F_paths(:,9-i,9-i));
    CF = TwoPeriodSwationCF(V_1,V_2,OptimalH(i),DF);
    V_2 = CF;
end
SwaptionPrice = mean(CF);






%% Supporting Functions

function [F_paths,V]=PathSimulation(Nsim,N_x,N,dt,sigma,T_s,K)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
% INPUT:
%%%%%%%%%%%%%%%%%%
% N_sim : Number of simulations
% N_x   : Number of time steps before maturity T_x (T_0=0 not included)
% N     : Total number of time steps (T_0=0 not included)
% dt    : Time interval
% sigma : Volatility, model parameter (dim N vector)
% T_s   : First exercise date
% K     : Strike
%%%%%%%%%%%%%%%%%%
% OUTPUT:
%%%%%%%%%%%%%%%%%%
% F_paths : Simulated Forward rate paths dim(F_paths)= Nsim X (N_x+1) X N
% V       : 2d-Matrix of intrinsic swap values for each simulated path,
%           during exercise time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F_paths = zeros(Nsim,N_x+1,N);

%Define matrix of intrinsic values
V = zeros(Nsim, N_x+1);

% Define theta as a function of k
theta = @(k) (pi*(k-2))/28;

% Flat term structure in t=0, F_k(0) = O.05 \forall k
F_paths(:,1,:) = 0.05*ones(Nsim,1,N);



for t=2:N_x+1
    nu = 0;
    for k=t:N
        nu = nu + dt*cos(theta(k)-theta(t))*sigma(k)*F_paths(:,t-1,k)./...
            (1+dt*F_paths(:,t-1,k));
        mu = sigma(k)*nu;
        dW = sqrt(dt)*(cos(theta(k))*randn(Nsim,1) + sin(theta(k))*...
            randn(Nsim,1));
        F_paths(:,t,k) = F_paths(:,t-1,k).*exp((mu-0.5*sigma(k)^2)*dt+...
                         sigma(k)*dW);
    
    end
    % Compute intrinsic value of the swap between T_s and T_x
    if t>(T_s/dt)
        P = zeros(Nsim,N);
        P(:,:) = cumprod(1./(1+dt*F_paths(:,t,:)),3);
        P(:,1:t-1)=0;
        V(:,t) = dt*diag(P*(reshape(F_paths(:,t,:)-K,[Nsim,N]))');
             
    end
    
end

end




function V = TwoPeriodSwationCF(V_1,V_2,H,DF)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute vector of discounted cash flow
%%%%%%%%%%%%%%%%%%
% INPUT:
%%%%%%%%%%%%%%%%%%
% V_1 : Intrinsic value period 1
% V_2 : Intrinsic value period 2
% H   : Threshold
% DF  : Discount factor for the second period
%%%%%%%%%%%%%%%%%%
% OUTPUT:
%%%%%%%%%%%%%%%%%%
% V : Vector of discounted cash flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %Two periods scheme
 exercise_idx = zeros(length(V_1),2);
 %Cash flow matrix
 CF = zeros(length(V_1),2);
    
 %First period, Exercise if intrinsic value >H
 exercise_idx(:,end-1) = (V_1>H);
 %If we did not exercise at T_{X_1}, we exercise at T_X if ITM
 exercise_idx(:,end) = (V_2>0).*(~exercise_idx(:,end-1));
 exercise_idx = logical(exercise_idx);
 CF(exercise_idx(:,end-1),end-1) = V_1(exercise_idx(:,end-1),1);
 %Discount the second cash flow using F_x(T_{x-1}) as the libor rate
 %L(T_{x-1},T_x)
 CF(exercise_idx(:,end),end) = V_2(exercise_idx(:,end),1).*DF(exercise_idx(:,end),1);
 
 V = sum(CF,2);

end
