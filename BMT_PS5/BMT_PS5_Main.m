%=========================================================================%
%========================== Advanced Derivatives =========================% 
%============================== Problem Set 5 ============================%
%======== BRODARD Lionel, MARCHAL Antoine, TISSOT-DAGUETTE Valentin ======%
%=========================================================================%

close all
clear all
clc

%import the data
data = xlsread('SX5E_Impliedvols.xlsx');

S_0=2772.7; %initial price
K=data(2:end,1)*S_0; %strikes
T=data(1,2:end); %time to maturities

dK=K(2)-K(1); % Step size for the strikes 

vol = data(2:end,2:end); 

% Define the function f(K,t) = K * sigma(K,t)
f = diag(K)*data(2:end,2:end);

% Interest and dividend 
r=0; q=0;

%% Compute the known prices using BS and the known volatilities

C_0=max(S_0-K,0);

C_obs=zeros(length(K),length(T));

for i=1:length(K)
    for j=1:length(T)
        if vol(i,j)~=0
            C_obs(i,j)=CallBS(S_0,K(i),r,T(j),vol(i,j),q);
        else
            C_obs(i,j)=0;
        end
    end
end  

%% Estimates the volatilities of the A matices, then compute all the call prices

%%%%%%%%%%%%% Setup %%%%%%%%%%%%%

f_est = zeros(size(vol));

C_models = zeros(size(C_obs));

fvals=zeros(length(T),1);

dT=[T(1),diff(T)]; Cj_1 = C_0; len=length(C_0);

for i=1:length(T)
    
    %first guess corresponding to the observed volatilities
    para0 = nonzeros(f(:,i));
    
    %%%%%%%%%% Optimization %%%%%%%%%%%%%%%
    % We look for the optimal sigma_tilde = sigma(K,T) * K 
   
    % Increase the max number of iterations and turn off the displays
    opts = optimset('MaxIter',1e3,'Display','off');
    
    % Upper and lower bound
    ub = repmat(S_0,length(para0),1); lb = zeros(length(para0),1);
         
    % It turns out that imposing constraints for the volatility doesn't
    % change the result -> simply use fminsearch
    [f_est_tmp,fval,~,~] = fmincon(@(para)opti_fct(Cj_1,C_obs(:,i),dT(i),...
                           dK,para) ,para0,[],[],[],[],lb,ub,[],opts);
                       
    fvals(i) = fval; % store the fvals
    
    %%%%%%%%%% Predictions %%%%%%%%%%%%%%%%%
    
    positions= find(C_obs(:,i)); % positions of observed call price
    
    sigma = sigma_construct(f_est_tmp, positions,len);
    
    f_est(:,i) = sigma; % stock the estimated sigma
    
    % Use the function matA which constructs the matrix A
    C_model= max(matA(sigma,dT(i),dK,len)\Cj_1,0); 
    
    C_models(:,i)=C_model; % stock the modelled prices
    
    % The prices obtained becomes the starting prices C_{j-1}
    % for the next iteration
    Cj_1 = C_model; 
end

%%  Compute the price for T=1 and T=1.5

T_prime = [1,1.5];

C_models_prime=zeros(len,length(T_prime));

for i=1:length(T_prime)
    
    %T_prime between T_j and T_{j+1} for some j
    
    j = find(T>T_prime(i),1)-1; 
    
    dT = T_prime(i)-T(j);
    
    sigma = f_est(:,j);
    
    % Use the function matA again
    C_model = max(matA(sigma,dT,dK,len)\C_models(:,j),0);
    
    C_models_prime(:,i) = C_model;
    
end 
        
%% Insert the two new column into our previous matrix of call price%%%%%

%The new columns for T in {1,1.5} will be inserted in C_models_total
C_models_total=C_models; T_total = T;

for i=1:length(T_prime)
    
    j=find(T_total>T_prime(i),1)-1;
    
    T_total = [T_total(1:j),T_prime(i),T_total(j+1:end)];
    
    C_models_total=[C_models_total(:,1:j),C_models_prime(:,i),C_models_total(:,j+1:end)];
    
end
    
%% Compute the implied BS volatilities from the computed call prices

% Black-Scholes implied volatilities of call prices for the expirations in the spreadsheet
%and for the expirations T = 1 and T = 1.5

sig0 = 0.2;

[impv_total,fval3] = implvBS(S_0,K,r,T_total,q,C_models_total,sig0);

figure

surf(T_total,K,impv_total)

%%%%%%%%% figure %%%%%%%%%

%surf(T,K(2:end-1),impv(2:end-1,:))
title('Volatility surface')
xlabel('Maturities')
ylabel('Strikes')
zlabel('Volatility')
hold on

% Plot the observed values
[X,Y] = find(vol>0); 
Z = vol(vol>0); 
plot3(T(Y),K(X),Z,'.r','markersize',10)
view(120,15)

%% Functions

function C = CallBS(S,K,r,T,sigma,q)
% Computes the value of the Call using BS formula

d1 = (log(S/K)+(r-q+sigma^2/2)*T)/(sigma*sqrt(T));
d2 = d1-sigma*sqrt(T);

C =  S*exp(-q*T)*normcdf(d1) - K*exp(-r*T)*normcdf(d2);

end

function tridiag = matA(sigma,dT,dK,n)

Z = 0.5*sigma.^2*dT/(dK)^2;

nones = ones(n, 1);

tridiag = diag(Z)*(2*diag(nones) - diag(nones(1:end-1), -1)...
        - diag(nones(1:end-1), 1)) + diag(nones);
    
%replace first and last row by (1,0,...,0) and (0,...,0,1)
tridiag(1,:) = [1,zeros(1,n-1)];
tridiag(n,:) = [zeros(1,n-1),1];

end


function f=opti_fct(Cj_1,Cobsj,dT,dK,sigma_est)

len=length(Cj_1);

positions= find(Cobsj); %positions of observed call price

sigma=sigma_construct(sigma_est, positions,len);

C_model=pinv(matA(sigma,dT,dK,len))*Cj_1;

%compare the prediction for which it exists an observation
f=sum((C_model(positions)-Cobsj(positions,1)).^2);

end

function sigma=sigma_construct(sigma_est, positions,len)

% Interpolation of sigmas between the volatilities. 
% By using the 'nearest' option, the resulting function is piecewise
% constant. For a point where the volatility is now known yet, 
% it uses the closest observe value as image, which is exactly in line with
% the graphical example in the problem statement.

sigma = interp1(positions,sigma_est,(1:len)','nearest');

%sigmas before the first estimated sigma
sigma(1:positions(1))=sigma(positions(1)); 

%sigmas after the last estimated sigma
sigma(positions(end):end)=sigma(positions(end)); 

end

function [impv,fval]=implvBS(S,K,r,T,q,C,sig0)
impv=zeros(size(C));
for i = 1:length(K)
    for j = 1:length(T)
        
        fct = @(sigma) (CallBS(S,K(i),r,T(j),sigma,q) - C(i,j));
        
        % The implied vol is obtained as the zero of the function fct 
        [impv(i,j),fval(i,j)] = fzero(fct,sig0);
    end
end
end