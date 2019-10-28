clear all
clc

%import data
data=xlsread('SX5E_Impliedvols.xlsx');


S0=2772.7; %initial price

K=data(2:end,1)*S0; %strikes

T=[data(1,2:end)]; %time to maturities

dK=K(2)-K(1); % Step size for the strikes 

vol=data(2:end,2:end);

vol_tilde=diag(K)*data(2:end,2:end);

% Interest and dividend 
r=0; q=0;

%% Compute the known prices using BS and the known volatilities

C0=max(S0-K,0);

C_obs=zeros(length(K),length(T));

for i=1:length(K)
    for j=1:length(T)
        if vol(i,j)~=0
            C_obs(i,j)=CallBS(S0,K(i),r,T(j),vol(i,j),q);
        else
            C_obs(i,j)=0;
        end
    end
end

%% Estimates the volatilities of the A matices, then compute all the call prices

% Setup
vol_est=zeros(size(vol));

C_models=zeros(size(C_obs));

fvals=zeros(length(T),1); 

dT=[T(1),diff(T)]; Cj_1 = C0; len=length(C0);
for i=1:length(T)
    para0 = nonzeros(vol_tilde(:,i)); %first guess corresponding to the observed volatilities
    UB = ones(1,length(para0))*S0; % Upper bound for sigma<S0
    LB = zeros(1,length(para0)); %Lower bound sigma>0
    %%%%%%%%%%Optimization%%%%%%%%%%%%%%%
    [para,fval,~,~]=fminsearchcon(@(para)opti_fct(Cj_1,C_obs(:,i),dT(i),dK,para) ,para0,LB,UB);
    fvals(i)=fval; %stock the fvals
    sigma_est=para;
    %%%%%%%%%%Predictions%%%%%%%%%%%%%%%%%
    positions= find(C_obs(:,i)); %positions of observed call price
    sigma=sigma_construct(sigma_est, positions,len);
    vol_est(:,i)=sigma; %stock the estimated sigma
    C_model=pinv(matA(sigma,dT(i),dK,len))*Cj_1; %C_j=A^{-1}C_{j-1}
    C_models(:,i)=C_model; %stock the modelled prices
    Cj_1=C_model; %the prices obtained becomes the starting prices C_{j-1} for the next iteration
end


       

%%  Compute the price for T=1 and T=1.5
T_prime=[1,1.5];
C_models_prime=zeros(len,length(T_prime));
for i=1:length(T_prime)
    j=find(T>T_prime(i),1)-1; %T_prime between T_j and T_{j+1}
    dT=T_prime(i)-T(j);
    sigma=vol_est(:,j);
    C_model=pinv(matA(sigma,dT,dK,len))*C_models(:,j);
    C_models_prime(:,i)=C_model;
end
 
        
%% Insert the two new column into our previous matrix of call price%%%%%
C_models_total=C_models; %In the end the to new columns will be inserted in C_models_total
T_total=T;
for i=1:length(T_prime)
    j=find(T_total>T_prime(i),1)-1;
    T_total=[T_total(1:j),T_prime(i),T_total(j+1:end)];
    C_models_total=[C_models_total(:,1:j),C_models_prime(:,i),C_models_total(:,j+1:end)];
end
    
%% Compute the implied BS volatilities from the computed call prices
% Black-Scholes implied volatilities of call prices for the expirations in the spreadsheet
sig0=0.2;
[impv,fval2]=implvBS(S0,K,r,T,q,C_models,sig0);
% Black-Scholes implied volatilities of call prices for the expirations in the spreadsheet
%and for the expirations T = 1 and T = 1.5
[impv_total,fval3]=implvBS(S0,K,r,T_total,q,C_models_total,sig0);
%%%%%%%%%figure%%%%%%%%%
figure(1)
surf(T,K,impv)
title('Volatility surface')
xlabel('Maturities')
ylabel('Strikes')
zlabel('Volatility')
hold on
%Plot les observed values
[X,Y]=find(vol>0); 
Z=vol(vol>0); 
plot3(T(Y),K(X),Z,'.r','markersize',10)
hold off
%% Functions

function C = CallBS(S,K,r,T,sigma,q)
% Computes the value of the Call using BS formula

d1 = (log(S/K)+(r-q+sigma^2/2)*T)/(sigma*sqrt(T));
d2 = d1-sigma*sqrt(T);

C =  S*exp(-q*T)*normcdf(d1) - K*exp(-r*T)*normcdf(d2);
end

function tridiag = matA(sigma,dT,dK,len)
Z=0.5*sigma.^2*dT/(dK)^2;%dim(sigma)=(len(K)-2,1)
n=len;
nones = ones(n, 1);
tridiag = diag(Z)*(2*diag(nones, 0) - diag(nones(1:end-1), -1) - diag(nones(1:end-1), 1))+diag(nones, 0);
tridiag(1,:) = [1,zeros(1,n-1)];%replace first and last row by (1,0,...,0) and (0,...,0,1)
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
%interpolation of sigmas between the sigmas that will be estimated  
sigma=interp1(positions,sigma_est,[1:len]','nearest');
sigma(1:positions(1))=sigma(positions(1)); %sigmas before the first estimated sigma
sigma(positions(end):end)=sigma(positions(end)); %sigmas after the last estimated sigma
end


function [impv,fval]=implvBS(S,K,r,T,q,C,sig0)
impv=zeros(size(C));
for i = 1:length(K)
    for j = 1:length(T)
        fct = @(sigma) (CallBS(S,K(i),r,T(j),sigma,q) - C(i,j))^2;
        [impv(i,j),fval(i,j)] = fminsearch(fct,sig0);
    end
end
end

