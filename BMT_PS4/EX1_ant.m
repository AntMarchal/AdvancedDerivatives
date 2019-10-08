data =xlsread('Impvols_SPX_AMZN.xlsx');
%_1->SPX
%_2->AMZN
S_1=2921;
S_2=1971;
T=0.296;
r=0.024;
q_1=0.018;
q_2=0.019;
rho=0.5;

K_1 = data(:,1);
IV_1 = data(:,2);
%remove the nans
K_1 = K_1(~isnan(K_1));
IV_1 = IV_1(~isnan(K_1));

K_2 = data(:,5);
IV_2 = data(:,6);



%Step1: compute undiscounted call prices
C_1_=zeros(length(K_1),1);% the _ indicates that this is undiscount price
C_2_=zeros(length(K_2),1);
for i=1:length(K_1)
    C_1_(i,1)=exp(r*T)*BSCall(S_1,K_1(i,1),r,T,IV_1(i,1),q_1);
end
for i=1:length(K_2)
    C_2_(i,1)=exp(r*T)*BSCall(S_2,K_2(i,1),r,T,IV_2(i,1),q_2);
end  
%step2: get the implied CDF
CDF1=1+(C_1(2:end)-C_1(1:end-1))./(K_1(2:end)-K_1(1:end-1));
CDF2=1+(C_2(2:end)-C_2(1:end-1))./(K_2(2:end)-K_2(1:end-1));
%Step3: generate 10000 mvnrnd and obtain random CDF
mu=[0,0];
Cov=[1, rho; rho,1];
X=normcdf(mvnrnd(mu,Cov,10000));
%step4: Use inverse CDF to get 10000 simulated S_1T, S_2T

%interp1(x,v,xq) returns interpolated values of a 1-D function at specific query points using 
%linear interpolation. Vector x contains the sample points, and v contains the corresponding
%values, v(x). Vector xq contains the coordinates of the query points.
    
S_1T=interp1(CDF1,K_1(1:end-1),X(:,1));
S_2T=interp1(CDF2,K_2(1:end-1),X(:,2));


%step5: compute payoff
A=S_1T/S_1-S_2T/S_2;
%remove the nans
A=A(~isnan(A));
psi_T=max(A,0);

%step6: discounted average

price=mean(psi_T)*exp(-r*T)












function c = BSCall(S,K,r,T,sigma,q)
% Computes the value of the Call using BS formula

d1 = (log(S/K)+(r-q+sigma^2/2)*T)/(sigma*sqrt(T));
d2 = d1-sigma*sqrt(T);

c =  S*exp(-q*T)*normcdf(d1) - K*exp(-r*T)*normcdf(d2);
end
