function [theta_new,acc] = MH_educ(Y,X,N,K,theta_old,sigma)
%Draw Proposed Distribution From Gaussian (MV Normal)
theta_prop= mvnrnd(theta_old,sigma);

%Likelihood Function of Proposed Using CDF of Normal F(Y|Theta)
Y_hat = X*theta_prop(1:K)';
f_i=zeros(N,1);

for j=1:N
    if theta_prop(K+1)<=0
        f_i(j)=0;
    else
        f_i(j)= normpdf(Y(j),Y_hat(j),sqrt(theta_prop(K+1)));
    end 
end

% Likelihood Function is the Product of individual pdfs
% Log Likelihood is the sum of individual log pdfs
% Attach The Education Prior Distribution
educ_prior= normpdf(theta_prop(2),.06,.007);

LL_prop= sum(log(f_i))+log(educ_prior);

%Likelihood Function Using CDF of Normal F(Y|Theta)
Y_hat=X*theta_old(1:K)';
f_i=zeros(N,1);

for j=1:N
    if theta_old(K+1)<=0
        f_i(j)=0;
    else
        f_i(j)= normpdf(Y(j),Y_hat(j),sqrt(theta_old(K+1)));
    end 
end

% Attach The Education Prior Distribution
educ_prior= normpdf(theta_old(2),.06,.007);

LL_old= sum(log(f_i));

Acc_prob= exp(LL_prop -LL_old);

thresh= rand;

if Acc_prob >= thresh
    theta_new = theta_prop;
    acc=1;
else 
    theta_new= theta_old;
    acc=0;
end
end