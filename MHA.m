% PROGRAM NAME: ps6MHA.m
clear, clc
rng(99890);

Data = xlsread('MH_data.xlsx');
[N,K]=size(Data); %Note: K includes the Dep Var

Y=Data(:,1);  %Dep Var
X=[ones(N,1) Data(:,2:end)]; %Indep Var and Constant

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% (#1) Estimate the OLS Parameters and Residuals %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 B=inv(X'*X)*X'*Y; % Beta Matrix
 Res= Y-(X*B); % Residuals
 sig_sq=(Res'*Res/(N-K-1)); % Sample Variance of Residuals
 var_sig= 2/(N-K-1)*(sig_sq^2); % Variance of Sigma Squared
 
 varB= sig_sq*inv(X'*X); 
 varB=diag(varB); % Var of Betas
 seB=sqrt(varB); % SE of Betas
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% (#2A) Metropolis-Hastings Algo using Flat Prior to Generate      %    
 %%% Posterior Distrubtions for Parameters                            %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 iteration=50000;
 burnin=1000;
 
 % Scale Variance of Betas and Sigma Sq. when choosing the proposed
 % Distribution
 r=0.08;
 
 Z=[varB ; var_sig]';
 Z=r*Z;
 
 theta = [B ; sig_sq]';
 
 % Posterior Storage
 Theta_flat= zeros(K+1, iteration);
 Theta_flat(:,1)= theta;
 acc_b_flat = 0;
 
 % Burnin
 
 for i=1:burnin
     [theta, acc] = MH_flat(Y,X,N,K,theta,Z);
     acc_b_flat = acc_b_flat + acc ;
 end 
 
 acc_rate_b= (acc_b_flat/burnin)*100
 
 % Continue
 
 acc_flat=0;
 for i=2:iteration
     [theta, acc]= MH_flat(Y,X,N,K,theta,Z);
     acc_flat = acc_flat + acc ;
     
     Theta_flat(:,i)= theta;
     
     i
 end
 
 acc_rate= (acc_flat/iteration)*100
 
 % Plots
 figure
 histfit(Theta_flat(1,:),50,'kernel')
 title(['\beta_{constant} with Mean= ' ,num2str(mean(Theta_flat(1,:))), ' and Variance= ',num2str(var(Theta_flat(1,:)))])
 
 figure
 histfit(Theta_flat(2,:),50,'kernel')
 title(['\beta_{educ} with Mean= ' ,num2str(mean(Theta_flat(2,:))), ' and Variance ',num2str(var(Theta_flat(2,:)))])
  
 figure
 histfit(Theta_flat(3,:),50,'kernel')
 title(['\beta_{exp} with Mean= ' ,num2str(mean(Theta_flat(3,:))), ' and Variance= ',num2str(var(Theta_flat(3,:)))])
    
 figure
 histfit(Theta_flat(4,:),50,'kernel')
 title(['\beta_{MSA} with Mean= ' ,num2str(mean(Theta_flat(4,:))), ' and Variance= ',num2str(var(Theta_flat(4,:)))])
  
 figure
 histfit(Theta_flat(5,:),50,'kernel')
 title(['\beta_{black} with Mean= ' ,num2str(mean(Theta_flat(5,:))), ' and Variance= ',num2str(var(Theta_flat(5,:)))])
 
 figure
 histfit(Theta_flat(6,:),50,'kernel')
 title(['\beta_{south} with Mean= ' ,num2str(mean(Theta_flat(6,:))), ' and Variance= ',num2str(var(Theta_flat(6,:)))])
 
 figure
 histfit(Theta_flat(7,:),50,'kernel')
 title(['\sigma^{2} with Mean= ' ,num2str(mean(Theta_flat(7,:))), ' and Variance= ',num2str(var(Theta_flat(7,:)))])
 
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% (#2B) Metropolis-Hastings Algo using a Prior on Education Parameters %
 %%% to Generate Posterior Distrubtions                                   %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 iteration=50000;
 burnin=1000;
 
 
 % Re-gen Theta Vector
 theta = [B ; sig_sq]';
 
 % Posterior Storage
 Theta_educ= zeros(K+1, iteration);
 Theta_educ(:,1)= theta;
 acc_b_educ = 0;
 
 % Burnin
 for i=1:burnin
     [theta, acc] = MH_educ(Y,X,N,K,theta,Z);
     acc_b_educ = acc_b_educ + acc ;
 end 
 
 acc_rate_b_educ= (acc_b_educ/burnin)*100
 
 % Continue
 
 acc_educ=0;
 for i=2:iteration
     [theta, acc]= MH_educ(Y,X,N,K,theta,Z);
     acc_educ = acc_educ + acc ;
     
     Theta_educ(:,i)= theta;
     
     i
 end
 
 acc_rate_educ= (acc_educ/iteration)*100
 
 % Plots
 figure
 histfit(Theta_educ(1,:),50,'kernel')
 title(['\beta_{constant} with Mean= ' ,num2str(mean(Theta_educ(1,:))), ' and Variance= ',num2str(var(Theta_educ(1,:)))])
 
 figure
 histfit(Theta_educ(2,:),50,'kernel')
 title(['\beta_{educ} with Mean= ' ,num2str(mean(Theta_educ(2,:))), ' and Variance ',num2str(var(Theta_educ(2,:)))])
  
 figure
 histfit(Theta_educ(3,:),50,'kernel')
 title(['\beta_{exp} with Mean= ' ,num2str(mean(Theta_educ(3,:))), ' and Variance= ',num2str(var(Theta_educ(3,:)))])
    
 figure
 histfit(Theta_educ(4,:),50,'kernel')
 title(['\beta_{MSA} with Mean= ' ,num2str(mean(Theta_educ(4,:))), ' and Variance= ',num2str(var(Theta_educ(4,:)))])
    
 figure
 histfit(Theta_educ(5,:),50,'kernel')
 title(['\beta_{black} with Mean= ' ,num2str(mean(Theta_educ(5,:))), ' and Variance= ',num2str(var(Theta_educ(5,:)))])
 
 figure
 histfit(Theta_educ(6,:),50,'kernel')
 title(['\beta_{south} with Mean= ' ,num2str(mean(Theta_educ(6,:))), ' and Variance= ',num2str(var(Theta_educ(6,:)))])
    
 figure
 histfit(Theta_educ(7,:),50,'kernel')
 title(['\sigma^{2} with Mean= ' ,num2str(mean(Theta_educ(7,:))), ' and Variance= ',num2str(var(Theta_educ(7,:)))])
    
    