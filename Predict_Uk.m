function output=Predict_Uk(age,uk,pstd)
%
%BAYSPLINE: a Bayesian b-spline calibration for the alkenone
%paleothermometer. Predicts values of SST from measurements of uk37'.
%Please cite the source publication when using this calibration:
%
%Tierney, JE and Tingley, MP (2018). BAYSPLINE: A new calibration for the
%alkenone paleothermometer. Paleoceanography and Paleoclimatology, 33. https://doi.org/10.1002/2017PA003201
%
%INPUTS:
%age = vector of depth or age values of length N
%uk = uk37' values of length N
%pstd = prior standard deviation (scalar). Recommended value: 10. This 
%conservative value that works for most applications.
%At high UK, a more restrictive value of 5 may be desirable to
%assign a low probability to SSTs > 35C (Unlikely given that UK approaches 1 near this value). 

%OUTPUTS:
%output.prior_mean = Prior mean value, taken from the UK timeseries
%converted to SST with the Prahl '88 culture equation.
%
%output.prior_std = Prior standard deviation (user set).
%
%output.jump_dist = standard deviation of the jumping distribution.
%Values are chosen to achieve a acceptance rate of ca. 0.44
%(Gelman, 2003).
%
%output.rhat = Rhat statistic to assess convergence. Should be close to 1.
%
%output.SST = 5 x N vector of inferred SSTs, includes 5% level (lower 2sigma), 16% level
%(lower 1sigma), 50% level (median values), 84% level (upper 1sigma), and
%95% level (upper 2 sigma).
%
%output.ens = full ensemble (N = 2500) of posterior SSTs.
%
%PLOTS:
%1) Prior distribution vs posterior distribution. Prior should be wider
%unless a strong control is warranted.
%
%2) Time series converted to SST values, with 1-sigma confidence intervals.
%
%For a full explanation of the approach see Tierney & Tingley (2018).
%% load model parameters
load bayes_posterior_v2.mat;

%thin the posterior draws a bit
b_draws_final=b_draws_final(1:3:end,:);
tau2_draws_final=tau2_draws_final(1:3:end);

%confirm UK obs are column vector
uk=uk(:);

N_Ts=length(uk);
N_p=length(tau2_draws_final);
%number of samples
N=500;
%burnin to discard
burnin=250;

%set priors. Use Prahl 88 to create a prior mean vector
prior_mean=(uk-.039)./.034;
%prior variance comes from user-defined pstd
prior_var=pstd^2*ones(N_Ts,1);
%save priors to output
output.prior_mean=prior_mean;
output.prior_std=pstd;

%set the initial SST values to prior mean
init=prior_mean;

%create empty matrices of the correct size
MH_samps_t=NaN(N_Ts,N_p,N-burnin);
accepts_t=NaN(N_Ts,N_p,N-burnin);

%make a spline with set knots
order=3; %spline order
kn = augknt(knots,order); %knots

%assign width of jumping distribution based on average UK temp to obtain
%40-50% accept rate
pmm=median(prior_mean);
if pmm<20
    JW=3.5;
    elseif pmm>=20 && pmm <= 23.7
        JW=3.7;
    else
        JW=pmm*0.8092-15.1405;
end
output.jump_dist=JW;
%% MH loop
tic
parfor jj=1:N_p
    accepts = NaN(N_Ts,N);
    samps = NaN(N_Ts,N);
    b_now=b_draws_final(jj,:);
    tau_now=tau2_draws_final(jj);
    %use spmak to put together the b-spline
    bs_b=spmak(kn,b_now);
    %extrapolate function
    bs=fnxtr(bs_b);
    %intialize at starting value
    samps(:,1)=init;
    s_now=samps(:,1);
    %evaluate mean UK value at current SST
    mean_now=fnval(bs,s_now);
    %evaluate likelihood
    LL_now=normpdf(uk,mean_now,sqrt(tau_now));
    %evaluate prior
    pr_now=normpdf(s_now,prior_mean,sqrt(prior_var));
    %multiply to get initial proposal S0
    S0_now=LL_now.*pr_now;
    for kk=2:1:N
        %generate proposal using normal jumping distr.
        s_prop=normrnd(s_now,JW);
        %evaluate mean value at current SST
        mean_now=fnval(bs,s_prop);
           %evaluate likelihood
            LL_now=normpdf(uk,mean_now,sqrt(tau_now));
            %evaluate prior
            pr_now=normpdf(s_prop,prior_mean,sqrt(prior_var));
            %multiply to get proposal S0_p
            S0_p=LL_now.*pr_now;
            %calculate MH ratio
            MH_rat=S0_p./S0_now;
            success_rate=min(1, MH_rat);

        %make the draw:
            draw=rand(N_Ts,1);
            B = draw <= success_rate;
            s_now(B)=s_prop(B);
            S0_now(B)=S0_p(B);

            accepts(B,kk)=1;
            samps(:,kk)=s_now;
    end
    MH_samps_t(:,jj,:)=samps(:,burnin+1:end);
    accepts_t(:,jj,:)=accepts(:,burnin+1:end);
end
toc
%now let's calculate the Rhat statistic to assess convergence.
rhats=NaN(size(MH_samps_t,1),1);
for i=1:size(MH_samps_t,1)
    [rhats(i),~]=ChainConvergence(squeeze(MH_samps_t(i,:,:)), N_p);
end
output.rhat=median(rhats);
%reshape
MH_c=reshape(MH_samps_t,N_Ts,N_p*(N-burnin));
%calculate acceptance
output.accepts = nansum(accepts_t(:))./(N_Ts*N_p*(N-burnin));

%sort and assign to output
MH_s=sort(MH_c,2);
pers5=round([.05 .16 .5 .84 .95].*size(MH_c,2));
output.SST=MH_s(:,pers5);
%%
%take a subsample of MH to work with for ksdensity.
MH_sub=MH_c(:,1:50:end);
output.ens=MH_sub;
%plot prior vs post
f1=figure(1); clf;

set(f1,'pos',[50 700 400 400]);
xt=[0:.1:40]';
prior=normpdf(xt,prior_mean(1),sqrt(prior_var(1)));
post=ksdensity(MH_sub(:),xt);
pr=plot(xt,prior,'k--'); hold on;
pt=plot(xt,post,'b-');
legend([pr pt],'Prior','Posterior');
legend('boxoff');

%%
%plot timeseries with 1-sigma errors
f2=figure(2); clf;
set(f2,'pos',[550 700 500 400]);
 plot(age,output.SST(:,3),'color','k','linewidth',2);
 hold on;
 plot(age,output.SST(:,2),'color',[.4 .4 .4],'linewidth',1,'linestyle','--');
  hold on;
 plot(age,output.SST(:,4),'color',[.4 .4 .4],'linewidth',1,'linestyle','--');
 set(gca,'xlim',[min(age) max(age)]);
function [Rhat, Neff]=ChainConvergence(chains, M)
%
% function [Rhat, Neff]=ChainConvergence(chains, M)
% 
% calculate the R-hat statistic for
% monitoring  convergence of multiple parallel MCMC runs. Also outputs
% n_eff.. 
% chains: a matrix of MCMC chains, each of the same length. 
% M: the number of different chains - must be one of the dimensions of
% chains. 
%


if ismember(M, size(chains))==0
    disp('Error: Second input must be the number of chains, so one of the dimensions of the first input.')
    return
elseif length(chains(:,1))==M
    chains=chains';
end

% each column is a chain. 
N=length(chains(:,1));

%quantities from Gelman:

psi_bar_dot_j = mean(chains);
psi_bar_dot_dot=mean(psi_bar_dot_j);

B= (N/(M-1))*sum((psi_bar_dot_j-psi_bar_dot_dot).^2);

s2_j= (1/(N-1))*sum((chains-kron(psi_bar_dot_j, ones(N,1))).^2);

W=(1/M)*sum(s2_j);

var_hat_pos=((N-1)/N)*W + (1/N)*B;

Rhat=sqrt(var_hat_pos/W);
Neff=M*N*min(var_hat_pos/B, 1);