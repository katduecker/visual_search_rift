%% Simple simulation: power, time on task (tot) and reaction time (rt)
% Power depends on tot and on rt. Two cases:
%   case 1: rt is correlated with tot
%   case 2: rt is not correlated with tot
%   case 3: power and rt are both correlated with tot, but rt has no
%           effect on power
%
% The noise in rt and in power is autocorrelated over trials (AR(1)).
%
% A GLM is fitted to power (power ~ rt + tot) and the T-value of rt is
% tested against a null distribution, obtained by shuffling power over
% blocks (as in c1_GLM_spectrum_interactions.m)

% 9 Oct 2026

clear all; close all; clc;
%% settings

n_blocks = 20;      % number of blocks
block_len = 40;     % trials per block
n_trl = n_blocks*block_len;     % number of trials

num_perm = 500;     % number of permutations

b_tot = 1;          % effect of tot on power
b_rt = -.5;          % effect of rt on power
r_rt_tot = .6;      % correlation between rt and tot in case 1

ar_rt = .9;         % autocorrelation (lag 1) of the noise in rt; 0: independent trials
ar_pow = .9;        % autocorrelation (lag 1) of the noise in power; 0: independent trials

% figure 1
c_plot = 3;         % case shown in the time course and autocorrelation panels
max_lag = 200;      % maximum lag (in trials) in the autocorrelation panel

% figure 2
b_rt_sim = [0, -.05, -.25, -.5];    % effects of rt on power that are simulated
n_sim = 1000;                       % simulated datasets per effect

rng(80)                           % this seed gives a spurious correlation between RT and pow in MATLAB 2026a when not accounting for null-distribution (for illustration, overall effect quantified below)

%% time on task: straight line from 0 to 1

tot = linspace(0,1,n_trl)';

% z-scored copy, so that b_tot and b_rt are on the same scale
tot_z = (tot - mean(tot))/std(tot);

%% autocorrelated noise

% autocorrelated series: each value is ar (auto-correlation) * previous trial + new random value
% (achieved through filtering), the value is normalized by sqrt(1-ar^2) so
% that original white noise and filtered (ar-) noise have the same range.

white_noise = randn(n_trl,1);                    % random values
[noise_rt,zf] = filter(1,[1 -ar_rt], white_noise);      % filtered noise (ar)
noise_rt_norm = noise_rt*sqrt(1-ar_rt^2);               % normalize

[auto_white,lags_aw] = xcorr(white_noise);
[auto_rt,lags_arcn] = xcorr(noise_rt);
[auto_rt_norm,lags_arcn2] = xcorr(noise_rt_norm);


noise_pow = filter(1,[1 -ar_pow],randn(n_trl,1))*sqrt(1-ar_pow^2);

%% reaction time

% case 1: RT correlated with tot
rt_corr = r_rt_tot*(-1)*tot_z + sqrt(1-r_rt_tot^2)*noise_rt_norm;

% case 2: RT not correlated with tot
rt_uncorr = noise_rt_norm;

%% power: tot + rt + noise

pow_corr   = b_tot*tot_z + b_rt*rt_corr   + noise_pow;          % power correlated with RT and tot, and RT is correlated with tot
pow_uncorr = b_tot*tot_z + b_rt*rt_uncorr + noise_pow;          % power correlated with RT and tot, and RT is not correlated with tot
pow_uncorr_rt = b_tot*tot_z + noise_pow;                        % power only correlated with tot (to GLM design matrix, we add the rt that is also correlated with tot)

%% check correlations

names = {'1: pow x rt+tot, rt x tot','2: pow x rt+tot, rt ~x tot', '3: pow x tot, rt x tot'};
rt_all = {rt_corr, rt_uncorr, rt_corr};                         % RTs used for design matrix
pow_all = {pow_corr, pow_uncorr, pow_uncorr_rt};                % power values that are the Y's for the model
b_rt_true = [b_rt, b_rt, 0];    % true effect of rt on power in each case

% %% plot time series
% 
% figure('Position',[50 50 1000 550])
% for c = 1:3
%     subplot(3,1,c); hold on
%     plot(rt_all{c},'Color',[.1 .45 .7])
%     plot(pow_all{c},'Color',[.85 .33 .1])
%     xlabel('trial'); title(names{c})
% end
% legend({'rt','power'},'Location','best')

%% GLM and block-shuffle of power
num_sim = length(rt_all);
T_rt = zeros(1,num_sim);                  % observed T-value of rt
T_perm_rt = zeros(num_perm,num_sim);      % null distribution
z_score_T = zeros(1,num_sim);
p_perm = zeros(1,num_sim);

% block boundaries
block_boundaries = (1:block_len:n_trl+1)';
num_blocks = length(block_boundaries)-1;

for c = 1:num_sim

    % design matrix
    X = ones(n_trl,3);
    % rt
    X(:,2) = rt_all{c};
    % tot
    X(:,3) = tot - mean(tot);

    Y = pow_all{c};

    cov_forming_matrix = pinv(X'*X);    % covariance forming matrix
    var_forming_matrix = diag(cov_forming_matrix);
    Xp = pinv(X);
    dof_error = n_trl - rank(X);

    % fit model
    model_beta = Xp*Y;
    RSS = sum((Y - X*model_beta).^2);
    varres = RSS/dof_error;
    varcope = var_forming_matrix*varres;
    model_T = model_beta./sqrt(varcope);

    T_rt(c) = model_T(2);

    % null distribution: shuffle power over blocks
    for n = 1:num_perm

        shuffled_order = randperm(num_blocks);

        shuf_Y = zeros(size(Y));
        cnt = 1;
        for i = 1:num_blocks
            start_index = block_boundaries(shuffled_order(i));
            end_index = block_boundaries(shuffled_order(i)+1)-1;
            block = Y(start_index:end_index);
            shuf_Y(cnt:cnt+(end_index-start_index)) = block;
            cnt = cnt+(end_index-start_index)+1;
        end

        % fit model
        model_beta_perm = Xp*shuf_Y;
        RSS = sum((shuf_Y - X*model_beta_perm).^2);
        varres = RSS/dof_error;
        varcope = var_forming_matrix*varres;
        model_T_perm = model_beta_perm./sqrt(varcope);

        T_perm_rt(n,c) = model_T_perm(2);
    end

    z_score_T(c) = (T_rt(c) - mean(T_perm_rt(:,c)))/std(T_perm_rt(:,c));
    p_perm(c) = (sum(abs(T_perm_rt(:,c)) >= abs(T_rt(c))) + 1)/(num_perm + 1);

    fprintf('\n%s\n', names{c})
    fprintf('  beta rt = %5.2f (true: %.2f), T = %5.2f\n', model_beta(2), b_rt_true(c), T_rt(c))
    fprintf('  null (shuffle power over blocks): mean = %5.2f, SD = %5.2f\n', mean(T_perm_rt(:,c)), std(T_perm_rt(:,c)))
    fprintf('  z = %5.2f, p = %.3f\n', z_score_T(c), p_perm(c))
end

%% Figure 1: the current example

% figure settings
col_rt  = [.1 .45 .7];
col_pow = [.85 .33 .1];
col_tot = [0 0 0];
col_case = [[.2 .6 .3];[.5 .3 .7];[.9 .6 0]];       % one colour per case

% Student's T distribution with the degrees of freedom of the GLM
t_ax = linspace(-15,15,1000);
% using exp and gammaln to avoid overflow of gamma function beyond
% gamma(171) due to high DOF
t_pdf = exp(gammaln((dof_error+1)/2) - gammaln(dof_error/2))/sqrt(dof_error*pi)*(1+t_ax.^2/dof_error).^(-(dof_error+1)/2);

fig = figure('Position',[50 50 1100 800]);
theme(fig, "light")
% top left: time course over trials for power, RT, TOT
subplot(2,2,1); hold on
plot(1:n_trl,pow_all{c_plot},'Color',col_pow)
plot(1:n_trl,rt_all{c_plot},'Color',col_rt)
plot(1:n_trl,tot_z,'Color',col_tot,'LineWidth',2)
xlim([1 n_trl])
xlabel('trial'); ylabel('a.u.')
legend({'pow','rt','tot'},'Location','northwest')
title(names{c_plot})

% top right: autocorrelation of rt and pow
[auto_rt_plot,lags_plot] = xcorr(rt_all{c_plot} - mean(rt_all{c_plot}),max_lag,'coeff');
auto_pow_plot = xcorr(pow_all{c_plot} - mean(pow_all{c_plot}),max_lag,'coeff');
auto_white_plot = xcorr(white_noise - mean(white_noise),max_lag,'coeff');

subplot(2,2,2); hold on
plot(lags_plot,auto_white_plot,'Color',[.6 .6 .6])
plot(lags_plot,auto_pow_plot,'Color',col_pow,'LineWidth',2)
plot(lags_plot,auto_rt_plot,'Color',col_rt,'LineWidth',2)
plot([block_len block_len],[-.2 1],'k--')
xlim([0 max_lag]); ylim([-.2 1])
xlabel('lag (trials)'); ylabel('autocorrelation')
legend({'white noise','pow','rt','block length'},'Location','northeast')
title('autocorrelation')

% bottom left: observed T-values against Student's T distribution
subplot(2,2,3); hold on
plot(t_ax,t_pdf,'k','LineWidth',2)
for c = 1:num_sim
    plot([T_rt(c) T_rt(c)],[0 max(t_pdf)],'Color',col_case(c,:),'LineWidth',2)
end
plot([-1.96 -1.96],[0 max(t_pdf)],'k:')
plot([1.96 1.96],[0 max(t_pdf)],'k:')
xlim([min([T_rt,-5])-1, max([T_rt,5])+1])
xlabel('T-value of rt'); ylabel('probability density')
leg = [{'Student''s T'}, names];
legend(leg,'Location','best')
title('T vs Student''s T')

% bottom right: observed T-values against the null distribution (block-shuffle)
subplot(2,2,4); hold on
plot(t_ax,t_pdf,'Color',[.6 .6 .6],'LineWidth',1)
bin_ctr = linspace(min(T_perm_rt(:)),max(T_perm_rt(:)),50);
for c = 1:num_sim
    cnt_null = hist(T_perm_rt(:,c),bin_ctr);
    dens_null = cnt_null/(sum(cnt_null)*(bin_ctr(2)-bin_ctr(1)));
    [x_st,y_st] = stairs(bin_ctr-(bin_ctr(2)-bin_ctr(1))/2,dens_null);
    plot(x_st,y_st,'Color',col_case(c,:),'LineWidth',1.5)
end
for c = 1:num_sim
    plot([T_rt(c) T_rt(c)],[0 max(t_pdf)],'Color',col_case(c,:),'LineWidth',2)
end
xlim([min([T_rt,-5])-1, max([T_rt,5])+1])
xlabel('T-value of rt'); ylabel('probability density')
leg = {'Student''s T'};
for c = 1:num_sim
    leg{end+1} = sprintf('T & null %d', c);
end
legend(leg,'Location','best')
title('T vs empirical null')

exportgraphics(gcf, "distributions.pdf", ContentType='vector')

%% Figure 2: many simulations for different effects of rt on power

% For each effect of rt on power (b_rt_sim), n_sim datasets with new noise
% in rt and power are simulated. rt is correlated with tot (as in cases 1
% and 3); b_rt = 0 is case 3. For each dataset:
%   T_sim:      T-value of rt from the GLM
%   z_t_sim:    T converted to z using Student's T distribution
%   z_null_sim: T z-scored against the null distribution (block-shuffle of power)

n_b = length(b_rt_sim);
T_sim = zeros(n_sim,n_b);
z_null_sim = zeros(n_sim,n_b);

tic
for ib = 1:n_b
    disp(['b_rt = ', num2str(b_rt_sim(ib))])

    for is = 1:n_sim

        % simulate
        noise_rt_sim = filter(1,[1 -ar_rt],randn(n_trl,1))*sqrt(1-ar_rt^2);
        noise_pow_sim = filter(1,[1 -ar_pow],randn(n_trl,1))*sqrt(1-ar_pow^2);

        rt_sim = r_rt_tot*tot_z + sqrt(1-r_rt_tot^2)*noise_rt_sim;
        pow_sim = b_tot*tot_z + b_rt_sim(ib)*rt_sim + noise_pow_sim;

        % GLM
        X = ones(n_trl,3);
        X(:,2) = rt_sim;
        X(:,3) = tot - mean(tot);
        Y = pow_sim;

        var_forming_matrix = diag(pinv(X'*X));
        Xp = pinv(X);
        dof_error = n_trl - rank(X);

        model_beta = Xp*Y;
        RSS = sum((Y - X*model_beta).^2);
        T_sim(is,ib) = model_beta(2)/sqrt(var_forming_matrix(2)*RSS/dof_error);

        % null distribution: shuffle power over blocks. Same as loop above
        % in matrix format.
        [~,shuffled_order] = sort(rand(num_blocks,num_perm),1);
        block_start = reshape(block_boundaries(shuffled_order),1,num_blocks,num_perm);
        % permutation index matrix, with randomised values in each column
        perm_idx = repmat(block_start,block_len,1) + repmat((0:block_len-1)',[1 num_blocks num_perm]);
        perm_idx = reshape(perm_idx,n_trl,num_perm);

        shuf_Y = Y(perm_idx);               % shuffle Y
        beta_perm = Xp*shuf_Y;              % fit GLM to all shuffled Ys
        RSS_perm = sum((shuf_Y - X*beta_perm).^2,1);    % Residual sum of squares
        T_perm = beta_perm(2,:)./sqrt(var_forming_matrix(2)*RSS_perm/dof_error);    % permuted T-distribution

        z_null_sim(is,ib) = (T_sim(is,ib) - mean(T_perm))/std(T_perm);              % z-scores
    end
end
toc


% plot
col_b = [[.5 .5 .5];[.99 .75 .3];[.9 .4 .2];[.6 .1 .2];[.3 0 .3];[0 0 0]];
col_b = col_b(1:n_b,:);

dat_all = {T_sim, z_null_sim};
tit_all = {'GLM T vs Student','z(GLM T) vs normal'};
xlab_all = {'T-value of rt','z(T-value of rt)'};

% reference distributions: Student's T for T, standard normal for z
z_pdf = exp(-t_ax.^2/2)/sqrt(2*pi);
ref_all = {t_pdf, z_pdf};
ref_name = {'Student''s T','standard normal','standard normal'};

figure('Position',[50 50 1400 800])
fig = gcf;
theme(fig, "light")
for ip = 1:2
    dat = dat_all{ip};

    % top: distribution over simulations against the reference distribution
    subplot(2,2,ip); hold on
    lim = max(abs(dat(:)));
    bin_ctr = linspace(-lim,lim,60);
    bin_w = bin_ctr(2)-bin_ctr(1);
    for ib = 1:n_b
        cnt_sim = hist(dat(:,ib),bin_ctr);
        [x_st,y_st] = stairs(bin_ctr-bin_w/2,cnt_sim/(sum(cnt_sim)*bin_w));
        plot(x_st,y_st,'Color',col_b(ib,:),'LineWidth',1.5)
    end
    plot(t_ax,ref_all{ip},'k','LineWidth',2)
    xlim([-lim lim])
    xlabel(xlab_all{ip}); ylabel('probability density')
    title(tit_all{ip})
    if ip == 1
        leg = cell(1,n_b+1);
        for ib = 1:n_b
            leg{ib} = ['b_{rt} = ', num2str(b_rt_sim(ib))];
        end
        leg{end} = 'reference distribution';
        legend(leg,'Location','northwest')
    end

    % bottom: one dot per simulated dataset
    subplot(2,2,2+ip); hold on
    for ib = 1:n_b
        x_jit = ib + (rand(n_sim,1)-.5)*.6; % jitter on x-axis
        plot(x_jit,dat(:,ib),'.','Color',col_b(ib,:),'MarkerSize',4)     % T or z-values
        text(ib,lim*1.08,sprintf('%.0f%%',100*mean(dat(:,ib) < - 1.96)),'HorizontalAlignment','center')
    end
    plot([.4 n_b+.6],[1.96 1.96],'k--')
    plot([.4 n_b+.6],[-1.96 -1.96],'k--')
    xlim([.4 n_b+.6]); ylim([-lim lim*1.18])
    set(gca,'XTick',1:n_b,'XTickLabel',b_rt_sim)
    xlabel('b_{rt}'); ylabel(xlab_all{ip})
    title(['% of effects < -1.96 (', ref_name{ip}, ')'])
end
exportgraphics(gcf, "FDR.pdf", ContentType='vector')
