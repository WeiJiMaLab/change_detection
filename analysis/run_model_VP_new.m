function [L P_C_HAT max_P_C_HAT] = run_model_VP_new(Data,meanvec,thetavec,alphavec,priorvec,samples)

% Set up bessel lookup function
bessel_break_point = 1e-3;
k_break_point = 1e-3;
Lookup_small = linspace(0,bessel_break_point,5000000)';
LookupY_small = besseli(0,Lookup_small(:,1),1);
lookup_spacing_small = 1/(Lookup_small(2)- Lookup_small(1));

Lookup = linspace(0,5000,5000000)';
Highest_J = 5000;
LookupY = besseli(0,Lookup(:,1),1);
lookup_spacing = 1/(Lookup(2)- Lookup(1));

% get priors
logprior = log(priorvec./(1-priorvec));

% initialize output vectors
L = zeros(length(meanvec),length(thetavec),length(alphavec),length(priorvec));
P_C_HAT = zeros([size(L) size(Data,1)]);

% precompute inverse Fisher information function (inverse J)
k = linspace(0,5000,10000000)';
J = k.*(besseli(1,k,1)./besseli(0,k,1));
J_new = linspace(min(J),max(J),100000)';
k_new = interp1(J,k,J_new);
J_spacing = 1/(J_new(2)-J_new(1));
J_max_high = max(J);

k_small = linspace(0,k_break_point,10000000)';
J_small = k_small.*(besseli(1,k_small,1)./besseli(0,k_small,1));
J_new_small = linspace(min(J_small),max(J_small),1000000)';
k_new_small = interp1(J_small,k_small,J_new_small);
J_spacing_small = 1/(J_new_small(2)-J_new_small(1));
J_break_point = max(J_new_small);
nSamples = samples;

% get set sizes
SetSizes = unique(Data(:,5));

for m = 1:length(SetSizes)
    
    N = SetSizes(m);
    
    % Initialize priors and Likelihood matrix
    curr_N_idx = find(Data(:,5)==N);
    nTrials = length(curr_N_idx);
    phi = pi/90*Data(curr_N_idx,56:(55+N));
    theta = pi/90*Data(curr_N_idx,64:(63+N));
        
    % keep track of NaNs to remove them
    nancount = zeros(1,length(meanvec));
    
    time_idx = length(meanvec)*length(thetavec)+1;
    
    % loop over mean J
    for jind = 1:length(meanvec)
        
        % loop over theta parameter
        for tind = 1:length(thetavec)
            
            time_idx = time_idx - 1;
            tic;
            
            % loop over alpha parameter
            for aind = 1:length(alphavec)
                
                % draw J's and compute corresponding kappa's (may take NaN value..)
                curr_alpha = alphavec(aind);
                J_bar = meanvec(jind)*(N^(-curr_alpha));
                J_sample = gamrnd(J_bar/thetavec(tind),thetavec(tind),[nSamples 2*N]);
                J_sample(J_sample>J_max_high) = J_max_high;
                kappa = zeros(size(J_sample));
                J_high = J_sample>J_break_point;
                J_low = J_sample<J_break_point;
                kappa(J_high) = myBessel(J_sample(J_high),J_spacing,k_new); % not actually bessel
                kappa(J_low) = myBessel(J_sample(J_low),J_spacing_small,k_new_small); % not actually bessel
                                
                % resample NaNs
                nanidx = find(isnan(kappa));
                
                % while there are NaNs
                while ~isempty(nanidx)
                    nancount(jind) = nancount(jind) + length(nanidx);
                    J_sample = gamrnd(J_bar/thetavec(tind),thetavec(tind),1,length(nanidx));
                    kappa(nanidx)  = interp1(J',k',J_sample);
                    nanidx = find(isnan(kappa));
                    fprintf('Number of NaN: %d\n',nancount(jind));
                end
                
                % prepare kappas for decision rule
                kappa1 = kappa(:,1:N);
                kappa2 = kappa(:,(N+1):(2*N));

                noise1 = circ_vmrnd(0,kappa1,1);
                noise2 = circ_vmrnd(0,kappa2,1);

                kappa1_sq = kappa1.^2;
                kappa2_sq = kappa2.^2;
                                
                kappa1_bessel_log = log(besseli(0,kappa1,1));
                kappa2_bessel_log = log(besseli(0,kappa2,1));
                
                % for each trial
                for i=1:nTrials
                    
                    % recycle noise and kappas, compute decision rule
                    noise1 = noise1([2:end 1],:);
                    kappa1 = kappa1([2:end 1],:);
                    kappa1_sq = kappa1_sq([2:end 1],:);
                    kappa1_bessel_log = kappa1_bessel_log([2:end 1],:);
                    
                    costerm = cos(bsxfun(@plus,phi(i,:)-theta(i,:),noise1-noise2));
                    K3 = sqrt(kappa1_sq + kappa2_sq + 2*kappa1.*kappa2.*costerm);
                    if ~isreal(K3)
                        fprintf('\nimaginary!\n');
                    end
                    
                    K3(K3>Highest_J) = Highest_J;
                    K3_bessel = myBessel(K3,lookup_spacing,LookupY);
                    K3_bessel(K3<bessel_break_point) = myBessel(K3(K3<bessel_break_point),lookup_spacing_small,LookupY_small);
                    d = log(row_sum(exp(kappa1_bessel_log+kappa2_bessel_log-log(K3_bessel)+kappa1+kappa2-K3)))-log(N);
                    
                    d = bsxfun(@plus,d,logprior);
                    
                    % compute p(C_hat = 1 | x,y) and save
                    pC_hat = sum(d>0)/nSamples;
                    
                    P_C_HAT(jind,tind,aind,:,curr_N_idx(i)) = shiftdim(pC_hat,-1);
                end
            end
            
            fprintf('N=%g, J1=%2.2f, theta=%2.2f, ETL=%2.2f minutes\n' ,N,meanvec(jind),thetavec(tind),time_idx*toc/60);
            
        end
        
    end
    
end

% old versions returned max log-likelihood predictions, now just a filler
max_P_C_HAT = Data;
fprintf('\nDone!\n')