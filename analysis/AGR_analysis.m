% compute the apparentl guessing rate given predictions from other models
% using the ML fits from real data

function AGR_analysis(skip_analysis)

% number of times to run the analysis (to eliminate bias)
total_runs = 10;
load subjCell

% revert to old data format
subjCell = revert_data(subjCell);

load parameter_ranges_new
set_sizes = unique(subjCell{1}(:,5))';

% get the actual ML parameters
[A_SA B_SA C_SA D_SA] = get_max_params(4,0);
[A_SR B_SR C_SR D_SR] = get_max_params(17,0);
[A_VP B_VP C_VP D_VP] = get_max_params(15,0);
[A_IL B_IL C_IL D_IL] = get_max_params(5,0);
[A_EP B_EP C_EP D_EP] = get_max_params(16,0);

% number of times to replicate each subject's data (to eliminate bias)
reps = 25;

% if not already generated the fake data, do it
if ~skip_analysis
for ii = 1:total_runs
    
    % for each subject, generate fake data using their ML parameters.
    % Uncomment the model currently being tested
    for jj = 1:length(subjCell)
%             subjFakeCell{jj} = make_fake_SA(repmat(subjCell{jj},reps,1),slotJvec(A_SA(jj)),kvec(B_SA(jj)),priorvec(C_SA(jj)));
%             subjFakeCell{jj} = make_fake_SR(repmat(subjCell{jj},reps,1),meanvec(A_SR(jj)),kvec(B_SR(jj)),priorvec(C_SR(jj)));
            subjFakeCell{jj} = make_fake_VP_new(repmat(subjCell{jj},reps,1),meanvecVR(A_VP(jj)),thetavec(B_VP(jj)),alphavec(C_VP(jj)),priorvec(D_VP(jj)));
%             subjFakeCell{jj} = make_fake_IL(repmat(subjCell{jj},reps,1),kvec(A_IL(jj)),lapsevec(B_IL(jj)),qvec(C_IL(jj)));
%             subjFakeCell{jj} = make_fake_EP(repmat(subjCell{jj},reps,1),meanvec(A_EP(jj)),alphavec(B_EP(jj)),0.5);
%             subjFakeCell{jj} = make_fake_EG(repmat(subjCell{jj},reps,1),meanvec(A_EP(jj)),alphavec(B_EP(jj)),0.5);
            
    end
    ii
    save(['subjFakeCell' num2str(ii)],'subjFakeCell')
end

% compute guessing rate as determined by fitting EP + set size independent 
% guessing model
all_guess = zeros(total_runs,length(subjCell),4);
all_J = zeros(total_runs,length(subjCell),4);
all_prior = zeros(total_runs,length(subjCell));
for ii = 1:total_runs
    
    [a1 b1 c1 d1] = get_max_params_cluster(7,2,ii);
   
    all_guess(ii,:,:) = lapsevec(b1);
    all_J(ii,:,:) = meanvecNPER(a1);
    all_prior(ii,:,:) = priorvec(c1);
end
save all_guess all_guess all_J all_prior


end

% plot results
load AGR_all
AGR = AGR_all;
figure;
set(gcf,'Units','inches','OuterPosition',[4 4 9.78 6.15/2])
ctr = .5;

% function defined below
draw_AGR_plot(1,'IL',AGR,set_sizes,ctr)
ctr = ctr+.75*1.81;
draw_AGR_plot(2,'SA',AGR,set_sizes,ctr)
ctr = ctr+.75*1.81;
draw_AGR_plot(3,'SR',AGR,set_sizes,ctr)
ctr = ctr+.75*1.81;
draw_AGR_plot(4,'EP',AGR,set_sizes,ctr)
ctr = ctr+.75*1.81;
draw_AGR_plot(5,'VP',AGR,set_sizes,ctr)

end

function draw_AGR_plot(model_idx,model_id,AGR,set_sizes,ctr)
[a1 b1 c1 d1] = get_max_params(7,0);
load parameter_ranges_new

fig_size = get(gcf,'Position');
subplot('Position',[ctr/fig_size(3) .5/fig_size(4) .75*1.44/fig_size(3) .75*1.3/fig_size(4)])

hold on
AGR_tmp = squeeze(mean(AGR(model_idx,:,:,:),2));

lapse = mean(AGR_tmp,1);
err = std(AGR_tmp,1)/sqrt(size(AGR,3));
polyX = [set_sizes set_sizes(end:-1:1)];
polyY = zeros(size(polyX));
polyY(1:length(set_sizes)) = lapse-err;
polyY((2*length(set_sizes)):-1:(length(set_sizes)+1)) = lapse+err;
fill(polyX,polyY,[204 204 204]/255,'LineStyle','none');
errorbar(set_sizes,mean(lapsevec(b1),1),std(lapsevec(b1),[],1)/sqrt(size(b1,1)),'ko','LineWidth',1,'MarkerSize',5)
b1_tmp = lapsevec(b1);
rmse = sqrt(mean((AGR_tmp(:)-b1_tmp(:)).^2));
title(model_id)
ylim([0 1])
set(gca,'XTick',set_sizes,'YTick',0:.2:1)
xlabel('Set size')
text(2,.9,['RMSE=' num2str(rmse,'%.2f')])
set(gca,'FontName','Arial','FontSize',11)

end