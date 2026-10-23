function plot_bam_results(out)
% Redraw from saved MATLAB data; no simulation is repeated.
if nargin<1,out=fullfile(fileparts(mfilename('fullpath')),'..','outputs','reproduced_matlab');end
data=load(fullfile(out,'simulation.mat'));
t=data.t;zr=data.zr;trajectories=data.trajectories;errors=data.errors;p=data.p;
set(groot,'defaultFigureVisible','off');colors=lines(3);names={'x_1','x_2','y_1','y_2'};
fig=figure('Position',[100 100 1100 700]);tiledlayout(2,2,'TileSpacing','compact');
for j=1:4,nexttile;plot(t,zr(j,:),'LineWidth',1.1);xlabel('t');ylabel(names{j});grid on;end
finish(fig,out,'fig1_states');
fig=figure('Position',[100 100 1100 700]);tiledlayout(2,2,'TileSpacing','compact');
for j=1:4
 nexttile;hold on;hh=gobjects(1,4);
 for k=1:3,hh(k)=plot(t,trajectories(j,:,k),'Color',colors(k,:),'LineWidth',1.1);end
 hh(4)=plot(t,zr(j,:),'k--','LineWidth',1.2);
 xlim([0 15]);xlabel('t');ylabel(names{j});grid on;box on;
 if j==1,legend(hh,{'history 1','history 2','history 3','reference'},'Location','best');end
end
finish(fig,out,'fig2_histories');
fig=figure('Position',[100 100 1050 420]);tiledlayout(1,2,'TileSpacing','compact');
nexttile;hold on;
for k=1:3,plot(t,errors(k,:),'Color',colors(k,:),'LineWidth',1.1);end
xlim([0 15]);xlabel('t');ylabel('State sup error');grid on;box on;
nexttile;hold on;set(gca,'YScale','log');
for k=1:3,plot(t,max(errors(k,:),1e-16),'Color',colors(k,:),'LineWidth',1.1);end
plot(t,max(errors(:,1))*exp(-p.lambda*t),'k--','LineWidth',1.1);
xlim([0 25]);ylim([1e-14 10]);xlabel('t');ylabel('State sup error');grid on;box on;
set(gca,'YMinorGrid','off');
legend('history 1','history 2','history 3','slope guide e^{-0.2t}','Location','southwest');
finish(fig,out,'fig3_errors');
fig=figure('Position',[100 100 1050 550]);hold on;xi=linspace(.01,5,1500);
fill([.01 xi 5],[0 min(12,1./(.105125*xi)) 0],[1 .89 .85],'EdgeColor','none');
fill([.01 xi 5],[0 min(12,1./(.92*xi)) 0],[.83 .91 .97],'EdgeColor','none');
fill([0 1/2.3 1/2.3 0],[0 0 2.5 2.5],[.85 .85 .85],'EdgeColor','none');
h1=plot(xi,1./(.105125*xi),'Color',[.8 .2 .05],'LineWidth',1.3);
h2=plot(xi,1./(.92*xi),'Color',[0 .4 .7],'LineWidth',1.3);
h3=plot([1/2.3 1/2.3 0],[0 2.5 2.5],'Color',[.4 .4 .4],'LineWidth',1.1);
plot([1 1],[1 2],'ko','MarkerFaceColor','k');text(1.1,.65,'case 1');text(1.1,2,'case 2');
xlim([0 5]);ylim([0 12]);xlabel('\xi');ylabel('\eta');grid on;box on;
legend([h1 h2 h3],{'Separated spectral boundary','Merged spectral boundary','Unweighted row boundary'},'Location','northeast');
finish(fig,out,'fig4_regions');
% Derivative error: E_D(t) = max_j |z_j'(t) - z_ref,j'(t)|.
% Skip this plot if the input data contain only state trajectories.
if isfield(data,'derivative_errors')
 derivative_errors=data.derivative_errors;
 nh=size(derivative_errors,1);dc=lines(nh);
 labels=arrayfun(@(k)sprintf('history %d',k),1:nh,'UniformOutput',false);
 fig=figure('Position',[100 100 1050 420]);tiledlayout(1,2,'TileSpacing','compact');
 nexttile;hold on;
 for k=1:nh,plot(t,derivative_errors(k,:),'Color',dc(k,:),'LineWidth',1.1);end
 xlim([0 min(15,t(end))]);xlabel('t');ylabel('Derivative sup error');grid on;box on;
 nexttile;hold on;set(gca,'YScale','log');
 for k=1:nh,plot(t,max(derivative_errors(k,:),1e-16),'Color',dc(k,:),'LineWidth',1.1);end
 guide=max(derivative_errors(:,1))*exp(-p.lambda*t);
 plot(t,guide,'k--','LineWidth',1.1);
 xlim([0 min(25,t(end))]);ylim([1e-15 max(10,2*max(derivative_errors(:)))]);
 xlabel('t');ylabel('Derivative sup error');grid on;box on;set(gca,'YMinorGrid','off');
 legend([labels,{sprintf('slope guide e^{-%gt}',p.lambda)}],'Location','southwest');
 finish(fig,out,'fig5_derivative_errors');
else
 warning('BAM:MissingDerivativeErrors','Rerun run_bam to generate derivative-error data and figure 5.');
end
disp('BAM_PLOTS_COMPLETE');
end
function finish(fig,out,name)
exportgraphics(fig,fullfile(out,[name '.png']),'Resolution',300);
exportgraphics(fig,fullfile(out,[name '.pdf']),'ContentType','vector');
savefig(fig,fullfile(out,[name '.fig']));close(fig);
end
