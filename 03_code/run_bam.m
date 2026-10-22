function run_bam(output_tag)
% Uses base MATLAB functions. Archived run: MATLAB R2026b, dde23.
if nargin<1,output_tag='matlab';end
tic;disp('BAM full MATLAB simulation started');
here=fileparts(mfilename('fullpath')); p=jsondecode(fileread(fullfile(here,'parameters.json')));
% New runs never overwrite the archived paper data.
out=fullfile(here,'..','outputs',output_tag); if ~exist(out,'dir'),mkdir(out);end
p.c=p.c(:);p.d=p.d(:);p.input_sin=p.input_sin(:);
p.input_cos=p.input_cos(:);p.input_ergodic=p.input_ergodic(:);
rng(p.seed); set(groot,'defaultFigureVisible','off');
opts=ddeset('RelTol',p.relative_tolerance,'AbsTol',p.absolute_tolerance,'MaxStep',p.max_step);
lags=[p.delay_linear,p.delay_product_1,p.delay_product_2];
rhs=@(t,z,Z) equation(t,z,Z,p);
ref=dde23(rhs,lags,zeros(8,1),[-p.burn_in,p.end_time],opts);
disp('Reference completed');
fineopts=ddeset(opts,'RelTol',p.relative_tolerance/10,'AbsTol',p.absolute_tolerance/10,'MaxStep',p.max_step/2);
refFine=dde23(rhs,lags,zeros(8,1),[-p.burn_in,p.end_time],fineopts);
disp('Refined reference completed');
refLong=dde23(rhs,lags,zeros(8,1),[-p.burn_in_check,p.end_time],fineopts);
disp('Long reference completed');
t=0:p.output_step:p.end_time;
zr=deval(ref,t); zf=deval(refFine,t); zl=deval(refLong,t);
% Derivatives are evaluated from the DDE right-hand side, not finite differences.
% At an initial-history junction this gives the right-hand derivative.
reference_derivatives=state_derivatives(ref,t,zeros(8,1),rhs,lags);
nh=size(p.histories,1); trajectories=zeros(4,numel(t),nh); errors=zeros(nh,numel(t));
trajectory_derivatives=zeros(4,numel(t),nh);
derivative_errors=zeros(nh,numel(t));
for k=1:nh
 hist=[p.histories(k,:)';zeros(4,1)];
 sol=dde23(rhs,lags,hist,[0,p.end_time],opts); z=deval(sol,t);
 trajectories(:,:,k)=z(1:4,:); errors(k,:)=max(abs(z(1:4,:)-zr(1:4,:)),[],1);
 trajectory_derivatives(:,:,k)=state_derivatives(sol,t,hist,rhs,lags);
 derivative_errors(k,:)=max(abs(trajectory_derivatives(:,:,k)-reference_derivatives),[],1);
 if k==1
  sf=dde23(rhs,lags,hist,[0,p.end_time],fineopts); zz=deval(sf,t);
  validation.trajectory_refinement=max(max(abs(zz(1:4,:)-z(1:4,:))));
  zzp=state_derivatives(sf,t,hist,rhs,lags);
  validation.trajectory_derivative_refinement=max(max(abs(zzp-trajectory_derivatives(:,:,k))));
 end
end
validation.reference_refinement=max(max(abs(zr(1:4,:)-zf(1:4,:))));
validation.burn_in_difference=max(max(abs(zf(1:4,:)-zl(1:4,:))));
zfp=state_derivatives(refFine,t,zeros(8,1),rhs,lags);
zlp=state_derivatives(refLong,t,zeros(8,1),rhs,lags);
validation.reference_derivative_refinement=max(max(abs(reference_derivatives-zfp)));
validation.derivative_burn_in_difference=max(max(abs(zfp-zlp)));
joint_errors=max(errors,derivative_errors);
validation.matlab_version=version; validation.solver='dde23';
validation.elapsed_seconds=toc;
validation.initial_errors=errors(:,1); validation.error_t10=errors(:,find(t>=10,1));
validation.error_t20=errors(:,find(t>=20,1));
validation.derivative_initial_errors=derivative_errors(:,1);
validation.derivative_error_t10=derivative_errors(:,find(t>=10,1));
validation.derivative_error_t20=derivative_errors(:,find(t>=20,1));
fid=fopen(fullfile(out,'matlab_validation.json'),'w'); fprintf(fid,'%s',jsonencode(validation)); fclose(fid);
writematrix([t',zr(1:4,:)',errors'],fullfile(out,'trajectories.csv'));
save(fullfile(out,'simulation.mat'),'p','t','zr','trajectories','errors','validation', ...
 'reference_derivatives','trajectory_derivatives','derivative_errors','joint_errors');
writematrix([t',reference_derivatives',derivative_errors'],fullfile(out,'derivative_errors.csv'));
names={'x_1','x_2','y_1','y_2'};
fig=figure('Position',[100 100 1100 700]); tiledlayout(2,2,'TileSpacing','compact');
for j=1:4,nexttile;plot(t,zr(j,:),'LineWidth',1.1);xlabel('t');ylabel(names{j});grid on;end
exportgraphics(fig,fullfile(out,'fig1_states.png'),'Resolution',220);close(fig);
fig=figure('Position',[100 100 1100 700]); tiledlayout(2,2,'TileSpacing','compact');
for j=1:4
 nexttile;plot(t,zr(j,:),'k--','LineWidth',1.3);hold on;
 for k=1:nh,plot(t,trajectories(j,:,k),'LineWidth',0.9);end
 xlim([0 15]);xlabel('t');ylabel(names{j});grid on;
 if j==1,legend('reference','history 1','history 2','history 3','Location','best');end
end
exportgraphics(fig,fullfile(out,'fig2_histories.png'),'Resolution',220);close(fig);
fig=figure('Position',[100 100 1050 400]);tiledlayout(1,2,'TileSpacing','compact');
nexttile;plot(t,errors,'LineWidth',1.1);xlim([0 15]);xlabel('t');ylabel('State sup error');grid on;
nexttile;semilogy(t,max(errors,1e-16),'LineWidth',1.1);hold on;
semilogy(t,max(errors(:,1))*exp(-p.lambda*t),'k--','LineWidth',1.1);
xlim([0 25]);ylim([1e-12 10]);xlabel('t');ylabel('State sup error');grid on;
legend('history 1','history 2','history 3','slope guide e^{-0.2t}','Location','southwest');
exportgraphics(fig,fullfile(out,'fig3_errors.png'),'Resolution',220);close(fig);
fig=figure('Position',[100 100 1050 550]);hold on;
xi=linspace(0.01,5,1500);
fill([0.01 xi 5],[0 min(12,1./(.105125*xi)) 0],[1 .89 .85],'EdgeColor','none');
fill([0.01 xi 5],[0 min(12,1./(.92*xi)) 0],[.83 .91 .97],'EdgeColor','none');
fill([0 1/2.3 1/2.3 0],[0 0 2.5 2.5],[.85 .85 .85],'EdgeColor','none');
h1=plot(xi,1./(.105125*xi),'Color',[.8 .2 .05],'LineWidth',1.3);
h2=plot(xi,1./(.92*xi),'Color',[0 .4 .7],'LineWidth',1.3);
h3=plot([1/2.3 1/2.3 0],[0 2.5 2.5],'Color',[.4 .4 .4],'LineWidth',1.1);
plot([1 1],[1 2],'ko','MarkerFaceColor','k');text(1.1,1,'case 1');text(1.1,2,'case 2');
xlim([0 5]);ylim([0 12]);xlabel('\xi');ylabel('\eta');grid on;box on;
legend([h1 h2 h3],{'Separated spectral boundary','Merged spectral boundary','Unweighted row boundary'},'Location','northeast');
exportgraphics(fig,fullfile(out,'fig4_regions.png'),'Resolution',220);close(fig);
for eta=[1 2]
 A=[zeros(2),p.P+p.c*ones(1,2);eta*(p.Q+p.d*ones(1,2)),zeros(2)];
 B=[zeros(2),p.B;eta*p.D,zeros(2)];
 J=A/p.leakage+2*B;G=2*(A+B);
 Al=[zeros(2),p.P*exp(p.lambda*p.delay_linear)+p.c*[exp(p.lambda*p.delay_product_1),exp(p.lambda*p.delay_product_2)];...
 eta*(p.Q*exp(p.lambda*p.delay_linear)+p.d*[exp(p.lambda*p.delay_product_1),exp(p.lambda*p.delay_product_2)]),zeros(2)];
 Bl=B*p.kernel_beta/(p.kernel_beta-p.lambda);
 R=eye(4)/(p.leakage-p.lambda);S=eye(4)+p.leakage*R;
 H=[R*Al,R*Bl;S*Al,S*Bl];Jl=Al*R+Bl*S;
 cert(eta).eta=eta;cert(eta).row_norm=norm(G,inf);
 cert(eta).rho_G=max(abs(eig(G)));cert(eta).rho_J=max(abs(eig(J)));
 cert(eta).rho_J_lambda=max(abs(eig(Jl)));
 cert(eta).M_lambda=norm((eye(8)-H)\[2*ones(4,1);(1+p.leakage)*ones(4,1)],inf);
 cert(eta).positive_vector=(eye(4)-J)\ones(4,1);
 cert(eta).residual=cert(eta).positive_vector-J*cert(eta).positive_vector;
end
fid=fopen(fullfile(out,'condition_matlab.json'),'w');fprintf(fid,'%s',jsonencode(cert));fclose(fid);
disp(validation)
plot_bam_results(out);
disp('BAM_MATLAB_COMPLETE');
end

function dz=equation(t,z,Z,p)
 x=z(1:2);y=z(3:4);r=z(5:8);
 u=p.input_sin.*sin(t)+p.input_cos.*cos(sqrt(2)*t)+p.input_ergodic/(1+t*t);
 dx=-p.leakage*x+(p.coefficient_base+p.coefficient_amplitude*sin(t))*p.P*tanh(Z(3:4,1)) ...
    +p.B*r(3:4)+p.c*tanh(Z(3,2))*tanh(Z(4,3))+u(1:2);
 dy=-p.leakage*y+p.scale_y*((p.coefficient_base+p.coefficient_amplitude*cos(sqrt(2)*t))*p.Q*tanh(Z(1:2,1)) ...
    +p.D*r(1:2)+p.d*tanh(Z(1,2))*tanh(Z(2,3)))+u(3:4);
 dz=[dx;dy;p.kernel_beta*(tanh([dx;dy])-r)];
end

function zp=state_derivatives(sol,t,history,rhs,lags)
% Reconstruct all delayed states from this solution and its constant history.
% Only the first four components are x_1', x_2', y_1', y_2'.
z=deval(sol,t); zp=zeros(4,numel(t)); t0=sol.x(1);
for q=1:numel(t)
 Z=zeros(8,numel(lags));
 for j=1:numel(lags)
  td=t(q)-lags(j);
  if td<=t0
   Z(:,j)=history(:);
  else
   Z(:,j)=deval(sol,td);
  end
 end
 dz=rhs(t(q),z(:,q),Z);
 zp(:,q)=dz(1:4);
end
end
