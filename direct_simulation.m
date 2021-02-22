% code for simulating optoelectronic oscillators on temporal networks (reproducing Fig.6 in the paper)
function direct_simulation(Amp,b,period)
% Amp -- Temporal activity
% b -- standard deviation of the fluctuation (before normalization)
% period -- switching period of the base temporal network

n = 99; % network size
m = 50;
t = 1e5; % total iteration steps

% coupling strength
sigma = 1.05;
% oscillator parameter
beta = 2.8;

% initial condition
x(1:n,1) = ones(n,1) + 1e-3*rand(n,1);

% noise
eta = 1e-16;

% graph Laplacian for t=1
adj = (ones(n) - eye(n))/n;
rowsum = sum(adj,2);
lap = diag(rowsum) - adj;

for ih = 2:t
  % evolve the system
  for i = 1:n
  	x(i,ih) = f(i,ih);
  end
  % adjacency matrix for the discrete-switching network
  if mod(floor(ih/period),2)==0
    adj = (ones(n)-eye(n) + [(6-8/(n+1))*Amp*(ones(m)-eye(m)),-2*Amp*ones(m,n-m);-2*Amp*ones(n-m,m),-2*Amp*(ones(n-m)-eye(n-m))])/n;
  else
    adj = (ones(n)-eye(n) - [(6-8/(n+1))*Amp*(ones(m)-eye(m)),-2*Amp*ones(m,n-m);-2*Amp*ones(n-m,m),-2*Amp*(ones(n-m)-eye(n-m))])/n;
  end
  % add random fluctuations to the network
  fluc = b*randn(n)/n;
  fluc = fluc - diag(diag(fluc));
  adj = adj + fluc;
  % generate graph Laplacian for the next time step
  rowsum = sum(adj,2);
  lap = diag(rowsum) - adj;
end

% dynamical equation
function z = f(i,ih)
	z = beta*I(x(i,ih-1)) - sigma*lap(i,:)*I(x(:,ih-1)) + normrnd(0,eta);
end

% node dynamics
function y = I(x)
	y = sin(x+pi/4).^2;
end

% sync error
syncerr = std(x);


% Spacetime plot
set(0,'DefaultAxesFontSize',40)
fonttype = 'Times';
fsize = 50;
txtattrib = {'FontName',fonttype,'FontSize',fsize,...
         'FontWeight','normal'};
txtattrib2 = {txtattrib{:},'Interpreter','Latex'};

figure
hAxis(1)=subplot(1,1,1);
pos = get( hAxis(1), 'Position');
pos(1)=.1;
pos(2)=.25;
pos(3)=.8;
pos(4)=.72;
set(hAxis(1), 'Position', pos);
set(gcf, 'PaperPosition', [0 0 20 8]);
set(gcf, 'PaperSize', [20 8]);
h = imagesc(1:t,1:n,x);
set(gca, 'YDir', 'normal');
colormap(flipud([lbmap(201,'brownblue')]));
caxis([0 3])
cbh = colorbar('Location','eastoutside');
set(cbh,'YTick',[0,1,2,3])

xlim([0,t])
ylim([0,n]+.5)
yticks([1,50,99])
xlabel('$t$', txtattrib2{:})
ylabel('$i$', txtattrib2{:})
saveas(gcf,'space_time.pdf')

end
