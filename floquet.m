% code for calculating the Floquet exponents of Stuart-Landau systems (reproducing Fig.3b in the paper)
function floquet(sigma,c1,c2)
% inputs for reproducing Fig.3b
% coupling strength -- sigma = 2.9;
% oscillator parameters -- c1 = -1.8; c2 = 4;

fonttype = 'Times';
fsize = 40;
txtattrib = {'FontName',fonttype,'FontSize',fsize,...
         'FontWeight','normal'};
txtattrib2 = {txtattrib{:},'Interpreter','Latex'};

tableau20 = [ 31, 119, 180; 174, 199, 232; 255, 127,  14; 255, 187, 120;...    
              44, 160,  44; 152, 223, 138; 214,  39,  40; 255, 152, 150;...    
             148, 103, 189; 197, 176, 213; 140,  86,  75; 196, 156, 148;...    
             227, 119, 194; 247, 182, 210;  23, 190, 207; 158, 218, 229;...    
             188, 189,  34; 219, 219, 141; 127, 127, 127; 199, 199, 199]/255;

n = 5; % number of temporal activities
m = 1000; % number of switching rates for each temporal activity
Amp = linspace(0,.2,n); % temporal activity
Omega = linspace(.15,40,m); % switching rate
Lambda = zeros(n,m); % Maximum Lyapunov exponent

% Calculate the Floquet exponents by numerically solving the ODE and constructing the principal fundamental matrix
for i = 1:n
	A = Amp(i);
	for j = 1:m
		w = Omega(j);
		time = [0 2*pi/w]; % one full period
		% constructing the principal fundamental matrix
		initial = [1,0];
		options = odeset('RelTol',1e-11,'AbsTol',1e-11); % the tolerance really matters here, not enough accuracy is the source of the artificial surge close to w=0 
		[ta,xa] = ode45(@(t,x)rhs(t,x,A,w,sigma,c1,c2),time,initial,options);
		initial = [0,1];
		[tb,xb] = ode45(@(t,x)rhs(t,x,A,w,sigma,c1,c2),time,initial,options);
		% extracting the maximum Lyapunov exponent
		B = [xa(end,:);xb(end,:)]';
		b = trace(B);
		c = exp(-2*(1+sigma)*2*pi/w);
		rho = (b+sqrt(b^2-4*c))/2;
		Lambda(i,j) = real(log(rho))/(2*pi/w);
		%Lambda(i,j) = max(log(abs(eig(B)))/(2*pi/w));
	end
end

% Calculate the averaged master stability function for comparison at the slow-switching limit
for i = 1:n
	A = Amp(i);
	alpha = linspace(sigma*(1-2*A),sigma*(1+2*A),m);
	MSF = -alpha - 1 + sqrt(1-2*c1*c2*alpha-c1^2*alpha.^2);
	idx = linspace(-1,1,m+1);
	dif = asin(idx(2:end))-asin(idx(1:end-1));
	w = (dif)/norm(dif,1);
	MLE(i) = sum(MSF.*w);
end
MLE

% Plotting Fig.3b
figure
set(0,'DefaultAxesFontSize',25)
hold on
for i = 1:n
	plot(Omega,Lambda(i,:),'color',tableau20(n+1-i,:),'LineWidth',4)
end
for i = 2:n
	plot(Omega,MLE(i)*ones(1,m),'--','color',tableau20(n+1-i,:),'LineWidth',2)
end
hold off
box on
xlim([0,40])
ylim([-.31,.1])
yticks([-.3 -.2 -.1 0 .1])
yticklabels({'-0.3','-0.2','-0.1','0.0','0.1'})
xlabel('$\omega$',txtattrib2{:})
ylabel('$\Gamma$',txtattrib2{:})

legend({'$A=0.00$','$A=0.05$','$A=0.10$','$A=0.15$','$A=0.20$'},'Location','southeast','NumColumns',2,txtattrib2{:},'FontSize',20)
legend('boxoff')

set(gcf, 'PaperPosition', [0 0 8 5]);
set(gcf, 'PaperSize', [8 5]);
saveas(gcf,'Lambda.pdf');

end