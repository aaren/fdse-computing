% This script sets up and solves the linear stability for the bioconvection problem, outlined in Childress, Levandowsky, and Spiegel, JFM, 63, pp591-613, 1975
% Here, we assume that the vertical velocity, U0, and the diffusivity, kappa are both constant
% The problem is solved in non-dimensional variables, where
% h=kappa/U0 is the lengthscale
% h/U0 is the timescale

% User input variables
LY=5; % The y-domain size, H/h
LX=10; % The x-domain size (sets the range of wavenumbers to search)
NY=100; % The number of gridpoints
dy=LY/NY; %The grid spacing - must be evenly spaced
Pr=1; % The Prandtl number, nu/kappa
% The Rayleigh number, Ra=g*alpha*C0*kappa^2/nu/U0^3
% The Rayleigh number can either be a scalar or a vector
% If a scalar, the script will plot the growth rate vs. horiz. wavenumber
% If a vector, the script will plot the max. growth rate vs. Ra
Ra=17.44; % The Rayleigh number, Ra=g*alpha*C0*kappa^2/nu/U0^3
%Ra=linspace(0,1,10);

if (length(Ra)==1) % If we are only calculating one value of the Rayleigh No.

% Create an x-wavenumber vector to explore solutions
kx=linspace(2*pi/LX,2*pi*20/LX,50);

% Create the y-cooridinate
y=(0:dy:LY);

% Use the equilibrium solution for the base state, C0=1
C=exp(y);

clear growth
for k=1:length(kx)
  [sigma(k),lambda_w(:,k),lambda_c(:,k)]=BCVT(y',C',kx(k),Pr,Ra,[1 0],[1 0],0);
end

[sigma_max,kmax]=max(sigma);

figure
subplot(1,2,1)
plot(kx,real(sigma(:)),'b-');
set(gca,'FontName','Times','FontSize',14);
xlabel('k_x');
ylabel('Growth rate');
title('Maximum growth rate vs. horizontal wavenumber');
line([kx(1) kx(end)],[0 0]);
axis tight
subplot(1,2,2)
plot(lambda_w(:,kmax),y,'b-');
hold on
plot(lambda_c(:,kmax),y,'r-');
set(gca,'FontName','Times','FontSize',14);
xlabel('Eigenvector');
ylabel('Nondimensional fluid depth, y^*');
line([0 0],[y(1) y(end)]);
legend('Vertical velocity','Concentration');
title('Eigenvectors associated with most unstable mode');

else

figure
% Create an x-wavenumber vector to explore solutions
kx=linspace(2*pi/LX,2*pi*20/LX,50);

% Create the y-cooridinate
y=(0:dy:LY);
% Use the equilibrium solution for the base state, C0=1
C=exp(y);

clear growth
disp(['Looping through all ' num2str(length(Ra)) ' Rayleigh numbers']);
for i=1:length(Ra)
i
for k=1:length(kx)
  [sigma(:,k),lambda_w(:,:,k),lambda_c(:,:,k)]=BCVT(y',C',kx(k),Pr,Ra(i),[1 0],[1 0],0);
end

plot(Ra(i),max(sigma(1,:)),'.');
hold on
end
[sigma_max,kmax]=max(sigma);

set(gca,'FontName','Times','FontSize',14);
xlabel('Rayleigh Number');
ylabel('Maximum growth rate');
line([Ra(1) Ra(end)],[0 0]);

end
