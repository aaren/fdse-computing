% This script analyzes the linear stability of a viscous stratified shear flow
% using the SSF.m solver written by William Smyth
% ********** User input parameters **********
LY=1; % The y-domain size
LX=10; % The x-domain size (sets the range of wavenumbers to search)
h=1/10; % Shear layer width
NY=100; % The number of gridpoints
dy=LY/NY; %The grid spacing - must be evenly spaced
nu=1/5000; % Kinematic viscosity (or 1/Re)
kappa=nu; % Diffusivity
S0=10; % Maximum shear
N0=sqrt(10); % Maximum buoyancy frequency
% ********** End user input parameters **********

% Create the y-cooridinate
y=(0:dy:1);

% Hyperbolic tangent velocity and buoyancy profiles
buoy=(N0^2/h)*tanh((y-LY/2)/(h));
vel=(S0/h)*tanh((y-LY/2)/(h));

% Create an x-wavenumber vector to explore solutions
kx=linspace(2*pi/LX,2*pi*20/LX,100);

clear growth
for k=1:length(kx)
  [sigma(:,k),lambda_w(:,:,k),lambda_b(:,:,k)]=SSF(y',vel',buoy',kx(k),0,nu,kappa,[0 0],[0 0],0);
end

plot(kx,real(sigma(1,:)),'b-');
set(gca,'FontName','Times','FontSize',14);
xlabel('k_x');
ylabel('Growth rate');
