% This script creates a new flow field
% We need to initialize U1, U2, U3, TH and P

% Since we are starting with a new flow, set TIME_STEP=0
% If you are continuing an existing simulation, don't do this (unless you want to reset the timestep)
TIME_STEP=0;

U1=zeros(NX,NY); % U1 is the velocity in the X direction
U2=zeros(NX,NY); % U2 is the velocity in the Y direction
U3=zeros(NX,NY); % U3 is the velocity in the Z direction
P=zeros(NX,NY); % P is the pressure
PHI=zeros(NX,NY); % PHI is a temporary variable used to make the velocity divergence free

TH=zeros(NX,NY,N_TH)+1;

% ********* User Input ********

% Start with a linear, unstable buoyancy profile
% for i=1:NX
% for j=1:NY
  % TH(i,j,:)=1-GYF(j);
% end
% end

% Example: For internal waves, create a stable buoyancy profile
% for i=1:NX
% for j=1:NY
    % b=1;
 % TH(i,j,:)= b / 2 * tanh((GXF(j) - LX / 8) / (LX / 50));
% And, optionally, add an initial perturbation
 % TH(i,j,:)=TH(i,j,:)+0.2*exp(-(GXF(i)-LX/2)^2/0.2^2-(GYF(j)-LY/2)^2/0.2^2);
% end
% end

% Create tanh buoyancy and velocity profiles
h=1/10; % Shear layer width
% NY=100; % The number of gridpoints
% dy=LY/NY; %The grid spacing - must be evenly spaced
% nu=1/5000; % Kinematic viscosity (or 1/Re)
% kappa=nu; % Diffusivity
S0=10; % Maximum shear
N0=sqrt(10); % Maximum buoyancy frequency

% Hyperbolic tangent velocity and buoyancy profiles
for i=1:NX
    for j=1:NY
        TH(i,j,:) = (N0^2/h) * tanh((GYF(j) - LY / 2) / (h));
        U1(i,j,:) = (S0/h) * tanh((GYF(j) - LY / 2) / (h));
    end
end

% Add a random perturbation to the velocity
U1=U1+0.001*(rand(NX,NY)-0.5);
U2=U2+0.001*(rand(NX,NY)-0.5);


% ********* End of User input *********

% Make sure that the new flow field satisfies the boundary conditions
rk_apply_bc_vel

