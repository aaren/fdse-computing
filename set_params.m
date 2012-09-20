% This file contains definitions for various physical and runtime paramters.

% ********* User Input **********

% Physical paramters:
Re=5000.0;			% Reynolds number
LX=1;				% x - Domain size
LY=1;				% y - Domain size

% In some applications, we don't want to allow velocity in/out of the 2D plane
% Specify whether to timestep U3 using this flag: (0=FALSE, 1=TRUE)
SOLVE_U3=0;

% Number of gridpoints 
% Note - Since there are 2 ghost cells in each direction, there will be (NX-2)*(NY-2) computational gridpoints
NX=100;
NY=100;

% The number of active/passive scalars
N_TH=1;			

% Set the following parameters for each of the scalars

% The buoyancy term, added to the RHS of the momentum equation is RI(n)*(GRAV_X,GRAV_Y,GRAV_Z)*TH(n)
% (GRAV_X,GRAV_Y,GRAV_Z) is intended to be a unit vector, giving the direction of the gravitational acceleration
% and RI(n) controls the amplitude of the buoyant force.
% Using this convention RI(n) should be positive if TH(n) is density and negative if TH(n) is buoyancy

for n=1:N_TH
  if (n==1)
    PR(n)=1.0;		% Prandtl (or Schmidt) number
    RI(n)=-1.0;		% Richardson number (normalized g/rho0)
  end
end 

% Specify angle of rotation for the domain (to look at slopes)
alpha = 30;
% Specify the direction of the gravitational unit vector 
% (Do for each scalar, but these are probably the same!) 
% create a vector, which will be normalised.
g = 1
g_x = -g * sind(alpha);
g_y = -g * cosd(alpha) ;
g_z = 0;
for n=1:N_TH
  GRAV_X(n) = g_x;
  GRAV_Y(n) = g_y;
  GRAV_Z(n) = g_z;
end

% Optionally, specify background horizontal gradients for each scalar
% specify total magnitude
drho = 1;
for n=1:N_TH
  DTHDX(n)= -drho * sin(alpha);
  DTHDY(n)= drho * cos(alpha);
  DTHDZ(n)=0;
end

% External body force  (can be scalars or NX,NY matrices)
% For more complicated forces (involving time, or derivatives),
%  use the script user_rhs.m instead of the body force terms
FORCE_X=0;
FORCE_Y=0;
FORCE_Z=0;
FORCE_TH=0;

% Parameters for rotating flows - for no rotation, set I_RO=0;
% The magnitude of the Coriolis parameter is set using 1/Rossby number
% The direction of the rotation vector is set by a unit vector CORI
% The arrays CORI_X, etc. can be spatially varying to account for the Beta effect.
I_RO=0;                     % 1/Rossby number
CORI_X=zeros(NX,NY);		% x-component of the Coriolis unit vector
CORI_Y=zeros(NX,NY);	   	% x-component of the Coriolis unit vector
CORI_Z=zeros(NX,NY);		% x-component of the Coriolis unit vector

% Runtime parameters:
RESTART=0; 		  	% A flag determining whether to run create_flow.m
					% If RESTART=0, run create_flow.m to set up a new flow field
					% If RESTART=1, restart an old simulation the variables in memory
% Timestepping parameters:
N_TIME_STEPS=1000;	% Number of timesteps
VARIABLE_DT=1;		% 1=Use CFL number to set variable timestep
					% 0=Use DELTA_T for constant timestep
DELTA_T=0.1;		% When using a fixed timestep (VARIABLE_DT=0), this is the timestep size
                    % When using a variable timestep (VARIABLE_DT=1), this is the maximum timestep
CFL=1.0;			% When using a variable timestep (VARIABLE_DT=1), the timestep is calculated from this CFL number

% Saving and display intervals
N_DISP_FLOW=50;		% Specify how often (in timesteps) to call display_flow
N_SAVE_FLOW=1;	% Specify how often (in timesteps) to save flow fields
N_SAVE_STATS=5;		% Specify how often (in timesteps) to calculate and save flow statistics

% MATLAB parameters for solving the implicit linear systems
iter_tol=1e-6;		% The error tolerance for the iterative solvers
iter_max=100;		% The maximum number of iterations  

% ******** End of User input **********

% The following are the sizes of each Runge-Kutta substep - don't change these values!
H_BAR(1)=DELTA_T*(8/15);
H_BAR(2)=DELTA_T*(2/15);
H_BAR(3)=DELTA_T*(5/15);
BETA_BAR(1)=1.0;
BETA_BAR(2)=25.0/8.0;
BETA_BAR(3)=9.0/4.0;
ZETA_BAR(1)=0.0;
ZETA_BAR(2)=-17.0/8.0;
ZETA_BAR(3)=-5.0/4.0;



