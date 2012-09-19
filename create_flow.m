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

% Create a tilted buoyancy profile to go with the tilted gravity unit vector
for i=1:NX
    for j=1:NY
        TH(i,j,:) = (1/sqrt(LX^2+LY^2)) * (GRAV_Y(1) * LY - GRAV_X(1) * LX) - (GRAV_Y(1) * GYF(j) - GRAV_X(1) * GXF(i));
        % And, optionally, add an initial perturbation
        % TH(i,j,:)=TH(i,j,:)+0.2*exp(-(GXF(i)-LX/2)^2/0.2^2-(GYF(j)-LY/2)^2/0.2^2);
    end
end


% Add a random perturbation to the velocity
U1=U1+0.001*(rand(NX,NY)-0.5);
U2=U2+0.001*(rand(NX,NY)-0.5);


% ********* End of User input *********

% Make sure that the new flow field satisfies the boundary conditions
rk_apply_bc_vel

