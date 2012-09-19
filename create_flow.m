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

% set the buoyancy difference
b=1;
% set the lock width as fraction of lock length
w = 1/16;
% Initialise lock fluid
for i=1:NX
    for j=1:NY
    TH(i,j,1)= b / 2 * tanh((GXF(i) - LX * w) / (LX / 50));
    % And, optionally, add an initial perturbation
    % TH(i,j,:)=TH(i,j,:)+0.2*exp(-(GXF(i)-LX/2)^2/0.2^2-(GYF(j)-LY/2)^2/0.2^2);

    end
end


% Add linear perturbation to rest of fluid
% calculate lock width in integers
wnx = int16(NX*w);

% set dye
for i=wnx:NX
    for j=1:NY
        TH(i,j,2)=0;
    end
end

% apply stratification
for i=wnx:NX
    for j=1:NY
        TH(i,j,1) = TH(i,j,1) + b / 2 * (1 - GYF(j));
    end
end




% Add a random perturbation to the velocity
% U1=U1+0.001*(rand(NX,NY)-0.5);
% U2=U2+0.001*(rand(NX,NY)-0.5);


% ********* End of User input *********

% Make sure that the new flow field satisfies the boundary conditions
rk_apply_bc_vel

