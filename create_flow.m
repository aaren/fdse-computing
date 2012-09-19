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

% set the buoyancy difference for the stratification
b=4;

% set the lock width as fraction of lock length
w = 1/16;
% calculate lock width in integers
wnx = int16(NX*w);

% Initialise lock fluid

% set buoyancy of lock fluid compared to stratifcation difference (lb=1) 
lb=0.5;

for i=1:NX
    for j=1:NY
    TH(i,j,1)= lb*b / 2 * (tanh(-(GXF(i) - LX * w) / (LX / 50))+1);
    % And, optionally, add an initial perturbation
    % TH(i,j,:)=TH(i,j,:)+0.2*exp(-(GXF(i)-LX/2)^2/0.2^2-(GYF(j)-LY/2)^2/0.2^2);

    end
end

% apply stratification
for i=wnx:NX
    for j=1:NY
        TH(i,j,1) = b  * (GYF(j));
    end
end

% set dye
for i=wnx:NX
    for j=1:NY
        TH(i,j,2)=0;
    end
end


% for non tanh lock fluid:
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % set the lock width as fraction of lock length
% lw = 1/16;
% lh=1;
% % calculate lock width and height in integers
% wnx = int16(NX*lw);
% hny = int16(NY*lh);
% 
% % set dye
% for i=wnx:NX
%     for j=hny:NY
%         TH(i,j,2)=0;
%     end
% end

% % Initialise lock fluid
% TH(1:wnx,1:hny,1)=b / 2 ;
% % 
% 




% Add a random perturbation to the velocity
% U1=U1+0.001*(rand(NX,NY)-0.5);
% U2=U2+0.001*(rand(NX,NY)-0.5);


% ********* End of User input *********

% Make sure that the new flow field satisfies the boundary conditions
rk_apply_bc_vel

