% Here, add code to calculate and save any flow statistics that you might need later
% This script is called every N_SAVE_STATS timesteps

n_save=floor((TIME_STEP-1)/N_SAVE_STATS)+1; % Calculate the index for the saved array

TIME_save_stats(n_save)=TIME; % Save the current simulation time

% For example, calculate the mean and RMS of the velocity and scalar fields
U1_mean(n_save)=mean(mean(U1));
U1_rms(n_save)=sqrt(mean(mean((U1-U1_mean(n_save)).^2)));

U2_mean(n_save)=mean(mean(U2));
U2_rms(n_save)=sqrt(mean(mean((U2-U2_mean(n_save)).^2)));

if (SOLVE_U3) 
  U3_mean(n_save)=mean(mean(U3));
  U3_rms(n_save)=sqrt(mean(mean((U3-U3_mean(n_save)).^2)));
end

for n=1:N_TH
  TH_mean(n_save,n)=mean(mean(TH(:,:,n)));
  TH_rms(n_save,n)=sqrt(mean(mean((TH(:,:,n)-TH_mean(n_save,n)).^2)));
end

% When calculating derivatives, it is important to keep the finite difference
% stencil in mind, and be careful about where your gradients are defined
% For example, here we calculate y-derivatives using centered differences
% with the gradient defined on the GY grid (see the grid diagram in the docs folder)
% Shortcut arrays have been defined to do this without using loops:
%   ii=(2:NX-1), ip=(3:NX), im=(1:NX-2), etc.
% Calculate the x-average of the vertical derivative of U1
  dU1dy(:,n_save)=mean((U1(ii,jj)-U1(ii,jm))./DY(ii,jj),1);
% Same for each scalar
for n=1:N_TH
  dTHdy(:,n_save)=mean((TH(ii,jj)-TH(ii,jm))./DY(ii,jj),1);
end

% ********** Add User defined functions here ***********



% ********** End of User defined functions ************

