% Here, add code to save flow variables during the simulation for post-processing and visualization
n_save_flow=floor((TIME_STEP-1)/N_SAVE_FLOW)+1; % Calculate the index for making a movie file

U1_save(:,:,n_save_flow)=U1(:,:);
U2_save(:,:,n_save_flow)=U2(:,:);
if (SOLVE_U3) % If we are solving for U3, save it too
  U3_save(:,:,n_save_flow)=U3(:,:);
end
for n=1:N_TH % Save each scalar
  TH_save(:,:,n_save_flow,n)=TH(:,:,n);
end

% ****** Add User defined functions here *********


% ****** End of User defined functions *********

% Save the simulation time 
TIME_save(n_save_flow)=TIME;


