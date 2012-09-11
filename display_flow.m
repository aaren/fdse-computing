% This script is called every N_DISP_FLOW timesteps, and writes something to the screen

n_display=floor((TIME_STEP-1)/N_DISP_FLOW)+1; % Calculate the index for making a movie file
clf;

% ******** User Input *********

% Change the following lines to plot whatever you like during the simulation
surf(GXF(ii),GYF(jj),zeros(length(ii),length(jj))',TH(ii,jj,1)','EdgeColor','none'),view(0,90); shading interp; 
set(gca,'FontName','Times','FontSize',14);
xlabel('X');
ylabel('Y');
title('Buoyancy');
caxis([0 1]);
axis([GX(ii(1)) GX(ii(end)) GYF(jj(1)) GYF(jj(end))]);
colorbar
set(gca,'FontName','Times','FontSize',14);

% ******** End of User Input *********

% Here, we save the frame to a moviefile to replay after the simulation is done
mov(n_display)=getframe(gcf);
% And pause for a short amount of time to make sure that there is enough time for the screen to refresh
pause(0.01);

