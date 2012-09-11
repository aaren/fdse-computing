% Diablo -> DNS In A Box, Laptop Optimized                MATLAB/Octave version 
%
% This code was written by John R. Taylor, June 2012, Cambridge, UK
% and is based on an algorithm developed by Thomas Bewley at UC San Diego
%
% This code computes incompressible flow in a box.
% Primative variables (U1,U2,U3,P) are used, and continuity is enforced with a
% fractional step algorithm.
%
% SPATIAL DERIVATIVES
% Spatial derivatives are calculated using second order finite-differences
% on a staggered grid with a momentum and energy-conserving scheme.
%
% TIME ADVANCEMENT
% Time advancement is accomplished with a low-storage 3rd order accurate Runge-Kutta-Wray method
% All viscous terms are treated with semi-implicit Crank-Nicolson
% All other terms are treated explicitly
%
% This code is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version. This code is distributed in the hope that it
% will be useful, but WITHOUT ANY WARRANTY; without even the implied
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details. You should have received a
% copy of the GNU General Public License along with this code; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place - Suite
% 330, Boston, MA 02111-1307, USA.
%

disp('****** WELCOME TO DIABLO ******');

%***************************************************************
% Here, we go through the usual steps to set up a new simulation
%***************************************************************

% Check to make sure that the code directory is in our current path
if (~exist('rk_step.m'))
  path(path,'./code'); % Add the code directory to our path
end

set_params;  % Define various physical and computational parameters

define_shortcuts;  % Define a few useful arrays for shorthand notation

create_grid; % Define a computational grid

set_bcs;  % Prescribe the boundary conditions

if (~RESTART) % If we aren't restarting from an old simulation...
TIME=0; create_flow;  % Set the TIME variable to zero, and specify the initial flow field
end

precondition; % Create preconditioner matrices, and perform an LU decomposition for the pressure solve

courant; % Make sure the initial timestep obeys the CFL criteria

%********************************************************
% Now, timestep the equations, repeatedly calling rk_step
%********************************************************

disp(['Ready to run ' int2str(N_TIME_STEPS) ' timesteps ']);

tic;
for TIME_STEP=TIME_STEP+1:N_TIME_STEPS  % This is the main time-stepping loop

% Every N_DISP_FLOW timesteps, display something...
  if (mod(TIME_STEP-1,N_DISP_FLOW)==0) display_flow;  end 
% Every N_SAVE_FLOW timesteps, save the flow field...
  if (mod(TIME_STEP-1,N_SAVE_FLOW)==0) save_flow;  end 
% Every N_SAVE_FLOW timesteps, calculate and save statistics...
  if (mod(TIME_STEP-1,N_SAVE_STATS)==0) save_stats;  end 
% Display a message to show our progress...
  if (mod(TIME_STEP,10)==0) disp(['Now beginning TIME_STEP = ' num2str(TIME_STEP) ', TIME = ' num2str(TIME)]); end

  for RK_STEP=1:3  % Perform 3 Runge-Kutta substeps, calling rk_step each time
    rk_step;
  end
   
% We're done with a full timestep, update the simulation time
  if (VARIABLE_DT) 
    TIME=TIME+dt;
  else
    TIME=TIME+DELTA_T; 
  end

% If we are using a variable timestep, update it now
  if (VARIABLE_DT) courant; end  

end
toc;

% Replay the movie frames saved in display_flow.m
movie(mov)

disp('****** Hello world!  Have a nice day! ******');  % We're done!

