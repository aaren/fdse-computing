% Set the upper bound on dt:
% Sometimes this is useful when starting with very small perturbations,
% or when a diffusion adds a limit to dt

% First, set the time-step to the maximum value
dt=DELTA_T;

% First, set a time-step constraint based on diffusion

dt=min(dt,0.5*min(min(min(DX(ii,jj),DY(ii,jj))))/(1/Re));
for n=1:N_TH
  dt=min(dt,dt*(1/Re)/(1/Re/PR(n))); % Take the larger of nu and kappa
end
% Make sure that we capture the inertial period (for rotating flows)
if (I_RO~=0)
  dt=min(dt,2*pi/I_RO/2);
end
% Make sure that we capture the buoyancy period (for stratified flows)
for n=1:N_TH
  if (RI(n)~=0)
    Nmax=sqrt(abs((GRAV_X(n)^2+GRAV_Y(n)^2+GRAV_Z(n)^2)*RI(n)*max(max((TH(ii,jp)-TH(ii,jj))./DY(ii,jj)))));
    dt=min(dt,0.05*2*pi/Nmax);
  end
end


% Now, enforce the CFL critiria
dt_x=CFL*min(min(DX(ii,jj)./abs(U1(ii,jj))));
dt_y=CFL*min(min(DY(3:NX-1,3:NY-1)./abs(U2(3:NX-1,3:NY-1))));  % Avoid j=2 which is a ghost cell

dt=min([dt dt_x dt_y]);

if (dt<=0)  error('DT<=0 in courant'); end;

% Set the size of the RK substeps (H_BAR)
H_BAR(1)=dt*(8/15);
H_BAR(2)=dt*(2/15);
H_BAR(3)=dt*(5/15);

