% This script performs one Runge-Kutta substep

% Set a few useful constants
% H_BAR, BETA_BAR, and ZETA_BAR are constants associated with the Runge-Kutta algorithm
TEMP1=(1/Re)*H_BAR(RK_STEP)/2;
TEMP2=H_BAR(RK_STEP)/2;
TEMP3=ZETA_BAR(RK_STEP)*H_BAR(RK_STEP);
TEMP4=H_BAR(RK_STEP);
TEMP5=BETA_BAR(RK_STEP)*H_BAR(RK_STEP);

% Store the old velocity field in the RHS array
R1=U1;
R2=U2;
R3=U3;
% Store the old scalar fields in the RHS array
if (N_TH~=0)
  RTH(ii,jj,1:N_TH)=TH(ii,jj,1:N_TH);
end

% As long as this isn't the first R-K step, add the previous RHS array
if (RK_STEP>1)
  R1=R1+TEMP3*F1;
  R2=R2+TEMP3*F2;
  R3=R3+TEMP3*F3;
end
% Do the same thing for the scalar arrays
if ((RK_STEP>1)&&(N_TH~=0))
  RTH(ii,jj,1:N_TH)=RTH(ii,jj,1:N_TH)+TEMP3*FTH(ii,jj,1:N_TH);
end

% Add the pressure gradient to the RHS arrays
R1(ii,jj)=R1(ii,jj)-TEMP4*(P(ii,jj)-P(im,jj))./DX(ii,jj);
R2(ii,jj)=R2(ii,jj)-TEMP4*(P(ii,jj)-P(ii,jm))./DY(ii,jj);

% Now start to build up the explicit RHS vectors
F1=zeros(size(U1));
F2=zeros(size(U2));
F3=zeros(size(U3));
FTH=zeros(size(TH));

% Optionally, add user-defined terms to the RHS forcing arrays (Fi,FTH)
user_rhs;

% Add the body force terms
F1=F1+FORCE_X;
F2=F2+FORCE_Y;
F3=F3+FORCE_Z;
FTH=FTH+FORCE_TH;

% Add the Coriolis force, interpolate onto the appropriate grid
if (I_RO~=0)
  F1(ii,jj)=F1(ii,jj)-I_RO*0.5*(CORI_Y(ii,jj).*U3(ii,jj)-CORI_Z(ii,jj).*(U2(ii,jp)+U2(ii,jj))/2 ...
                               +CORI_Y(im,jj).*U3(im,jj)-CORI_Z(im,jj).*(U2(im,jp)+U2(im,jj))/2);

  F2(ii,jj)=F2(ii,jj)-I_RO*0.5*(-CORI_X(ii,jj).*U3(ii,jj)+CORI_Z(ii,jj).*(U1(ip,jj)+U1(ii,jj))/2 ...
                                -CORI_X(ii,jm).*U3(ii,jm)+CORI_Z(ii,jm).*(U1(ip,jm)+U1(ii,jm))/2);

  F3(ii,jj)=F3(ii,jj)-I_RO*(CORI_X(ii,jj).*(U2(ii,jp)+U2(ii,jj))/2 ...
                           -CORI_Y(ii,jj).*(U1(ip,jj)+U1(ii,jj))/2);
% This is the interpolated version of:  
%  F1=F1-I_RO*(CORI_Y.*U3-CORI_Z.*U2);
%  F2=F2-I_RO*(CORI_Z.*U1-CORI_X.*U3);
%  F3=F3-I_RO*(CORI_X.*U2-CORI_Y.*U1);
end

% Add the buoyancy force for each scalar
% Use second order accurate interpolation
for n=1:N_TH
  F1(ii,jj)=F1(ii,jj)+RI(n)*GRAV_X(n)*(TH(ii,jj,n).*DXF(im,jj)+TH(im,jj,n).*DXF(ii,jj))./(2*DX(ii,jj));
  F2(ii,jj)=F2(ii,jj)+RI(n)*GRAV_Y(n)*(TH(ii,jj,n).*DYF(ii,jm)+TH(ii,jm,n).*DYF(ii,jj))./(2*DY(ii,jj));
  F3(ii,jj)=F3(ii,jj)+RI(n)*GRAV_Z(n)*TH(ii,jj,n);
end

% Add the nonlinear terms (explicit time-stepping)
% d/dx(U1*U1) (overbar(U1)*overbar(U1))
F1(ii,jj)=F1(ii,jj) - ((0.5*(U1(ip,jj)+U1(ii,jj))).^2 - (0.5*(U1(im,jj)+U1(ii,jj))).^2) ...
                                                 ./DX(ii,jj);
% d/dy(U1*U2) (breve(U2)*overbar(U1))
F1(ii,jj)=F1(ii,jj) - ( (U1(ii,jj)+U1(ii,jp))/2.*(DXF(ii,jp).*U2(ii,jp)+DXF(im,jp).*U2(im,jp))/2./DXF(ii,jp) ...
                       -(U1(ii,jj)+U1(ii,jm))/2.*(DXF(ii,jj).*U2(ii,jj)+DXF(im,jj).*U2(im,jj))/2./DXF(ii,jj) ) ...
                                                       ./DYF(ii,jj); 
% y-momentum equation
% d/dy(U2*U2) - (overbar(U2)*overbar(U2))
F2(ii,jj)=F2(ii,jj) - ((0.5*(U2(ii,jj)+U2(ii,jp))).^2 - (0.5*(U2(ii,jj)+U2(ii,jm))).^2)...
                                    ./DY(ii,jj);
% d/dx(U1*U2) - (breve(U1)*overbar(U2))
F2(ii,jj)=F2(ii,jj) - ( (DYF(ip,jj).*U1(ip,jj)+DYF(ip,jm).*U1(ip,jm))/2./DY(ip,jj).*(U2(ip,jj)+U2(ii,jj))/2 ...
                       -(DYF(ii,jj).*U1(ii,jj)+DYF(ii,jm).*U1(ii,jm))/2./DY(ii,jj).*(U2(ii,jj)+U2(im,jj))/2 ) ...
                                                       ./DXF(ii,jj);

% z-momentum equation
% d/dx(U1*U3) (breve(U1)*overbar(U3))
F3(ii,jj)=F3(ii,jj) - (U1(ip,jj).*(U3(ip,jj)+U3(ii,jj))/2 - U1(ii,jj).*(U3(ii,jj)+U3(im,jj))/2) ...
                                                     ./DXF(ii,jj);
% d/dy(U2*U3) (brebe(U2)*overbar(U3))
F3(ii,jj)=F3(ii,jj) - (U2(ii,jp).*(U3(ii,jp)+U3(ii,jj))/2 - U2(ii,jj).*(U3(ii,jj)+U3(ii,jm))/2) ...
                                                     ./DYF(ii,jj); 

for n=1:N_TH % loop over all scalars

% Add a term owing to a background scalar gradient
FTH(ii,jj,n)=FTH(ii,jj,n)-0.5*(U1(ii,jj)+U1(ip,jj))*DTHDX(n);
FTH(ii,jj,n)=FTH(ii,jj,n)-0.5*(U2(ii,jj)+U2(ii,jp))*DTHDY(n);
FTH(ii,jj,n)=FTH(ii,jj,n)-U3(ii,jj)*DTHDZ(n);

% Add the scalar advection terms (explicit timestepping)
% d/dx(U1*TH) (breve(U1)*overbar(TH))
FTH(ii,jj,n)=FTH(ii,jj,n) - (U1(ip,jj).*(TH(ip,jj,n)+TH(ii,jj,n))/2 - U1(ii,jj).*(TH(ii,jj,n)+TH(im,jj,n))/2) ...
                                                     ./DXF(ii,jj);
% d/dy(U2*TH) (breve(U2)*overbar(TH))
FTH(ii,jj,n)=FTH(ii,jj,n) - (U2(ii,jp).*(TH(ii,jp,n)+TH(ii,jj,n))/2 - U2(ii,jj).*(TH(ii,jj,n)+TH(ii,jm,n))/2) ...
                                                     ./DYF(ii,jj);                                             
end % End loop over scalars           
           
% We are done computing nonlinear advection terms
% Add Fi to Ri
R1(ii,jj)=R1(ii,jj)+TEMP5*F1(ii,jj);
R2(ii,jj)=R2(ii,jj)+TEMP5*F2(ii,jj);
R3(ii,jj)=R3(ii,jj)+TEMP5*F3(ii,jj);

% Add the explicit part of the Crank-Nicolson viscous terms
R1(ii,jj)=R1(ii,jj)+TEMP1 * ...
          ( ((U1(ip,jj)-U1(ii,jj))./DXF(ii,jj) - (U1(ii,jj)-U1(im,jj))./DXF(im,jj)) ...
                                           ./DX(ii,jj) ...
           +((U1(ii,jp)-U1(ii,jj))./DY(ii,jp) - (U1(ii,jj)-U1(ii,jm))./DY(ii,jj)) ...
                                           ./DYF(ii,jj) );
R2(ii,jj)=R2(ii,jj)+TEMP1 * ...
          ( ((U2(ip,jj)-U2(ii,jj))./DX(ip,jj) - (U2(ii,jj)-U2(im,jj))./DX(ii,jj)) ...
                                           ./DXF(ii,jj) ...
           +((U2(ii,jp)-U2(ii,jj))./DYF(ii,jj) - (U2(ii,jj)-U2(ii,jm))./DYF(ii,jm)) ...
                                           ./DY(ii,jj) ); 

R3(ii,jj)=R3(ii,jj)+TEMP1 * ...
          ( ((U3(ip,jj)-U3(ii,jj))./DX(ip,jj) - (U3(ii,jj)-U3(im,jj))./DX(ii,jj)) ...
                                           ./DXF(ii,jj) ...
           +((U3(ii,jp)-U3(ii,jj))./DY(ii,jp) - (U3(ii,jj)-U3(ii,jm))./DY(ii,jj)) ...
                                           ./DYF(ii,jj) );                                                 
                
for n=1:N_TH % loop over all scalars
    
% Add FTH to RTH
RTH(ii,jj,n)=RTH(ii,jj,n)+TEMP5*FTH(ii,jj,n);

% Add the explicit part of the Crank-Nicolson diffusive term
RTH(ii,jj,n)=RTH(ii,jj,n)+(TEMP1/PR(n)) * ...
          ( ((TH(ip,jj,n)-TH(ii,jj,n))./DX(ip,jj) - (TH(ii,jj,n)-TH(im,jj,n))./DX(ii,jj)) ...
                                           ./DXF(ii,jj) ...
           +((TH(ii,jp,n)-TH(ii,jj,n))./DY(ii,jp) - (TH(ii,jj,n)-TH(ii,jm,n))./DY(ii,jj)) ...
                                           ./DYF(ii,jj) );
end % End loop over all scalars

                                       
% At this point, the RHS arrays have all been fully constructed
% Solve the implicit systems of equations, starting with TH

for n=1:N_TH % loop over all scalars

% Create some arrays that will be useful for constructing the implicit
% matrices
ii2=2:NXM;
jj2=2:NYM;
ip2=3:NX;
im2=1:NX-2;
jp2=3:NY;
jm2=1:NY-2;
row2=floor((0:(NX-2)*(NY-2)*5-1)/5);
row2=row2+NX+1; % Skip the j=1 rows, which will be handeled with boundary conditions
row2=row2+floor((0:(NX-2)*(NY-2)*5-1)/(5*(NX-2)))+1; % Skip the i=1 rows
row2=row2+floor((0:(NX-2)*(NY-2)*5-1)/(5*(NX-2)))+0; % Skip the i=0 rows

ri=zeros((NX-2)*(NY-2)*5,1);
ci=zeros((NX-2)*(NY-2)*5,1);
vi=zeros((NX-2)*(NY-2)*5,1);
MAT=zeros(NX,NY);
VEC=zeros(NX,NY);

% Here, ri and ci are the column and row indeces of the sparse matrix,
% and vi is the value at that location, i.e. MAT(ri,ci)=vi;
ri(1:(NX-2)*(NY-2)*5,1)=row2';

ci(1:5:(NX-2)*(NY-2)*5,1)=row2(1:5:(NX-2)*(NY-2)*5);
ci(2:5:(NX-2)*(NY-2)*5,1)=row2(1:5:(NX-2)*(NY-2)*5)-1;
ci(3:5:(NX-2)*(NY-2)*5,1)=row2(1:5:(NX-2)*(NY-2)*5)+1;
ci(4:5:(NX-2)*(NY-2)*5,1)=row2(1:5:(NX-2)*(NY-2)*5)-NX;
ci(5:5:(NX-2)*(NY-2)*5,1)=row2(1:5:(NX-2)*(NY-2)*5)+NX;

vi(1:5:(NX-2)*(NY-2)*5,1)=reshape(1-TEMP1*(-1./DX(ii2,jj2)./DXF(ii2,jj2)-1./DX(ip2,jj2)./DXF(ii2,jj2) ...
                          -1./DY(ii2,jp2)./DYF(ii2,jj2)-1./DY(ii2,jj2)./DYF(ii2,jj2)),(NX-2)*(NY-2),1);
vi(2:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DXF(ii2,jj2)./DX(ii2,jj2),(NX-2)*(NY-2),1);
vi(3:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DXF(ii2,jj2)./DX(ip2,jj2),(NX-2)*(NY-2),1);
vi(4:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DY(ii2,jj2)./DYF(ii2,jj2),(NX-2)*(NY-2),1);
vi(5:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DY(ii2,jp2)./DYF(ii2,jj2),(NX-2)*(NY-2),1);

VEC=zeros(NX*NY,1); % Create a blank vector to hold the RHS of the linear system
VEC(row2(1:5:end),1)=reshape(RTH(ii,jj,n),(NX-2)*(NY-2),1);

MAT=sparse(ri,ci,vi,NX*NY,NX*NY);

rk_apply_bc_th; % Apply the boundary conditions


[TEMP,FLAG]=bicgstab(MAT,VEC,iter_tol,iter_max,L_TH,U_TH,reshape(TH(1:NX,1:NY,n),NX*NY,1));
if (FLAG==0) 
    TH(1:NX,1:NY,n)=reshape(TEMP,NX,NY);
else %The system failed to converge with the desired accuracy solve directly instead
    disp('Convergence failed');
    TH(1:NX,1:NY,n)=reshape(MAT\VEC,NX,NY); 
end

end  % END FOR n=1:N_TH

% Now solve for U1
% Here, ri and ci are the same, we only need to change the vi vector
vi(1:5:(NX-2)*(NY-2)*5,1)=reshape(1-TEMP1*(-1./DX(ii2,jj2)./DXF(ii2,jj2)-1./DX(ii2,jj2)./DXF(im2,jj2) ...
                          -1./DY(ii2,jp2)./DYF(ii2,jj2)-1./DY(ii2,jj2)./DYF(ii2,jj2)),(NX-2)*(NY-2),1);
vi(2:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DXF(im2,jj2)./DX(ii2,jj2),(NX-2)*(NY-2),1);
vi(3:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DXF(ii2,jj2)./DX(ii2,jj2),(NX-2)*(NY-2),1);
vi(4:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DY(ii2,jj2)./DYF(ii2,jj2),(NX-2)*(NY-2),1);
vi(5:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DY(ii2,jp2)./DYF(ii2,jj2),(NX-2)*(NY-2),1);

VEC=zeros(NX*NY,1); % Create a blank vector to hold the RHS of the linear system
VEC(row2(1:5:end),1)=reshape(R1(ii,jj),(NX-2)*(NY-2),1);

MAT=sparse(ri,ci,vi,NX*NY,NX*NY);
rk_apply_bc_u1; % Apply the boundary conditions

[TEMP,FLAG]=bicgstab(MAT,VEC,iter_tol,iter_max,L_U1,U_U1,reshape(U1(1:NX,1:NY),NX*NY,1));
if (FLAG==0) 
    U1=reshape(TEMP,NX,NY);
else %The system failed to converge with the desired accuracy solve directly instead
    U1=reshape(MAT\VEC,NX,NY); 
end

% Solve for U2
vi(1:5:(NX-2)*(NY-2)*5,1)=reshape(1-TEMP1*(-1./DX(ii2,jj2)./DXF(ii2,jj2)-1./DX(ip2,jj2)./DXF(ii2,jj2) ...
                          -1./DY(ii2,jj2)./DYF(ii2,jm2)-1./DY(ii2,jj2)./DYF(ii2,jj2)),(NX-2)*(NY-2),1);
vi(2:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DXF(ii2,jj2)./DX(ii2,jj2),(NX-2)*(NY-2),1);
vi(3:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DXF(ii2,jj2)./DX(ip2,jj2),(NX-2)*(NY-2),1);
vi(4:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DY(ii2,jj2)./DYF(ii2,jm2),(NX-2)*(NY-2),1);
vi(5:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DY(ii2,jj2)./DYF(ii2,jj2),(NX-2)*(NY-2),1);

VEC=zeros(NX*NY,1); % Create a blank vector to hold the RHS of the linear system
VEC(row2(1:5:end),1)=reshape(R2(ii,jj),(NX-2)*(NY-2),1);
 
MAT=sparse(ri,ci,vi,NX*NY,NX*NY);
rk_apply_bc_u2; % Apply the boundary conditions

[TEMP,FLAG]=bicgstab(MAT,VEC,iter_tol,iter_max,L_U2,U_U2,reshape(U2(1:NX,1:NY),NX*NY,1));
if (FLAG==0) 
    U2=reshape(TEMP,NX,NY);
else %The system failed to converge with the desired accuracy solve directly instead
    U2=reshape(MAT\VEC,NX,NY); 
end

if (SOLVE_U3) % If we are solving for U3

vi(1:5:(NX-2)*(NY-2)*5,1)=reshape(1-TEMP1*(-1./DX(ii2,jj2)./DXF(ii2,jj2)-1./DX(ip2,jj2)./DXF(ii2,jj2) ...
                          -1./DY(ii2,jp2)./DYF(ii2,jj2)-1./DY(ii2,jj2)./DYF(ii2,jj2)),(NX-2)*(NY-2),1);
vi(2:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DXF(ii2,jj2)./DX(ii2,jj2),(NX-2)*(NY-2),1);
vi(3:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DXF(ii2,jj2)./DX(ip2,jj2),(NX-2)*(NY-2),1);
vi(4:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DY(ii2,jj2)./DYF(ii2,jj2),(NX-2)*(NY-2),1);
vi(5:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DY(ii2,jp2)./DYF(ii2,jj2),(NX-2)*(NY-2),1);

VEC=zeros(NX*NY,1); % Create a blank vector to hold the RHS of the linear system
VEC(row2(1:5:end),1)=reshape(R3(ii,jj),(NX-2)*(NY-2),1);

MAT=sparse(ri,ci,vi,NX*NY,NX*NY);
rk_apply_bc_u3; % Apply the boundary conditions

[TEMP,FLAG]=bicgstab(MAT,VEC,iter_tol,iter_max,L_U3,U_U3,reshape(U3(1:NX,1:NY),NX*NY,1));
if (FLAG==0) 
    U3=reshape(TEMP,NX,NY);
else %The system failed to converge with the desired accuracy solve directly instead
    U3=reshape(MAT\VEC,NX,NY); 
end

end

% Now, we have the new velocity field, but it is not guarenteed to be divergence free

% Perform the projection step of the fractional step method to update the pressure, and remove any divergence in U

% Here, we solve del^2(phi)=div(U)
% Construct the linear system using second order finite differences
% Note that since the pressure and scalars are defined on the same grid, this is very similar to solvin for TH (or U3)
 
VEC=zeros(NX*NY,1); % Create a blank vector to hold the RHS of the linear system
VEC(row2(1:5:end),1)=reshape(((U1(ip2,jj2)-U1(ii2,jj2))./DXF(ii2,jj2) + ...
                           +(U2(ii2,jp2)-U2(ii2,jj2))./DYF(ii2,jj2))/TEMP4,(NX-2)*(NY-2),1);

rk_apply_bc_p; % Apply boundary conditions
 
% The coefficients in the pressure equation don't change in time, so we
% will use a direct solve based on an LU decomposition.  This takes a while
% to set up when we first run the code, but is fast afterwards.
PHI=reshape(U_P\(L_P\(VEC(p_P,:))),NX,NY);

% Update the velocity based on PHI - this `projects' U onto a divergence-free field
% U=U-GRAD(PHI)
U1(ii,jj)=U1(ii,jj)-TEMP4*(PHI(ii,jj)-PHI(im,jj))./DX(ii,jj);
U2(ii,jj)=U2(ii,jj)-TEMP4*(PHI(ii,jj)-PHI(ii,jm))./DY(ii,jj);
% U3 remains the same since d/dz=0

% Apply boundary conditions to the updated velocity field
rk_apply_bc_vel;

% Finally, update the pressure field
% Note that here we divide by H_BAR since it was absorbed into the definition of PHI
P=P+PHI;

% We are now done with this RK STEP

