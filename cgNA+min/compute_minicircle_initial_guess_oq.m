function [pos, quat] = compute_minicircle_initial_guess_oq(helix_uv,m,n,phimin)
% Initialization of output variables (origin and quaternion of initial minicircle shape for initial input register angle)
pos  = zeros(n+1,3);
quat = zeros(n+1,4);

phi       =  phimin;        % initial register angle
v         =  helix_uv(4:6); % helix inter translation in A
u         =  helix_uv(1:3); % helix inter rotation (Cayley vector) in rad/5

dpsi      =  2*atan(norm(u)/10);    % in radian
gamma     =  u/norm(u);
dz        =  dot(v,gamma);

rh        =  norm(v - dz*gamma)/(2*sin(dpsi/2));   % helix radius
rc        =  n*dz/(2*pi);                          % torus/circle radius

rotQ      =  [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1]; % rotation matrix corresponding to register angle phi
qphi      =  [0; 0; sin(phi/2);cos(phi/2)];   % quaternion of rotQ
qphiT     =  [-0;-0;-sin(phi/2);cos(phi/2)];  % quaternion coresponding to transpose of rotQ

angphi = 2*pi*m/n;
angal  = 2*pi/n;

for i=0:n
    rci    =  (rc + rh*cos(angphi*i));
    posc   =  [rh*sin(angphi*i); rci*cos(angal*i)-rc; rci*sin(angal*i)];  % bp origins on toroid

    qc1 =   sin(angal*i/2)*cos(angphi*i/2);
    qc2 = - sin(angal*i/2)*sin(angphi*i/2);
    qc3 =   cos(angal*i/2)*sin(angphi*i/2);
    qc4 =   cos(angal*i/2)*cos(angphi*i/2);
    qc  =   [qc1;qc2;qc3;qc4];                      % bp quaternions on toroid

    % Including register angle
    poscf   =  rotQ*posc;
    [qstar] =  quaternion_product(qc,qphiT);
    [qcf]   =  quaternion_product(qphi,qstar);

    pos(i+1,1:3)   =  poscf';  % minicircle: origins of initial guess bp frames for given register angle
    quat(i+1,1:4)  =  qcf';    % minicircle: quaternions of initial guess bp frames for given register angle
end

end

