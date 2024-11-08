clc
clear all
close all

%% define any input sequence
seq      = 'TAATGCATGTACATTTATGGTGGCCCACGTGAGACCCAAACGGGTGTGATCGTGTGACTTGTACATAGGGTGAAATATGTGGGCCGTGCCGCAG';   % any seq S
seq_name = "random (94bp)";
nbp      = length(seq);     % length of sequence in base-pairs
link     = round(nbp/10.5); % standard link computation
dellk    = 0;               % adding/subtracting extra link (optional input)
Lk       = link + dellk;    % initial link for which computation will start

disp(sprintf('Initializing minimization of sequence: "%s"', seq_name));
disp(sprintf('Number of base pairs:                   %d', nbp));
disp(sprintf('Linking number of the initial guess:    %d', Lk));

%% defining global variables
global whats
global mmats
global nmats
global omats
global pmats
global qmats
global q4_at_1

%% reconstructing cgNA+ periodic groundstate & periodic stiffness matrix
disp('Preparing cgNA+ model periodic groundstate vector and stiffness matrix');

addpath('./cgDNAp_periodic')
params = load('Prmset_cgDNA+_CGF_10mus_int_12mus_ends.mat'); % cgNA+ model parameter set

[whats, stiff] = constructSeqParmsPeriodic(seq, params);

%% Decomposing sparse stiffness matrix into nonzero sub-blocks
mmats = zeros(nbp,  6,  6);
nmats = zeros(nbp, 18, 18);
omats = zeros(nbp,  6, 18);
pmats = zeros(nbp, 18,  6);
qmats = zeros(nbp, 18, 18);
for i=1:nbp-1
    for j=1:6
        for k=1:6
            mmats(i,j,k) = stiff(18+24*(i-1)+j, 18+24*(i-1)+k);
        end
    end
    for j=1:6
        for k=1:18
            omats(i,j,k) = stiff(18+24*(i-1)+j, 24+24*(i-1)+k);
            pmats(i,k,j) = stiff(   24*(i-1)+k, 18+24*(i-1)+j);
        end
    end
    for j=1:18
        for k=1:18
            nmats(i,j,k) = stiff(  24*(i-1)+j,   24*(i-1)+k);
            qmats(i,j,k) = stiff(  24*(i-1)+j,24+24*(i-1)+k);
        end
    end
end
i=nbp;
for j=1:6
    for k=1:6
        mmats(i,j,k) = stiff(18+24*(i-1)+j, 18+24*(i-1)+k);
    end
end
for j=1:6
    for k=1:18
        omats(i,j,k) = stiff(18+24*(i-1)+j, k);
        pmats(i,k,j) = stiff(  24*(i-1)+k, 18+24*(i-1)+j);
    end
end
for j=1:18
    for k=1:18
        nmats(i,j,k) = stiff(  24*(i-1)+j,   24*(i-1)+k);
        qmats(i,j,k) = stiff(  24*(i-1)+j, k);
    end
end

%%  computing initial guess minicircle of input seq and link
disp('Computing an initial guess minicircle');

lambda  = [-2500:0.5:2500]; % range for multiplier
tol     = 1.0e-2;           % Set the tolerance to find the given link

% Create matrix P so that Pz = whats where z is helicoidal configuration of a given S with prescribed Lk
P = sparse(24*nbp,18*nbp+6);
for i=1:2*nbp
    if mod(i,2)==0
        st=i*9+(i/2-1)*6;
        P(st+1:st+6,end-5:end) = eye(6);
    else
        st=(i-1)*12;
        en=(i-1)*9;
        P(st+1:st+18,en+1:en+18) = eye(18);
    end
end
A = (P')*stiff*P;
b = (P')*stiff*whats;

% create matrix E required for linear solve with multiplier
E                                      = sparse(18*nbp+6,18*nbp+6);
E(18*nbp+1:18*nbp+3,18*nbp+1:18*nbp+3) = eye(3);

% Create a vector of link numbers associated with each lambda for the full range of multiplier
m = zeros(1,length(lambda));

for i=1:length(lambda)
    % Define the new matrix A_new to solve constraint problem
    A_new = (A-lambda(i)*E);
    % Check it is positive definite to solve the linear system
    [R,posdef(i)]=chol(A_new); % posdef=0 if A_new is positive definite
    if posdef(i)~=0
        break
    end
    % Solve the linear system using the Cholesky decomposition
    tmp = R' \ b;
    z = R \ tmp;  % Helicoidal configuration vector
    w = P*z;      % Reconstruct the full internal vector
    % Evaluate the squared norm of the Caylay vector u
    u = norm(z(18*nbp+1:18*nbp+3));
    % Compute the angle of rotation in radians from nrom of  uniform Caley vector u
    theta = 2*atan(u/10);
    % Compute the number of links associated to each norm(u) corresponding to each lambda
    m(i) = theta*nbp/(2*pi);
    % Storing the unform inter coordinates u & v
    uv = z(end-5:end);
    % Check if we found the correct number of links
    if m(i)<Lk+tol && m(i)>Lk-tol
        whelix   = w;
        helix_uv = uv;
        break;
    end
end

% Define value of initial register angle (0<phimin>2pi)
phimin    = 0;
% Compute initial position and quaterions corresponding to minicircle using uniform helical inters of prescribed link for any initial register angle
[os, qs] = compute_minicircle_initial_guess_oq(helix_uv,Lk,nbp,phimin);

%% minimizing cgDNA+ energy
disp('Starting minimization');

% Compute shifted coordinates for base-pair label coordinates
whelix_shifted = whelix - whats;  % required just for intras and phosphates orders or [pw intra pc inters]
intras_with_p  = zeros(nbp, 18);
for i=1:nbp
    intras_with_p(i,1:18) = whelix_shifted(24*(i-1)+1:24*(i-1)+18);
end
%  Construct zvec vector: complete vector (except first o & q) containing all the variables (absolute origins and quaternions, plus relative base-pair lavel coordinates) of initial guess minicircle
zvec       = zeros(25*nbp-7,1);
zvec(1:18) = transpose(intras_with_p(1,:));
for i=2:nbp
    zvec(18+25*(i-2)+1:18+25*(i-2)+25) = transpose([intras_with_p(i,:) os(i,:) qs(i,:)]);
end
q4_at_1  = qs(nbp+1,4);

% Using fminunc with trust-region algorithm and asociated parameters
options                     =  optimset('Display','iter','MaxIter',5000,'GradObj','on','Hessian','on','TolFun',1e-16,'Algorithm','trust-region','TolX',1e-16);
[zeq,eneq,exitflag,output]  =  fminunc('discrete_dna_penalty_en_grad_hess',zvec,options); % here function "discrete_dna_penalty_en_grad_hess" returns energy, gradient and hessian

disp('Minimization stops with following details');
output

% Computing energy, gradient and hessian at stoped configration
[en,grad,hess]    =  discrete_dna_penalty_en_grad_hess(zeq);
% Computing eigenvalues of Hessian
[V,D]    = eigs(hess,length(hess));
[d,ind]  = sort(diag(D));
Ds       = D(ind,ind);
for ij=1:length(hess)
    eigval(ij) = Ds(ij,ij);
end
disp('Five smallest eigenvalues of the Hessian at found energy configration:');
disp(sprintf('                                          %e\n', eigval(1:5)));

%% Compute cgNA+ model internal coordinates, linking number and cgNA+ energy of initial and final minicircle shape
[inters1,intras1] = compute_internal_variables(zvec);  % inters and base-pare lavel coordinates corresponding to initial guess minicircle shape
[inters2,intras2] = compute_internal_variables(zeq);   % inters and base-pare lavel coordinates corresponding to corresponding to final minicircle shape
%  computing internal coordinates of minicircles  %%%%%%%%%
initial_guess_minicircle_shape = zeros(24*nbp,1);
final_minicircle_shape         = zeros(24*nbp,1);
for i=1:nbp
    initial_guess_minicircle_shape(24*(i-1)+1 :24*(i-1)+18) = intras1(:,i);
    initial_guess_minicircle_shape(24*(i-1)+19:24*(i-1)+24) = inters1(:,i);
    final_minicircle_shape(24*(i-1)+1 :24*(i-1)+18)         = intras2(:,i);
    final_minicircle_shape(24*(i-1)+19:24*(i-1)+24)         = inters2(:,i);
end

% Computing the minicircle cgNA+ model energy associated to initial minicircle shape
cgNA_initial_minicircle_energy = 0.5*(initial_guess_minicircle_shape-whats)'*stiff*(initial_guess_minicircle_shape-whats);
disp(sprintf('cgNA+ initial guess minicircle energy:      %f RT', 0.593*cgNA_initial_minicircle_energy));
% Computing the minicircle cgNA+ model energy associated to final minicircle shape
cgNA_minicircle_energy = 0.5*(final_minicircle_shape-whats)'*stiff*(final_minicircle_shape-whats);
disp(sprintf('cgNA+ energy of found minicircle:      %f RT', 0.593*cgNA_minicircle_energy));

% Computing the linking number of final minicircle shape
bp_level = framesPeriodic(final_minicircle_shape);
% Take phosphate position and close the curves
rpw       = [bp_level([1:nbp, 1]).rpw];
rpc       = [bp_level([1:nbp, 1]).rpc];
Lk_final  = compute_lk(rpw, rpc);
disp(sprintf('Linking number of the found minicircle:     %d', Lk_final));

%% Plotting
do_title = true; do_save = false;is_diff=false; guess_nb=0;
% Comparing periodic groundstate of linear fragment, initial guess minicircle shape and found minicircle shape
plot2D_cgDNAp_coordinates(whats, initial_guess_minicircle_shape, final_minicircle_shape, do_title, seq_name, guess_nb, do_save, is_diff);
% 3d plot of found/final minicircle shape
plot3D_cgNAp_minicircle(final_minicircle_shape, seq, do_title, seq_name, guess_nb, do_save)

%% END

