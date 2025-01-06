function [] = cgDNApviewer( shape, sequence , bp_id)

%--------------------------------------------------------------------------
% cgDNA function: [] = cgDNAviewer( shape, sequence  )
%--------------------------------------------------------------------------
% 3D viewer of cgDNA ground-states, or in generale, 3D viewer of a coarse
% grained DNA fragment given in cgDNA-like internal coordinates. 
%
% Input : 
%
%  shapes    ground-state coordinate vector 
%            in non-dimensional Curves+ form (see Note 1)
%            [size N x 1].
%
%  sequence  DNA sequence in the {A,T,G,C} alphabet corresponding to the
%            cofiguration 'shapes'
%
%  bp_id     Index of the basepair to be used as center frame (i.e
%            baspair(bp_id).R = eye(3) and basepair(bp_id).r = zeros(3,1). 
%
% Note 1:
%
%    The entries of the input variable shapes must be 
%    consistent with the following ordering of the 
%    structural coordinates
%
%     y_1, z_1, ..., y_{nbp-1}, z_{nbp-1}, y_{nbp}
%
%    where for each a=1,2,3,... we have
%
%     y_a = (Buckle,Propeller,Opening,Shear,Stretch,Stagger)_a
%
%     z_a = (Tilt,Roll,Twist,Shift,Slide,Rise)_a.
%
%    For example
%
%     shapes((a-1)*12+1) = Buckle at basepair a
%     shapes((a-1)*12+2) = Propeller at basepair a
%      ...
%     shapes((a-1)*12+6) = Stagger at basepair a
%     shapes((a-1)*12+7) = Tilt at junction a, a+1
%     shapes((a-1)*12+8) = Roll at junction a, a+1
%      ...
%     shapes((a-1)*12+12) = Rise at junction a, a+1.
%
% If you find this code useful, please cite:
%
% D. Petkeviciute, M. Pasi, O. Gonzalez and J.H. Maddocks. 
%  cgDNA: a software package for the prediction of sequence-dependent 
%  coarse-grain free energies of B-form DNA. Submitted (2014). 
%
%--------------------------------------------------------------------------

% extract the 2 extra phosphates group coordinates
% upadte shapes to let rest of viewer unchanged
% shape = shape(1:end-6); % remove the last inter
% p_5 = shape(1:6);
% p_3 = shape(end-5:end);
% shape = shape(7:end-6); % remove extra phosphates

% Compute the absolute coordinates from the internal coordinates
bp_level = frames(shape);

% If the third argument is provided, change the reference frame
if nargin > 2
    if bp_id ~=1
        % Reframe the base and basepair frames with respect to the
        % reference basepair index
        bp_level = reframe(bp_level,bp_id);
    end
end

% Compute 
[ RB, RBc, Pw, Pc ] = sequence3DConf( bp_level, sequence ) ;

% compute positions of extra phosphates groups
scale = 1.5;
dims = size(Pw);
nbp = dims(3)+1;

Pw_temp = zeros(4,3,nbp);
Pw_temp(:,:,1:nbp-1) = Pw;
D = cay(p_5(1:3)) ;
O = p_5(4:6);
Pw_temp(:,:,nbp) = [O' ; (O+scale*D(:,1))' ; (O+scale*D(:,2))' ; (O+scale*D(:,3))'];
Pw = Pw_temp;

Pc_temp = zeros(4,3,nbp);
Pc_temp(:,:,1:nbp-1) = Pc;
Rc = bp_level(nbp).Rc*diag([1,-1,-1]) ;
D = Rc*cay(p_3(1:3)) ;
O = bp_level(nbp).rc + Rc*(p_3(4:6)) ;
Pc_temp(:,:,nbp) = [O' ; (O+scale*D(:,1))' ; (O+scale*D(:,2))' ; (O+scale*D(:,3))'] ;
Pc = Pc_temp;

% Display rigid bodies
display3Dstructure(RB, RBc, Pw, Pc)

% Fix the point of view
view(3)

% Display the axis labels
xlabel('x')
ylabel('y')
zlabel('z')

% Visualize the grid 
grid on 

% Fix the axis properties 
axis equal on

end

function Q = cay(k)

    I = eye(3) ;
    alpha = 1/10 ;
    k = alpha*k ;
    X = [   0   -k(3)  k(2) ;
           k(3)   0   -k(1) ;
          -k(2)  k(1)   0 ] ;
    Q = (I+X)/(I-X) ;

end