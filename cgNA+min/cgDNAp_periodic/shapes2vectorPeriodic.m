function y = shapes2vectorPeriodic(eta,w,etapW,wpW,u,v,etapC,wpC)

%-------------------------------------------------------
% cgDNA function: y = shapes2vector(eta,w,u,v)
%-------------------------------------------------------
% This function re-orders the ground-state coordinates.
%
% Input: 
%
%    eta  list of intra-basepair rotational coords 
%         (Buckle,Propeller,Opening) along molecule 
%         [size nbp x 3]
%
%    w    list of intra-basepair translational coords 
%         (Shear,Stretch,Stagger) along molecule 
%         [size nbp x 3]
%
%    u    list of inter-basepair rotational coords 
%         (Tilt,Roll,Twist) along molecule 
%         [size (nbp-1) x 3]
%
%    v    list of inter-basepair translational coords 
%         (Shift,Slide,Rise) along molecule 
%         [size (nbp-1) x 3].
%
%
% Output:
%
%    y    overall coordinate vector 
%         [size N x 1]
%
%    where N = 12*nbp - 6 and nbp is the length 
%    of the DNA sequence (number of basepairs).
%
%
% Note:
%
%    The entries in the vector y are consistent with the
%    following ordering of the structural coordinates
%
%     x_1, z_1, ..., x_{nbp-1}, z_{nbp-1}, x_{nbp}
%
%    where for each a=1,2,3,... we have
%
%     x_a = (Buckle,Propeller,Opening,Shear,Stretch,Stagger)_a
%
%     z_a = (Tilt,Roll,Twist,Shift,Slide,Rise)_a.
%
%    For example
%
%     y((a-1)*12+1) = Buckle at basepair a
%     y((a-1)*12+2) = Propeller at basepair a
%      ...
%     y((a-1)*12+6) = Stagger at basepair a
%     y((a-1)*12+7) = Tilt at junction a, a+1
%     y((a-1)*12+8) = Roll at junction a, a+1
%      ...
%     y((a-1)*12+12) = Rise at junction a, a+1.
%
%
% If you find this code useful, please cite:
%
% D. Petkeviciute, M. Pasi, O. Gonzalez and J.H. Maddocks. 
%  cgDNA: a software package for the prediction of 
%  sequence-dependent coarse-grain free energies of B-form 
%  DNA. Nucleic Acids Research 2014; doi: 10.1093/nar/gku825.
%
%-------------------------------------------------------

    nbp = size(eta,1);
    N = 24*nbp;
    
    intra = [eta, w];
    inter = [u, v];
    pho_W = [etapW, wpW];
    pho_C = [etapC, wpC];
    
    y = [ pho_W intra pho_C inter ]';
    y = reshape(y, N, 1);
    
end
