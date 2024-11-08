function [energy,grad,hess] = discrete_dna_penalty_en_grad_hess(zvec)
global whats
global mmats
global nmats
global omats
global pmats
global qmats
global q4_at_1
% We have 18*nbp intra (base-pair level coordinates) unknowns and 7*(nbp-2) unknowns (o,q),
%    with the first (o,q)=(0,e4) and the last (0,\pm e4) depending on initial guess
% Ordering of zvec is intra1, then intra2,oq2,intra3,oq3,...oqn-2,intran-1,oqn-1 then intran, oqn
s              = size(zvec);
zlen           = s(1);
nbp            = (zlen+7)/25;
penalty_weight = 100;   % coefficient for the n-2 penalties on (|q|^2-1)^2

% Read (o,q) values from zvec -- and assign values at ends
q = zeros(4,nbp+1);
o = zeros(3,nbp+1);
for i=2:nbp
    o(:,i) = zvec(36+25*(i-2)+1:36+25*(i-2)+3);
    q(:,i) = zvec(36+25*(i-2)+4:36+25*(i-2)+7);
end
q(4,1)     = 1;
q(4,nbp+1) = q4_at_1;

% Read intra values from zvec
intras      = zeros(18,nbp);
intradiffs  = zeros(18,nbp);
for i=2:nbp-1
    intradiffs(:,i) = zvec(18+25*(i-2)+1:18+25*(i-2)+18);
end
intradiffs(:,1)   = zvec(1:18);
intradiffs(:,nbp) = zvec(zlen-24:zlen-7);
for i=1:nbp
    intras(:,i) = intradiffs(:,i)+whats(24*(i-1)+1:24*(i-1)+18);
end

% Set up key constant matrices
id = eye(4);
b  = zeros(3,4,4);
b(1,1,4) = 1; b(1,2,3)  = 1; b(1,3,2)=-1; b(1,4,1) = -1;
b(2,1,3) = -1; b(2,2,4) = 1; b(2,3,1) = 1; b(2,4,2) = -1;
b(3,1,2) = 1; b(3,2,1)  = -1; b(3,3,4) = 1; b(3,4,3) = -1;

% Store values of B1q,B2q,B3q in triply-indexed array bq
bq = zeros(3,4,nbp+1);
for i=1:3
    bi = zeros(4,4);
    for j1 = 1:4
        for j2 = 1:4
            bi(j1,j2) = b(i,j1,j2);
        end
    end
    for k = 1:nbp+1
        dum = bi*q(:,k);
        for j1=1:4
            bq(i,j1,k) = dum(j1);
        end
    end
end

% Compute the inters
cay    = zeros(3,nbp);
tr     = zeros(3,nbp);
inters = zeros(6,nbp);
dfac   = zeros(1,nbp);
for i=1:nbp
    dfac(i) = q(:,i+1)'*q(:,i);
    for k=1:3
        bkq = zeros(4,1);
        for j1 = 1:4
            bkq(j1)=bq(k,j1,i);
        end
        cay(k,i)    = 10/dfac(i)*q(:,i+1)'*bkq;
        inters(k,i) = cay(k,i);
    end
    d_o = o(:,i+1)-o(:,i);

    if dot(q(:,i+1),q(:,i))<0
        sign=-1;
    else
        sign=1;
    end
    dirs = compute_ds(q(:,i+1)+sign*q(:,i));

    for k=1:3
        tr(k,i)       = d_o'*dirs(:,k);
        inters(k+3,i) = tr(k,i);
    end
end
% Store values of inter-interhat
interdiffs = zeros(6,nbp);
for i=1:nbp
    interdiffs(:,i) = inters(:,i) - whats(18+24*(i-1)+1:18+24*(i-1)+6);
end

%% Compute the energy -- M contribs, then N, then O, then P, then Q
%   (no factor of 1/2 for O,P,Q since contribs from above and below diag)
energy       = 0.0;
mat_6x6      = zeros(6,6);
mat_6x18     = zeros(6,18);
mat_18x6     = zeros(18,6);
mat_18x18    = zeros(18,18);
vec_6x1      = zeros(6,1);
vec2_6x1     = zeros(6,1);
vec_18x1     = zeros(18,1);
vec2_18x1    = zeros(18,1);

% M
for i=1:nbp
    for j1=1:6
        vec_6x1(j1) = interdiffs(j1,i);
        for j2=1:6
            mat_6x6(j1,j2) = mmats(i,j1,j2);
        end
    end
    energy       =  energy+vec_6x1'*(mat_6x6*vec_6x1)/2;
end
% N
for i=1:nbp
    for j1=1:18
        vec_18x1(j1)  =  intradiffs(j1,i);
        for j2=1:18
            mat_18x18(j1,j2) = nmats(i,j1,j2);
        end
    end
    energy = energy+vec_18x1'*(mat_18x18*vec_18x1)/2;
end
% O
for i=1:nbp
    if i<nbp
        for j1=1:6
            vec_6x1(j1)   =  interdiffs(j1,i);
            for j2=1:18
                mat_6x18(j1,j2) = omats(i,j1,j2);
            end
        end
        for j1=1:18
            vec2_18x1(j1)  =  intradiffs(j1,i+1);
        end
    else % i = nbp
        for j1=1:6
            vec_6x1(j1)   =  interdiffs(j1,nbp);
            for j2=1:18
                mat_6x18(j1,j2) = omats(nbp,j1,j2);
            end
        end
        for j1=1:18
            vec2_18x1(j1)  =  intradiffs(j1,1);
        end
    end
    energy  =  energy+vec_6x1'*(mat_6x18*vec2_18x1);
end
% P
for i=1:nbp
    for j1=1:18
        vec_18x1(j1)  = intradiffs(j1,i);
        for j2=1:6
            mat_18x6(j1,j2) = pmats(i,j1,j2);
        end
    end
    for j1=1:6
        vec2_6x1(j1) = interdiffs(j1,i);
    end
    energy = energy+vec_18x1'*(mat_18x6*vec2_6x1);
end
% Q
for i=1:nbp
    if i<nbp
        for j1=1:18
            vec_18x1(j1)  = intradiffs(j1,i);
            vec2_18x1(j1) = intradiffs(j1,i+1);
            for j2=1:18
                mat_18x18(j1,j2) = qmats(i,j1,j2);
            end
        end
    else % i=nbp
        for j1=1:18
            vec_18x1(j1)  = intradiffs(j1,nbp);
            vec2_18x1(j1) = intradiffs(j1,1);
            for j2=1:18
                mat_18x18(j1,j2) = qmats(i,j1,j2);
            end
        end
    end
    energy = energy+vec_18x1'*(mat_18x18*vec2_18x1);
end

% Add penalty factors for each |q|^2
for i=2:nbp
    energy = energy+penalty_weight*(q(:,i)'*q(:,i)-1)^2;
end

%% Fill in gradient vector
grad = zeros(zlen,1);

% Intra slots
for i=1:nbp
    contrib = zeros(18,1);
    % N
    for j1=1:18
        vec_18x1(j1) = intradiffs(j1,i);
        for j2=1:18
            mat_18x18(j1,j2) = nmats(i,j1,j2);
        end
    end
    contrib = contrib+mat_18x18*vec_18x1;
    % Q
    if i<nbp
        for j1=1:18
            vec_18x1(j1)  =  intradiffs(j1,i+1);
            for j2=1:18
                mat_18x18(j1,j2)  =  qmats(i,j1,j2);
            end
        end
        contrib = contrib+mat_18x18*vec_18x1;
    else % i=nbp
        for j1=1:18
            vec_18x1(j1)  =  intradiffs(j1,1);
            for j2=1:18
                mat_18x18(j1,j2)  =  qmats(i,j1,j2);
            end
        end
        contrib = contrib+mat_18x18*vec_18x1;
    end
    if i>1
        for j1=1:18
            vec_18x1(j1) = intradiffs(j1,i-1);
            for j2=1:18
                mat_18x18(j1,j2)  =  qmats(i-1,j1,j2);
            end
        end
        contrib = contrib+transpose(mat_18x18)*vec_18x1;
    else % i=1
        for j1=1:18
            vec_18x1(j1) = intradiffs(j1,nbp);
            for j2=1:18
                mat_18x18(j1,j2)  =  qmats(nbp,j1,j2);
            end
        end
        contrib = contrib+transpose(mat_18x18)*vec_18x1;
    end
    % P
    for j1=1:6
        vec_6x1(j1) = interdiffs(j1,i);
        for j2=1:18
            mat_18x6(j2,j1) = pmats(i,j2,j1);
        end
    end
    contrib = contrib+mat_18x6*vec_6x1;
    % O
    if i>1
        for j1=1:6
            vec_6x1(j1) = interdiffs(j1,i-1);
            for j2=1:18
                mat_6x18(j1,j2)  =  omats(i-1,j1,j2);
            end
        end
        contrib  =  contrib+transpose(mat_6x18)*vec_6x1;
    else % i=1
        for j1=1:6
            vec_6x1(j1) = interdiffs(j1,nbp);
            for j2=1:18
                mat_6x18(j1,j2)  =  omats(nbp,j1,j2);
            end
        end
        contrib  =  contrib+transpose(mat_6x18)*vec_6x1;
    end
    %%%%%%%%%%%% store values in grad vector
    if i<2
        grad(1:18) = contrib;
    else
        grad(18+25*(i-2)+1:18+25*(i-2)+18)  =  contrib;
    end
end

% Vectors that contribute to inter slots of grad (and to Hessian)
dedz = zeros(6,nbp);

for i=1:nbp
    contrib = zeros(6,1);
    % P
    for j1=1:18
        vec_18x1(j1) = intradiffs(j1,i);
        for j2=1:6
            mat_18x6(j1,j2) = pmats(i,j1,j2);
        end
    end
    contrib  =  contrib+transpose(mat_18x6)*vec_18x1;
    % M
    for j1=1:6
        vec_6x1(j1)  =  interdiffs(j1,i);
        for j2=1:6
            mat_6x6(j1,j2)  =  mmats(i,j1,j2);
        end
    end
    contrib  =  contrib+mat_6x6*vec_6x1;
    % O
    if i<nbp
        for j1=1:18
            vec_18x1(j1)  =  intradiffs(j1,i+1);
            for j2=1:6
                mat_6x18(j2,j1)  =  omats(i,j2,j1);
            end
        end
        contrib   =  contrib+mat_6x18*vec_18x1;
    else % i=nbp
        for j1=1:18
            vec_18x1(j1)  =  intradiffs(j1,1);
            for j2=1:6
                mat_6x18(j2,j1)  =  omats(i,j2,j1);
            end
        end
        contrib   =  contrib+mat_6x18*vec_18x1;
    end
    dedz(:,i) =  contrib;
end

% Matrices d(z^a)/d(o,q)^a that contribute to inter slots of grad (and to Hessian)
% Note that a=1 is left as all zeroes; should not be used
dzdoq = zeros(6,7,nbp);

for a=2:nbp
    for j1=1:3
        for j2=1:4
            dzdoq(j1,j2+3,a)  =  (-10*bq(j1,j2,a+1)-cay(j1,a)*q(j2,a+1))/dfac(a);
        end
    end
    dirs  = compute_ds(q(:,a)+q(:,a+1));
    ddirs = compute_dds(q(:,a)+q(:,a+1));
    odiff = o(:,a+1)-o(:,a);
    for j1=1:3
        for j2=1:3
            dzdoq(j1+3,j2,a)  =  -dirs(j2,j1);
        end
    end
    for j1=1:3
        for j2=1:4
            dzdoq(j1+3,j2+3,a)  =  ddirs(j1,1,j2)*odiff(1)+ddirs(j1,2,j2)*odiff(2)+ddirs(j1,3,j2)*odiff(3);
        end
    end
end
% Matrices d(z^(a-1))/d(o,q)^a that contribute to inter slots of grad (and to Hessian)
% Note that a=1 is left as all zeroes; should not be used
dzprevdoq = zeros(6,7,nbp);

for a=2:nbp
    for j1=1:3
        for j2=1:4
            dzprevdoq(j1,j2+3,a)  =  (10*bq(j1,j2,a-1)-cay(j1,a-1)*q(j2,a-1))/dfac(a-1);
        end
    end
    dirs   =  compute_ds(q(:,a-1)+q(:,a));
    ddirs  =  compute_dds(q(:,a-1)+q(:,a));
    odiff  =  o(:,a) - o(:,a-1);
    for j1=1:3
        for j2=1:3
            dzprevdoq(j1+3,j2,a) = dirs(j2,j1);
        end
    end
    for j1=1:3
        for j2=1:4
            dzprevdoq(j1+3,j2+3,a)  =  ddirs(j1,1,j2)*odiff(1)+ddirs(j1,2,j2)*odiff(2)+ddirs(j1,3,j2)*odiff(3);
        end
    end
end


% Inter slots
mat67 = zeros(6,7);
for a=2:nbp   % The index of (o,q)
    for j1=1:6
        vec_6x1(j1)  =  dedz(j1,a-1);
        for j2=1:7
            mat67(j1,j2)  =  dzprevdoq(j1,j2,a);
        end
    end
    contrib  =  transpose(mat67)*vec_6x1;
    for j1=1:6
        vec_6x1(j1) = dedz(j1,a);
        for j2=1:7
            mat67(j1,j2) = dzdoq(j1,j2,a);
        end
    end
    contrib                            =  contrib+transpose(mat67)*vec_6x1;
    % Add term coming from penalty weight
    qthis                              =  q(:,a);
    fac                                =  4*penalty_weight*(qthis'*qthis-1);
    contrib                            =  contrib+fac*[0;0;0;qthis];
    grad(36+25*(a-2)+1:36+25*(a-2)+7)  =  contrib;
end


%% Setup for building a sparse Hessian (nentries = number of nonzeros
%    we'll put in the Hessian)
nintraintra = 324*nbp+324*(nbp)+324*(nbp);  % (a,a)  (a,a-1)  (a,a+1)
nintrainter = 108*(nbp-1)+108*(nbp-1)+108*(nbp-1);  % (a,a)  (a,a-1)  (a,a+1)
ninterintra = nintrainter;
ninterinter = 49*(nbp-1)+49*(nbp-2)+49*(nbp-2); % (a,a) (a,a-1) (a,a+1)
nentries    = nintraintra+nintrainter+ninterintra+ninterinter;
I           = zeros(nentries,1);
J           = zeros(nentries,1);
X           = zeros(nentries,1);

nset = 0;   % tracks which nonzero entry we are inputting

% (intra,intra) entries of hessian
for a=1:nbp  % intra (a,a) slots
    if a==1
        rowrange = 1:18;
        colrange = rowrange;
    else
        rowrange = 19+25*(a-2):36+25*(a-2);
        colrange = rowrange;
    end
    for j=1:18
        for k=1:18
            nset    = nset+1;
            I(nset) = rowrange(j);
            J(nset) = colrange(k);
            X(nset) = nmats(a,j,k);
        end
    end
end
for a=1:nbp  % intra (a,a+1) slots (and (a,a-1) done by symmetry)
    if a==1
        rowrange = 1:18;
        colrange = rowrange+18;
    elseif a==nbp
        rowrange = 19+25*(a-2):36+25*(a-2);
        colrange = 1:18;
    else
        rowrange = 19+25*(a-2):36+25*(a-2);
        colrange = rowrange+25;
    end
    for j=1:18
        for k=1:18
            nset    = nset+1;
            I(nset) = rowrange(j);
            J(nset) = colrange(k);
            X(nset) = qmats(a,j,k);
            nset    = nset+1;
            I(nset) = J(nset-1);
            J(nset) = I(nset-1);
            X(nset) = X(nset-1);
        end
    end
end


mat1 = zeros(18,6);
mat2 = zeros(6,7);
mat3 = zeros(18,6);
mat4 = zeros(6,7);
% (intra,inter) and (inter,intra) entries of hessian
for a=2:nbp % (intra,inter) (a,a)   and (inter,intra) (a,a) by symmetry
    rowrange = 19+25*(a-2):36+25*(a-2);
    colrange = rowrange(1)+18:rowrange(1)+24;
    for j=1:18
        for k=1:6
            mat1(j,k) = pmats(a  ,j,k);
            mat3(j,k) = omats(a-1,k,j);
        end
    end
    for j=1:6
        for k=1:7
            mat2(j,k) = dzdoq(j,k,a);
            mat4(j,k) = dzprevdoq(j,k,a);
        end
    end
    result = mat1*mat2+mat3*mat4;
    for j=1:18
        for k=1:7
            nset    = nset+1;
            I(nset) = rowrange(j);
            J(nset)=colrange(k);
            X(nset)=result(j,k);
            nset    = nset+1;
            I(nset) = J(nset-1);
            J(nset) = I(nset-1);
            X(nset) = X(nset-1);
        end
    end
end


for a=1:nbp-1 % (intra,inter) (a,a+1)  and (inter,intra) (a+1,a) by symm
    if a==1
        rowrange = 1:18;
        colrange = rowrange(1)+36:rowrange(1)+42;
    else
        rowrange = 19+25*(a-2):36+25*(a-2);
        colrange = rowrange(1)+43:rowrange(1)+49;
    end
    for j=1:18
        for k=1:6
            mat1(j,k) = pmats(a,j,k);
        end
    end
    for j=1:6
        for k=1:7
            mat2(j,k) = dzprevdoq(j,k,a+1);
        end
    end
    result=mat1*mat2;
    for j=1:18
        for k=1:7
            nset    = nset+1;
            I(nset) = rowrange(j);
            J(nset) = colrange(k);
            X(nset) = result(j,k);
            nset    = nset+1;
            I(nset) = J(nset-1);
            J(nset) = I(nset-1);
            X(nset) = X(nset-1);
        end
    end
end



for a=3:nbp+1 % (intra,inter) (a,a-1) and (inter,intra) (a-1,a) by symm
    if a==nbp+1
        rowrange = 1:18;
        colrange = 12+25*(a-2):18+25*(a-2);
    else % a < nbp+1
        rowrange = 19+25*(a-2):36+25*(a-2);
        colrange = rowrange(1)-7:rowrange(1)-1;
    end
    for j=1:18
        for k=1:6
            mat1(j,k) = omats(a-1,k,j);
        end
    end
    for j=1:6
        for k=1:7
            mat2(j,k) = dzdoq(j,k,a-1);
        end
    end
    result = mat1*mat2;
    for j=1:18
        for k=1:7
            nset    = nset+1;
            I(nset) = rowrange(j);
            J(nset) = colrange(k);
            X(nset) = result(j,k);
            nset    = nset+1;
            I(nset) = J(nset-1);
            J(nset) = I(nset-1);
            X(nset) = X(nset-1);
        end
    end
end

% Some new matrices to prepare for (inter,inter) derivs
d2thdq2       = zeros(nbp  ,3,4,4);  % a=1 left as zero (will not be used)
d2thdqnext2   = zeros(nbp-1,3,4,4);
d2thdqnextdq  = zeros(nbp-1,3,4,4);  % a=1 left as zero (will not be used)
dodotd2s      = zeros(nbp  ,3,4,4);
for a=1:nbp
    qthis   =  q(:,a);
    qnext   =  q(:,a+1);
    dotprod =  transpose(qnext)*qthis;
    othis   =  o(:,a);
    onext   =  o(:,a+1);
    d2s     =  compute_drdotd2ds(qthis+qnext,onext-othis);
    for j=1:3
        caythis  =  cay(j,a);
        bqthis   =  zeros(4,1);
        bqnext   =  zeros(4,1);
        for k=1:4
            bqthis(k) = bq(j,k,a);
            bqnext(k) = bq(j,k,a+1);
        end
        %  (a,a) deriv
        mat = 2*caythis*(qnext*qnext')+10*(bqnext*qnext')+10*(qnext*bqnext');
        mat = mat/(dotprod*dotprod);
        if a>1
            for k1=1:4
                for k2=1:4
                    d2thdq2(a,j,k1,k2) = mat(k1,k2);
                end
            end
        end
        % (a+1,a+1) deriv
        mat  =  2*caythis*(qthis*qthis')-10*(bqthis*qthis')-10*(qthis*bqthis');
        mat  =  mat/(dotprod*dotprod);
        if a<nbp
            for k1=1:4
                for k2=1:4
                    d2thdqnext2(a,j,k1,k2) = mat(k1,k2);
                end
            end
        end
        % (a+1,a) deriv
        mat  =  2*caythis*(qnext*qthis')-10*(qnext*bqthis')+10*(bqnext*qthis');
        mat  =  mat/(dotprod*dotprod);
        if a<nbp && a>1
            for k1=1:4
                for k2=1:4
                    d2thdqnextdq(a,j,k1,k2) = mat(k1,k2)+(-10*b(j,k1,k2)-caythis*id(k1,k2))/dotprod;
                end
            end
        end
        % dot products of do with second derivs of dj
        for k1=1:4
            for k2=1:4
                dodotd2s(a,j,k1,k2) = d2s(j,k1,k2);
            end
        end
    end
end

% (inter,inter) entries of hessian
mat1 = zeros(7,7);
mat2 = zeros(7,6);
mat3 = zeros(6,6);
mat4 = zeros(6,7);
for a=2:nbp-1 % (a,a+1)   and (a+1,a) by symmetry
    rowrange = 37+25*(a-2):43+25*(a-2);
    colrange = rowrange(1)+25:rowrange(1)+31;
    vec=zeros(6,1);
    for j=1:6
        vec(j) = dedz(j,a);
    end
    ddirs = compute_dds(q(:,a)+q(:,a+1));
    for n=1:3
        mat=zeros(7,6);
        for k1=4:7
            for k2=4:6
                mat(k1,k2) = ddirs(k2-3,n,k1-3);
            end
        end
        mat1(:,n) = mat*vec;
    end
    for n=4:7
        mat=zeros(7,6);
        for k1=1:3
            for k2=4:6
                mat(k1,k2) = -ddirs(k2-3,k1,n-3);
            end
        end
        for k1=4:7
            for k2=4:6
                mat(k1,k2)=dodotd2s(a,k2-3,n-3,k1-3);
            end
        end
        for k1=4:7
            for k2=1:3
                mat(k1,k2) = d2thdqnextdq(a,k2,k1-3,n-3);
            end
        end
        mat1(:,n)=mat*vec;
    end
    for k1=1:6
        for k2=1:6
            mat3(k1,k2) = mmats(a,k1,k2);
        end
    end
    for k1=1:6
        for k2=1:7
            mat4(k1,k2) = dzprevdoq(k1,k2,a+1);
            mat2(k2,k1) = dzdoq(k1,k2,a);
        end
    end
    result = mat1+mat2*mat3*mat4;
    for j=1:7
        for k=1:7
            nset    = nset+1;
            I(nset) = rowrange(j);
            J(nset) = colrange(k);
            X(nset) = result(j,k);
            nset    = nset+1;
            I(nset) = J(nset-1);
            J(nset) = I(nset-1);
            X(nset) = X(nset-1);
        end
    end
end

for a=2:nbp % (a,a)
    rowrange = 37+25*(a-2):43+25*(a-2);
    colrange = rowrange;
    vec      = zeros(6,1);
    for j=1:6
        vec(j) = dedz(j,a-1);
    end
    ddirs = compute_dds(q(:,a-1)+q(:,a));
    for n=1:3
        mat=zeros(7,6);
        for k1=4:7
            for k2=4:6
                mat(k1,k2) = ddirs(k2-3,n,k1-3);
            end
        end
        mat1(:,n)=mat*vec;
    end
    for n=4:7
        mat=zeros(7,6);
        for k1=1:3
            for k2=4:6
                mat(k1,k2) = ddirs(k2-3,k1,n-3);
            end
        end
        for k1=4:7
            for k2=4:6
                mat(k1,k2) = dodotd2s(a-1,k2-3,n-3,k1-3);
            end
        end
        for k1=4:7
            for k2=1:3
                mat(k1,k2) = d2thdqnext2(a-1,k2,k1-3,n-3);
            end
        end
        mat1(:,n) = mat*vec;
    end
    for k1=1:6
        for k2=1:6
            mat3(k1,k2) = mmats(a-1,k1,k2);
        end
    end
    for k1=1:6
        for k2=1:7
            mat4(k1,k2) = dzprevdoq(k1,k2,a);
            mat2(k2,k1) = dzprevdoq(k1,k2,a);
        end
    end
    result = mat1+mat2*mat3*mat4;
    % Second half of formula
    vec = zeros(6,1);
    for j=1:6
        vec(j)=dedz(j,a);
    end
    if a<nbp
        ddirs = compute_dds(q(:,a)+q(:,a+1));
    else
        ddirs = compute_dds(q(:,a)+q(:,1));
    end
    for n=1:3
        mat=zeros(7,6);
        for k1=4:7
            for k2=4:6
                mat(k1,k2) = -ddirs(k2-3,n,k1-3);
            end
        end
        mat1(:,n) = mat*vec;
    end
    for n=4:7
        mat=zeros(7,6);
        for k1=1:3
            for k2=4:6
                mat(k1,k2) = -ddirs(k2-3,k1,n-3);
            end
        end
        for k1=4:7
            for k2=4:6
                mat(k1,k2) = dodotd2s(a,k2-3,n-3,k1-3);
            end
        end
        for k1=4:7
            for k2=1:3
                mat(k1,k2) = d2thdq2(a,k2,k1-3,n-3);
            end
        end
        mat1(:,n)=mat*vec;
    end
    for k1=1:6
        for k2=1:6
            mat3(k1,k2)=mmats(a,k1,k2);
        end
    end
    for k1=1:6
        for k2=1:7
            mat4(k1,k2) = dzdoq(k1,k2,a);
            mat2(k2,k1) = dzdoq(k1,k2,a);
        end
    end
    result = result+mat1+mat2*mat3*mat4;
    % Penalty weight part of formula
    qthis      =  q(:,a);
    vec        =  [0;0;0;qthis];
    mat5       =  eye(7);
    mat5(1,1)  =  0;
    mat5(2,2)  =  0;
    mat5(3,3)  =  0;
    mat5       =  4*penalty_weight*(qthis'*qthis-1)*mat5;
    mat5       =  mat5+8*penalty_weight*(vec*vec');
    result     =  result+mat5;
    for j=1:7
        for k=1:7
            nset    =  nset+1;
            I(nset) =  rowrange(j);
            J(nset) = colrange(k);
            X(nset) =  result(j,k);
        end
    end
end
hess = sparse(I,J,X,zlen,zlen);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


