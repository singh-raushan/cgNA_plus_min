function [inters,intras] = compute_internal_variables(zvec)
global whats
global q4_at_1
% We have 18*nbp intra (base-pair level coordinates) unknowns and 7*(nbp-2) unknowns (o,q),
%    with the first (o,q)=(0,e4) and the last (0,\pm e4) depending on initial guess
% Ordering of zvec is intra1, then intra2,oq2,intra3,oq3,...oqn-2,intran-1,oqn-1 then intran, oqn
s = size(zvec);
zlen = s(1);
nbp = (zlen+7)/25;

% Read (o,q) values from zvec -- and assign values at ends
q = zeros(4,nbp+1);
o = zeros(3,nbp+1);
for i=2:nbp
    o(:,i)=zvec(36+25*(i-2)+1:36+25*(i-2)+3);
    q(:,i)=zvec(36+25*(i-2)+4:36+25*(i-2)+7);
end
q(4,1)=1;
q(4,nbp+1)=q4_at_1;

% Read intra values from zvec
intras = zeros(18,nbp);
intradiffs = zeros(18,nbp);
for i=2:nbp-1
    intradiffs(:,i)=zvec(18+25*(i-2)+1:18+25*(i-2)+18);
end
intradiffs(:,1)=zvec(1:18);
intradiffs(:,nbp)=zvec(zlen-24:zlen-7);
for i=1:nbp
    intras(:,i)=intradiffs(:,i)+whats(24*(i-1)+1:24*(i-1)+18);
end

% Set up key constant matrices
id = eye(4);
b = zeros(3,4,4);
b(1,1,4) = 1; b(1,2,3) = 1; b(1,3,2)=-1; b(1,4,1) = -1;
b(2,1,3) = -1; b(2,2,4) = 1; b(2,3,1) = 1; b(2,4,2) = -1;
b(3,1,2) = 1; b(3,2,1) = -1; b(3,3,4) = 1; b(3,4,3) = -1;

% Store values of B1q,B2q,B3q in triply-indexed array bq
bq = zeros(3,4,nbp+1);
for i=1:3
    bi = zeros(4,4);
    for j1 = 1:4
        for j2 = 1:4
            bi(j1,j2)=b(i,j1,j2);
        end
    end
    for k = 1:nbp+1
        dum = bi*q(:,k);
        for j1=1:4
            bq(i,j1,k)=dum(j1);
        end
    end
end

% Compute the inters
cay = zeros(3,nbp);
tr = zeros(3,nbp);
inters = zeros(6,nbp);
dfac = zeros(1,nbp);
for i=1:nbp
    dfac(i) = q(:,i+1)'*q(:,i);
    for k=1:3
        bkq = zeros(4,1);
        for j1 = 1:4
            bkq(j1)=bq(k,j1,i);
        end
        cay(k,i) = 10/dfac(i)*q(:,i+1)'*bkq;
        inters(k,i)=cay(k,i);
    end
    do = o(:,i+1)-o(:,i);
    if dot(q(:,i+1),q(:,i))<0
        sign=-1;
    else
        sign=1;
    end
    dirs = compute_ds(q(:,i+1)+sign*q(:,i));
    for k=1:3
        tr(k,i)=do'*dirs(:,k);
        inters(k+3,i)=tr(k,i);
    end
end

end

