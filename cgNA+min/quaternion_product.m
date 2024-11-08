function [pq] = quaternion_product(p,q)
pq = zeros(4,1);
% p,q, pq are (4,1) vector
p1 = p(1);
p2 = p(2);
p3 = p(3);
p4 = p(4);

q1 = q(1);
q2 = q(2);
q3 = q(3);
q4 = q(4);

pq(1) = p1*q4 + p2*q3 - p3*q2 + p4*q1;
pq(2) = -p1*q3 + p2*q4 + p3*q1 + p4*q2;
pq(3) =  p1*q2 - p2*q1 + p3*q4 + p4*q3;
pq(4) = -p1*q1 - p2*q2 - p3*q3 + p4*q4;
end