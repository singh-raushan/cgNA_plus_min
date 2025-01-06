function Lk = compute_lk(rw, rc)
%--------------------------------------------------------------------------
% cgDNA function: [Lk] = compute_lk(rw, rc)
%--------------------------------------------------------------------------
% This function computes the linking number of two piecewise linear closed
% curves
%
% Input:
%
% rw                points of one of the curves (designed Watson) as [3, N]
%
% rc                points of the other curve (designed Crick) in the same form
%
% Output:
%
%   Lk              The linking number of the two curves
%
% Adaptation of method 1a for approximating writhe
% of Klenin, Langowski 2000
% It uses the two curves and so should only be done for closed solutions
% (at least one of the curves needs to be closed)

    % Indices M1 N1, M2 N2 of points r_{M1 N1} and r_{M2 N2}
    % To compute the normal vectors in equ. (15)
    nIndices = [ 1, 3, 1, 4;
                 1, 4, 2, 4;
                 2, 4, 2, 3;
                 2, 3, 1, 3 ];

    [dimRw, numPtsRw] = size(rw);
    [dimRc, numPtsRc] = size(rc);
    if (dimRw ~= 3 || dimRc ~= 3 || numPtsRw < 3 || numPtsRc < 3)
        error('Incorrect input dimensions')
    end

    omega = 0.0;

    % The points r1, r2, r3, r4
    r = zeros(3, 4);
    % The normals n1, n2, n3, n4
    n = zeros(3, 4);

    % Run the double summation analogous to equ. (13) of K&L2000
    % but here there are two curves and both sums need to be complete
    for i = 2 : numPtsRw
        % Get the points r1 and r2 (Watson)
        r(:, 1) = rw(:, i - 1);
        r(:, 2) = rw(:, i);
        % Run through all segments of the Crick curve
        for j = 2: numPtsRc
            % Get the points r3 and r4 (Crick)
            r(:, 3) = rc(:, j - 1);
            r(:, 4) = rc(:, j);
            % Compute the normals: equ. (15)
            for k = 1:4
                tmp1 = r(:, nIndices(k, 2)) - r(:, nIndices(k, 1));
                tmp2 = r(:, nIndices(k, 4)) - r(:, nIndices(k, 3));

                % Compute the cross product
                prod = cross(tmp1, tmp2);
                n(:, k) = prod / norm(prod);
            end

            % Equ. (16a)
            omegaStar = 0;
            for k = 1:4
                dotP = n(:, k)' * n(:, mod(k, 4) +1);
                % Take care not to exceed (-1.0, 1.0) even by a small fraction
                % as asin will return a complex number
                omegaStar = omegaStar + asin(max(min(dotP, 1), -1));
            end

            % Equ. (16b)
            r34 = r(:, 4) - r(:, 3);
            r12 = r(:, 2) - r(:, 1);
            r13 = r(:, 3) - r(:, 1);
            % (r34 x r12) . r13
			omega = omega + omegaStar * sign(cross(r34, r12)' * r13);
        end
    end

    % Analogue to Equ (13)
    Lk = 0.25 / pi * omega;
end

