function [shapes, stiff] = constructSeqParmsPeriodic(seq, ps)

% A modification of the original cgDNA constructSeqParms -> cgDNA+ periodic
%
% Constructs periodic cgDNA+ ground-state coordinate vector and stiffness
% matrix for a given sequence, using the specified parameter set
%
% This is used for minicircle computatins

%%%% dimer = params.dimer;
%%%% base  = params.base;
stiff_int = ps.stiff_int;
sigma_int = ps.sigma_int;


% Initialize variables
nbp = numel(seq);


nv    = 24;
nover = 18;
N     = nv * nbp;
%%%% nz    = (nv^2 + 2*nv*nover) * (nbp-1) + nover^2; % this is written by Glowacki but it wrong it corresponds to non periodic
nz    = (nv^2 + 2*nv*nover) * (nbp);  % written by Raushan for periodic
stiff = spalloc(N,N,nz);
sigma = zeros(N,1);


% Assemble stiffness matrix 'stiff' and sigma vector 'sigma'
for i = 1:nbp
    k = (i-1)*24 + 1;
    if(i < nbp)
        stiff_block_dimer  = stiff_int.(seq(i:i+1));
        sigma_block_dimer  = sigma_int.(seq(i:i+1));
        
        stiff(k:k+41,k:k+41) = stiff(k:k+41,k:k+41) + stiff_block_dimer;   % dimer(fsi(dimer,seq(i:i+1))).b18;
        sigma(k:k+41)        = sigma(k:k+41) + sigma_block_dimer;          % dimer(fsi(dimer,seq(i:i+1))).c18;
    end
    %   stiff(k:k+5,k:k+5) = stiff(k:k+5,k:k+5) + base(fsi(base,seq(i))).b;
    %   sigma(k:k+5)       = sigma(k:k+5) + base(fsi(base,seq(i))).c;
end

%%%% lastStep = [seq(end) seq(1)];
stiff_XNX1_dimer = stiff_int.(seq(end:-(nbp-1):1));
sigma_XNX1_dimer = sigma_int.(seq(end:-(nbp-1):1));

% 24 of the 42 parameters go directly at the "end"
stiff(k:k+23,k:k+23) = stiff(k:k+23,k:k+23) + stiff_XNX1_dimer(1:24,1:24); % dimer(fsi(dimer,lastStep)).b18(1:12,1:12);
sigma(k:k+23)        = sigma(k:k+23) + sigma_XNX1_dimer(1:24);             % dimer(fsi(dimer,lastStep)).c18(1:12);
% The periodicity implies extra overlap of the first intras...
stiff(1:18,1:18)       = stiff(1:18,1:18) + stiff_XNX1_dimer(25:42,25:42);     % dimer(fsi(dimer,lastStep)).b18(13:18,13:18);
sigma(1:18)           = sigma(1:18) + sigma_XNX1_dimer(25:42);               % dimer(fsi(dimer,lastStep)).c18(13:18);
% ... and couplings between begining and end
stiff(1:18,k:k+23)    = stiff(1:18,k:k+23) + stiff_XNX1_dimer(25:42,1:24);   % dimer(fsi(dimer,lastStep)).b18(13:18,1:12);
stiff(k:k+23,1:18)    = stiff(k:k+23,1:18) + stiff_XNX1_dimer(1:24,25:42);   % dimer(fsi(dimer,lastStep)).b18(1:12,13:18);

% Compute ground-state coord vector via matrix inversion
shapes = stiff\sigma;

end

%--------------------------------------------------------
% function i = fsi(struc, s)
%
% n = size(struc);
%
% for j = 1:n(2)
%     if(strcmpi(struc(j).S,s))
%         i=j;
%         return;
%     end;
% end
%
% end



