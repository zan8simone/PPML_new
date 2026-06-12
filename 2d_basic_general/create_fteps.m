function fteps = create_fteps(lattice,epsmask)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is free software distributed under the BSD licence (see the 
%  containing folder).
% However, shall the results obtained through this code be included 
%  in an academic publication, we kindly ask you to cite the source 
%  website.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates the Fourier transform of the epsilon "mask" 
%
% 2-d structure, arbitrary lattice
%
% Simone Zanotto, 2026
%
% ========= INPUT
%
% lattice: stucture as created by create_lattice(...)
% epsmask is a structure containing fields  epsmask.xx
%                                           epsmask.xy
%                                            etc.
%                        each of which is a [lattice.N1 x lattice.N2 complex matrix] 
%                        containing the values of permittivity 
%                        at the real space grid points 
%
%  ======== OUTPUT
%
%  fteps.xx, fteps.xy etc: [lattice.N x lattice.N complex matrix] containing the FTs 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fnv = {'xx','xy','yx','yy','z'}; % list of names from tensor components of permittivity

[N1,N2] = size(epsmask.(fnv{1}));
if (N1 ~= lattice.N1) || (N2 ~= lattice.N2)
    error('The fields in epsmask must have the same size as direct lattice grid');
end

for ifn = 1:5

    [N1t,N2t] = size(epsmask.(fnv{ifn}));
    if (N1t ~= N1) || (N2t ~= N2)
        error('The fields in epsmask must have the same size');
    end

fttmp = fftshift( fft2( epsmask.(fnv{ifn}) ) )/N1/N2;

ft = zeros(size(lattice.G,1));

for I = 1:size(lattice.G,1)
for J = 1:size(lattice.G,1)
    ft(I,J) = fttmp( lattice.mn(I,1) - lattice.mn(J,1) + lattice.N1/2 + 1, ...
                     lattice.mn(I,2) - lattice.mn(J,2) + lattice.N2/2 + 1); 
end
end

fteps.(fnv{ifn}) = ft;


end
