function Stelsel = Stelsel(t, y, u, delta_L, iter_L, size_L, V_in, delta_Z, C_0, C_sat, MASS_DENS, K_V, K_g, K_b, g, b)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initiation


% Initiate all possible equations. This optimizes the program.
Stelsel = zeros(iter_L + 3, 1);

% L_bins are all the bins in which a size of crystal can occur. Visually as
% a histogram, this would be the width of the bin.
L_bins = [0:delta_L: size_L]';


% Growth [m_crystal/ s] and Birth [(# of crystals)/ (m^3_space* s)]
G = K_g* (y(iter_L + 2) - C_sat)^(g);
B = K_b* (y(iter_L + 2) - C_sat)^(b);

% Test value 
D_in = 1.6*10^-3; %[m_space]
U_thermal = 500; % [J/(m^2_space *s *K)]
a = 4/D_in; %Area per unit of volume [1/m_space]
MASS_DENS_SOL = 1000; %[kg/m^3_space]
C_p = 4130; %[J/(kg *K)]
T_Wall = 20 + 273; %[K]

% The Flag is only 1 in the first iteration. In a future step the Flag is
% multiplied by birth, so this only comes into play in the first iteration. 
% Otherwise birth would always dominate over growth.
if C_0 == y(iter_L + 2)
    Flag = 1;
else
    Flag = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations


% This is the boundary condition and how it changes. Read: dn_0/dz. After 1
% iteration the birth term will be eliminated. It is missing a negative
% term at the moment, thus it will be a constant after the first iteration.
Stelsel(1) = B*Flag/ (u* delta_L) - G*y(1)/(u* delta_L);


% This for loop produces the system of equations
for i = 2:(iter_L + 1)
    % This is the Population Balance Equation (PBE). It is discretized backwards,
    % because we otherwise can't properly use the boundary condition. 
    % All equations from the 2 up and including the iter_L + 1'st line of
    % Stelsel are discretized PBE equations. This describes the whole tube.
    Stelsel(i) = -(G/ (u* delta_L))* (y(i) - y(i-1));

end


% The iter_L + 2'nd equation is the Mass Balance equation (MB). This is
% a discretized integral, which needs all values of the discretized
% population balance model except for the boundary. 
% As a result it is intuitive to build this equation every iteration
% instead of making it at once.
Stelsel(iter_L + 2) = - 3* K_V* MASS_DENS * G/u * (dot((L_bins.^2) , y(1:iter_L +1)) * delta_L);



Stelsel(iter_L + 3) = - U_thermal*a/(MASS_DENS_SOL* C_p* u) *(y(iter_L + 3) - T_Wall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end