function K = K(index_1,index_2, total_thickness_au, kane_energy_au, Eg_GaAs_au)
K = (pi/137)^2*kane_energy_au./(real(index_1.*index_2).*total_thickness_au*Eg_GaAs_au^4);
end

