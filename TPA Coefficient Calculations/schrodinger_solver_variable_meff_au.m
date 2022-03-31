function [wavefunctions, energies, dz] = schrodinger_solver_variable_meff_au(V_au, z_au, effective_mass_function_au, E_increment_au)
% Solves the Schrodinger equation for a given potential distribution using
% the shoothing method (see book  Quantum Wells, Wires and Dots by Harrison
% for a detailed description).
%   V: potential distribution, [meV]
%   z: z position, [Angstroms]. V should be defined relative to these coordinates
%   Meff: The effective mass at every given z [kg]
%   Eincrement: [meV] This parameter determines the step size in the
%       shooting method search for solutions. With many wells or complex
%       structures it is possible that two energy levels fit within one
%       Eincrement. In this case, the solutions will flip divergence twice,
%       leading the solver to miss both of them. It is very important to count
%       the states and make sure they match what you'd expect (i.e. 4 quantum
%       wells should have 4 states per n value except possibly at highest n
%       value). Adjust Eincrement accordingly.

    % Define the energy range to search for solutions. By default we will
    % search all values between the lowest potential and highest potential
    % ('bound' states, although they may not be truly bound when a bias is
    % applied)
    E = min(V_au):E_increment_au:max(V_au);

    % Find step size for the position. This assumes a constant step size
    dz = z_au(2) - z_au(1);
    
    % For the numerical solution we must define effective mass at midpoint
    % values between z steps. We define this as the average between two
    % nearest values    
    
    % div will tell us which direction the solution diverges for a given
    % energy guess
    div = zeros(size(E));

% This loop goes through for each energy guess, stepping by Eincrement, and
% calculates which way the wavefunction diverges. Later, we will know where
% to look more closely to find the exact solutions.

for jnd = 1:length(E)
    
    effective_mass = effective_mass_function_au(E(jnd));
    effective_mass_plus_half = 1/2*(effective_mass + [effective_mass(2:end) effective_mass(end)]);
    effective_mass_minus_half = 1/2*(effective_mass + [effective_mass(1) effective_mass(1:end-1)]);
    
    % X = wave function. Initialize the function and set X(1) to 0 and X(2)
    % to 1. Somehow these boundary conditions just work for this specific
    % solver, see the book for more details.
    X = zeros(size(z_au));
    X(1) = 0;
    X(2) = 1;
    
    % Build up the wavefunction X(z) index by index. This loop starts at 2.
    % Since X(3) = f(X(2), X(1)), we can solve for X(3) based on our
    % boundary conditions. Continue for the rest of the wavefunction.
    for ind = 2:length(z_au)-1
        X(ind + 1) = effective_mass_plus_half(ind)*((2*dz^2*(V_au(ind) - E(jnd)) + 1/effective_mass_plus_half(ind) + 1/effective_mass_minus_half(ind))*X(ind) - X(ind - 1)/effective_mass_minus_half(ind));
    end
    
    % We have now fully solved the schrodinger equation at the guessed
    % energy. Check which direction it diverges (+ or - is all we need)
    div(jnd) = sign(X(end));
    
end

% Using the vector of divergence values (which will be a vector filled with
% -1 and +1) we find the indices where we flip from -1 to +1. Basically,
% this means the difference between two adjacent values is nonzero.
crossing_indices = find(abs(diff(div)) > 0);

% Find the energies at which these crossings occured. Remember each element
% of div corresponds to the divergence for the energy with the same index
% in the E vector.
crossings = E(crossing_indices);

% Now we check if it crossed from - to + or + to -.
crossing_signs = sign(div(crossing_indices - 1));

% By the number of crossings we know how many solutions we have (assuming
% we did not miss any by having Eincrement too large. Define a matrix where
% the rows are the wavefunctions for each solution. Same for energies, 1
% for each solution
wavefunctions = zeros(length(crossings), length(z_au));
energies = zeros(length(crossings), 1);

% Now go through each crossing value and look more closely for the actual
% solution. To get a good wavefunction it is very important to have the
% energy very accurate.
for lnd = 1:length(crossings)
    
    % Recall the crossing energy
    En = crossings(lnd);
    
    % decrease the step size by a factor of 2
    dE = E_increment_au/2;
    
    % Initialize wavefunction with boundary conditions described earlier
    X = zeros(size(z_au));
    X(1) = 0;
    X(2) = 1;
    
    % Keep track of how many attempts we have made at finding a solution
    attempts = 0;
    
    % Some special considerations need to be made for the first guess at a
    % solution
    first_run = 1;
    
    % abs(X(end))/abs(max(X)) >= 1e-5 checks to make sure the wavefunction
    % goes to zero at z(end). This solver is made only for bound states. We
    % also time out after 100 attempts (this can be adjusted for
    % hard-to-fnd solutions). The first_run is not a very elegant solution,
    % but the while loop conditions are not met otherwise.
    while (abs(X(end))/abs(max(X)) >= 1e-5 && (attempts < 100) || first_run)

        effective_mass = effective_mass_function_au(En);
        effective_mass_plus_half = 1/2*(effective_mass + [effective_mass(2:end) effective_mass(end)]);
        effective_mass_minus_half = 1/2*(effective_mass + [effective_mass(1) effective_mass(1:end-1)]);

        % Build wavefunction step by step
        for ind = 2:length(z_au)-1
            X(ind + 1) = effective_mass_plus_half(ind)*((2*dz^2*(V_au(ind) - En) + 1/effective_mass_plus_half(ind) + 1/effective_mass_minus_half(ind))*X(ind) - X(ind - 1)/effective_mass_minus_half(ind));
        end
    
        % Keep track of the energy we just guessed
        En_prev = En;
        
        % Now increment that energy by dE in the direction of the crossing.
        % If we go from - to + we need to increase the energy, + to - we
        % need to decrease (or the other way around, 
        En = En + dE*crossing_signs(lnd)*sign(X(end));
        
        % Cut the energy increments in half again (like a binary search)
        dE = dE/2;
        
        % increment the number of attempts
        attempts = attempts + 1;
        
        % We are no longer in the first run
        first_run = 0;
    end
    
    % Normalize the wavefunctions and store them in the output variables
    % wavefunctions and energies
    X = X/norm(X)/sqrt(dz);
    wavefunctions(lnd, :) = X;
    energies(lnd) = En_prev;
        
end


% If we have negative effective masses, make the output energies and
% wavefunctions negative.
if (mean(effective_mass) < 0)
    energies = flipud(energies);
    wavefunctions = flipud(wavefunctions);
end