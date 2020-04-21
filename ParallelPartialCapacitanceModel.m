function [C] = ParallelPartialCapacitanceModel(firstHalfPlaneEpsilons, firstHalfPlaneHeights, secondHalfPlaneEpsilons, secondHalfPlaneHeights, fingers, fingerLength, interFingerWidth, fingerWidth)
%ParallelPartialCapacitanceModel Summary of this function goes here
%   According to R. Igreja and C. J. Dias, ?Extension to the analytical model of the interdigital electrodes capacitance for a multi-layered structure,? Sensors and Actuators A: Physical, vol. 172, no. 2, pp. 392?399, Dec. 2011, doi: 10.1016/j.sna.2011.09.033.

    %% CONSTANTS
    voidPermittivity = 8.85418782e-12;
    
    %% VARIABLES CHECKING
        % First half plane
%         assert(any(size(firstHalfPlaneEpsilons) ~= 1), 'firstHalfPlaneEpsilons need to be unidimentional');
%         assert(any(size(firstHalfPlaneHeights) ~= 1), 'firstHalfPlaneHeights need to be unidimentional');
        assert(length(firstHalfPlaneEpsilons) == length(firstHalfPlaneHeights), 'firstHalfPlaneEpsilons must be the same size as firstHalfPlaneHeights');

            % Check if epsillon increase monotonicaly
            for i = 2:length(firstHalfPlaneEpsilons)
                assert(not(firstHalfPlaneEpsilons(i-1) > firstHalfPlaneEpsilons(i)), 'Parallel Partial Capacitance model can be applied on monotonal decreasing layer only');
            end
                
        % Second half plane
%         assert(any(size(secondHalfPlaneEpsilons) ~= 1), 'secondHalfPlaneEpsilons need to be unidimentional');
%         assert(any(size(secondHalfPlaneHeights) ~= 1), 'secondHalfPlaneHeights need to be unidimentional');
        assert(length(secondHalfPlaneEpsilons) == length(secondHalfPlaneHeights), 'secondHalfPlaneEpsilons must be the same size as secondHalfPlaneHeights');

            % Check if epsillon increase monotonicaly
            for i = 2:length(secondHalfPlaneEpsilons)
                assert(not(secondHalfPlaneEpsilons(i-1) > secondHalfPlaneEpsilons(i)), 'Parallel Partial Capacitance model can be applied on monotonal decreasing layer only');
            end
            
    %% NORMALIZATION
    fingers = round(abs(fingers));
    firstHalfPlaneHeights= abs(firstHalfPlaneHeights);
    secondHalfPlaneHeights = abs(secondHalfPlaneHeights);
    fingerLength = abs(fingerLength);
    interFingerWidth = abs(interFingerWidth);
    fingerWidth = abs(fingerWidth);
    
    %% ADIMENTIONAL COEFFICIENT COMPUTATION
    lambda = 2*(fingerWidth+interFingerWidth); % independent on layer height
    nu = 2*fingerWidth/lambda; % metalization ratio % independent on layer height
    
    % Compute interior electrodes in void
    kIInfinite = sin(pi*nu/2);
    infiniteInteriorAirLayerCapacitance = voidPermittivity*fingerLength*ellipke(kIInfinite)/ellipke(primeComputation(kIInfinite));
    
    % Compute exterior electrodes in void
    kEInfinite = 2*sqrt(nu)/(1+nu);
    infiniteExteriorAirLayerCapacitance = voidPermittivity*fingerLength*ellipke(kEInfinite)/ellipke(primeComputation(kEInfinite));
   
    C = partialCapacitanceComputation(firstHalfPlaneEpsilons, firstHalfPlaneHeights) + partialCapacitanceComputation(secondHalfPlaneEpsilons, secondHalfPlaneHeights);
    
    function [partialCapacitance] = partialCapacitanceComputation(halfPlaneEpsilons, halfPlaneHeights)
        halfPlaneEpsilons = [halfPlaneEpsilons, 1]; % Add air layer
        
        % Compute coefficients dependent on correponding layer height
        r = halfPlaneHeights/lambda;
        q = exp(-4*pi*r);
        k = (theta(2, 0, q)/theta(3, 0, q))^2;
        
        % Compute internal capacitance
        [t2, ~, ~] = ellipj(ellipke(k)*nu, k);
        t4 = 1/k;
        kI = t2*sqrt((t4^2 - 1)/(t4^2 - t2^2));
        partialInteriorCapacitance = infiniteInteriorAirLayerCapacitance;
         for i = 1:(length(halfPlaneEpsilons) - 1)
            partialInteriorCapacitance = partialInteriorCapacitance + (halfPlaneEpsilons(i)-halfPlaneEpsilons(i+1))*voidPermittivity*fingerLength*ellipke(kI(i))/ellipke(primeComputation(kI(i)));
        end
        
        % Compute internal capacitance
        t3 = cosh(pi*(1-nu)/(8*r));
        t4 = cosh(pi*(nu+1)/(8*r));
        kE = sqrt((t4^2 - t3^2)/(t4^2 - 1)) / t3;
        partialExteriorCapacitance = infiniteExteriorAirLayerCapacitance;
        for i = 1:(length(halfPlaneEpsilons) - 1)
            partialExteriorCapacitance = partialExteriorCapacitance + (halfPlaneEpsilons(i)-halfPlaneEpsilons(i+1))*voidPermittivity*fingerLength*ellipke(kE(i))/ellipke(primeComputation(kE(i)));
        end
        
    
    partialCapacitance = (fingers-3)*partialInteriorCapacitance/2 + 2*partialInteriorCapacitance*partialExteriorCapacitance/(partialInteriorCapacitance+partialExteriorCapacitance);
    end
    
    function [coefficient] = primeComputation(coefficient)
        coefficient = sqrt(1-coefficient^2);
    end
end