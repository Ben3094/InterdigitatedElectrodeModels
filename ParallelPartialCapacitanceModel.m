function C = ParallelPartialCapacitanceModel(upperHalfPlaneEpsilons, upperHalfPlaneHeights, lowerHalfPlaneEpsilons, lowerHalfPlaneHeights, fingers, fingerLength, interFingerWidth, fingerWidth, upperHalfPlaneSurrondingEpsilon, lowerHalfPlaneSurrondingEpsilon)
%ParallelPartialCapacitanceModel Compute parallel partial capacitance model for interdigitated electrodes with multiple material layers
% According to R. Igreja and C. J. Dias, ?Extension to the analytical model of the interdigital electrodes capacitance for a multi-layered structure,? Sensors and Actuators A: Physical, vol. 172, no. 2, pp. 392?399, Dec. 2011, doi: 10.1016/j.sna.2011.09.033.
% Uses the package of Moiseev Igor (2020). Elliptic Integrals and Functions (https://www.mathworks.com/matlabcentral/fileexchange/8805-elliptic-integrals-and-functions), MATLAB Central File Exchange. Retrieved July 1, 2020.

    %% CONSTANTS
    voidPermittivity = 8.85418782e-12;
    
    %% VARIABLES CHECKING
    if nargin < 9
        upperHalfPlaneSurrondingEpsilon = 1;
        lowerHalfPlaneSurrondingEpsilon = 1;
    end
    
    % Check if half planes agree
    assert(length(upperHalfPlaneEpsilons) == length(upperHalfPlaneHeights), 'upperHalfPlaneEpsilons must be the same size as upperHalfPlaneHeights');
    assert(length(lowerHalfPlaneEpsilons) == length(lowerHalfPlaneHeights), 'lowerHalfPlaneEpsilons must be the same size as lowerHalfPlaneHeights');

    % Check if epsillon increase monotonicaly
    for i = 2:length(upperHalfPlaneEpsilons)
        assert(not(upperHalfPlaneEpsilons(i-1) >= upperHalfPlaneEpsilons(i)), 'Parallel Partial Capacitance model can be applied on monotonal decreasing layer only');
        assert(not(lowerHalfPlaneEpsilons(i-1) >= lowerHalfPlaneEpsilons(i)), 'Parallel Partial Capacitance model can be applied on monotonal decreasing layer only');
    end
            
    %% NORMALIZATION
    fingers = round(abs(fingers));
    upperHalfPlaneHeights= abs(upperHalfPlaneHeights);
    lowerHalfPlaneHeights = abs(lowerHalfPlaneHeights);
    fingerLength = abs(fingerLength);
    interFingerWidth = abs(interFingerWidth);
    fingerWidth = abs(fingerWidth);
    
    %% ADIMENTIONAL COEFFICIENT COMPUTATION
    lambda = 2*(fingerWidth+interFingerWidth); % independent on layer height
    nu = 2*fingerWidth/lambda; % metalization ratio % independent on layer height
    
    %% COMPUTE CAPACITANCE
    % Compute interior electrodes in void
    kIInfinite = sin(pi*nu/2);
    infiniteInteriorVoidLayerCapacitance = voidPermittivity*fingerLength*ellipke(kIInfinite)/ellipke(primeComputation(kIInfinite));
    
    % Compute exterior electrodes in void
    kEInfinite = 2*sqrt(nu)/(1+nu);
    infiniteExteriorVoidLayerCapacitance = voidPermittivity*fingerLength*ellipke(kEInfinite)/ellipke(primeComputation(kEInfinite));
   
    C = partialCapacitanceComputation(upperHalfPlaneEpsilons, upperHalfPlaneHeights, upperHalfPlaneSurrondingEpsilon) + partialCapacitanceComputation(lowerHalfPlaneEpsilons, lowerHalfPlaneHeights, lowerHalfPlaneSurrondingEpsilon);
    
    function [partialCapacitance] = partialCapacitanceComputation(halfPlaneEpsilons, halfPlaneHeights, halfPlaneSurrondingEpsilon)
        planeNumber = length(halfPlaneHeights);
        halfPlaneEpsilons = [halfPlaneEpsilons halfPlaneSurrondingEpsilon]; % Add surronding layer
        
        kI = zeros(planeNumber, 1);
        kE = zeros(planeNumber, 1);
        for i = 1:planeNumber
            % Compute coefficients dependent on correponding layer height
            r = halfPlaneHeights(i)/lambda;
        
            % Compute internal capacitance
            q = exp(-4*pi*r);
            k = (theta(2, 0, q)/theta(3, 0, q))^2;
            t4I = 1/k;
            [t2, ~, ~] = ellipj(ellipke(k)*nu, k^2);
            % kI(i) = sqrt((t4I^2 - 1)/(t4I^2 - t2^2)); % According to current article
            kI(i) = t2 * sqrt((t4I^2 - 1)/(t4I^2 - t2^2)); % According to R. Igreja and C. J. Dias, ?Analytical evaluation of the interdigital electrodes capacitance for a multi-layered structure,? Sensors and Actuators A: Physical, vol. 112, no. 2, pp. 291?301, May 2004, doi: 10.1016/j.sna.2004.01.040.
        
            % Compute external capacitance
            t4E = cosh(pi*(nu+1)/(8*r));
            t3 = cosh(pi*(1-nu)/(8*r));
            kE(i) = sqrt((t4E^2 - t3^2)/(t4E^2 - 1)) / t3;
        end
           
        % Compute internal capacitance 
        partialInteriorCapacitance = halfPlaneSurrondingEpsilon*infiniteInteriorVoidLayerCapacitance;
        partialExteriorCapacitance = halfPlaneSurrondingEpsilon*infiniteExteriorVoidLayerCapacitance;
        for i = 1:planeNumber
            partialInteriorCapacitance = partialInteriorCapacitance + (halfPlaneEpsilons(i)-halfPlaneEpsilons(i+1))*voidPermittivity*fingerLength*ellipke(kI(i))/ellipke(primeComputation(kI(i)));
            partialExteriorCapacitance = partialExteriorCapacitance + (halfPlaneEpsilons(i)-halfPlaneEpsilons(i+1))*voidPermittivity*fingerLength*ellipke(kE(i))/ellipke(primeComputation(kE(i)));
        end
        
        partialCapacitance = (fingers-3)*partialInteriorCapacitance/2 + 2*partialInteriorCapacitance*partialExteriorCapacitance/(partialInteriorCapacitance+partialExteriorCapacitance);
    end
    
    function [primeModulus] = primeComputation(modulus)
        primeModulus = sqrt(1-modulus^2);
    end
end
