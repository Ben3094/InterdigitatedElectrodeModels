% Finger width is equal to inter-finger width in this test
clear all
close all

substrateEpsillon = 3.225;
substrateThickness = 50e-6;

interlayerEpsillon = 3.225;
interlayerThickness = 100e-6;

minSensitiveLayerEpsillon = 5.020823;
maxSensitiveLayerEpsillon = 12;
deltaSensitiveLayerEpsillon = maxSensitiveLayerEpsillon - minSensitiveLayerEpsillon;

minSensitiveLayerThickness = 5e-6;
maxSensitiveLayerThickness = 50e-6;
deltaSensitiveLayerThickness = 1e-6;
sensitiveLayerThickness = minSensitiveLayerThickness:deltaSensitiveLayerThickness:maxSensitiveLayerThickness;
meanSensitiveLayerThickness = (maxSensitiveLayerThickness - minSensitiveLayerThickness)/2;

minFingers = 5;
maxFingers = 15;
fingers = minFingers:maxFingers;
meanFingers = (maxFingers - minFingers)/2;

minFingerLength = 10e-3;
maxFingerLength = 30e-3;
fingerLengthSteps = 2e-3;
fingerLengths = minFingerLength:fingerLengthSteps:maxFingerLength;
meanFingerLength = (maxFingerLength - minFingerLength) / 2;

minInterFingerWidth = 50e-6
maxInterFingerWidth = 300e-6
meanInterFingerWidth = (maxInterFingerWidth - minInterFingerWidth) / 2
interFingerWidthStep = (maxInterFingerWidth - minInterFingerWidth)/100;
interFingerWidths = minInterFingerWidth:interFingerWidthStep:maxInterFingerWidth;
fingerWidth = interFingerWidths;

%% Sensitive layer thickness sweep
minCapacitances = zeros(length(sensitiveLayerThickness), 1);
maxCapacitances = zeros(length(sensitiveLayerThickness), 1);
meanCapacitances = zeros(length(sensitiveLayerThickness), 1);
for i = 1:length(sensitiveLayerThickness)
    minCapacitances(i) = ParallelPartialCapacitanceModel(substrateEpsillon, substrateThickness, [interlayerEpsillon minSensitiveLayerEpsillon], [interlayerThickness interlayerThickness+sensitiveLayerThickness(i)], meanFingers, meanFingerLength, meanInterFingerWidth, meanInterFingerWidth);
    maxCapacitances(i) = ParallelPartialCapacitanceModel(substrateEpsillon, substrateThickness, [interlayerEpsillon maxSensitiveLayerEpsillon], [interlayerThickness interlayerThickness+sensitiveLayerThickness(i)], meanFingers, meanFingerLength, meanInterFingerWidth, meanInterFingerWidth);
    meanCapacitances(i) = mean([minCapacitances(i), maxCapacitances(i)]);
end
figure
hold on
title('Igreja and Dias PPC model')
xlabel('Sensitive layer thickness (m)');
ylabel('Capacitance (F)');
plot(sensitiveLayerThickness, minCapacitances);
plot(sensitiveLayerThickness, meanCapacitances, '--');
plot(sensitiveLayerThickness, maxCapacitances);
legend('Min', 'Mean', 'Max');

deltaCapacitances = zeros(length(minCapacitances), 1);
for i = 1:length(minCapacitances)
    deltaCapacitances(i) = maxCapacitances(i) - minCapacitances(i);
end
% figure
% plot(sensitiveLayerThickness, deltaCapacitances/deltaSensitiveLayerEpsillon);
% title('Igreja and Dias PPC model sensitivity')
% xlabel('Sensitive layer thickness (m)');
% ylabel('Sensibility (F/(F/m))');

sensitivities = zeros(length(minCapacitances)-1, 1);
for i = 1:length(minCapacitances)-1
    sensitivities(i) = (deltaCapacitances(i+1) - deltaCapacitances(i)) / (sensitiveLayerThickness(i+1) - sensitiveLayerThickness(i));
end
figure
plot(sensitiveLayerThickness(2:end), sensitivities);
title('Igreja and Dias PPC model sensitivity')
xlabel('Sensitive layer thickness (m)');
ylabel('Sensibility (F/m)');

%% Inter-finger width sweep
minCapacitances = zeros(length(interFingerWidths), 1);
maxCapacitances = zeros(length(interFingerWidths), 1);
meanCapacitances = zeros(length(interFingerWidths), 1);
for i = 1:length(interFingerWidths)
    minCapacitances(i) = ParallelPartialCapacitanceModel(substrateEpsillon, substrateThickness, [interlayerEpsillon minSensitiveLayerEpsillon], [interlayerThickness interlayerThickness+meanSensitiveLayerThickness], meanFingers, meanFingerLength, interFingerWidths(i), interFingerWidths(i));
    maxCapacitances(i) = ParallelPartialCapacitanceModel(substrateEpsillon, substrateThickness, [interlayerEpsillon maxSensitiveLayerEpsillon], [interlayerThickness interlayerThickness+meanSensitiveLayerThickness], meanFingers, meanFingerLength, interFingerWidths(i), interFingerWidths(i));
    meanCapacitances(i) = mean([minCapacitances(i), maxCapacitances(i)]);
end
figure
hold on
title('Igreja and Dias PPC model')
xlabel('Inter-finger widths (m)');
ylabel('Capacitance (F)');
plot(interFingerWidths, minCapacitances);
plot(interFingerWidths, meanCapacitances, '--');
plot(interFingerWidths, maxCapacitances);
legend('Min', 'Mean', 'Max');

deltaCapacitances = zeros(length(minCapacitances), 1);
for i = 1:length(minCapacitances)
    deltaCapacitances(i) = maxCapacitances(i) - minCapacitances(i);
end
% figure
% plot(interFingerWidths, deltaCapacitances/deltaSensitiveLayerEpsillon);
% title('Igreja and Dias PPC model sensitivity')
% xlabel('Inter-finger widths (m)');
% ylabel('Sensibility (F/(F/m))');

sensitivities = zeros(length(minCapacitances)-1, 1);
for i = 1:length(minCapacitances)-1
    sensitivities(i) = (deltaCapacitances(i+1) - deltaCapacitances(i)) / (interFingerWidths(i+1) - interFingerWidths(i));
end
figure
plot(interFingerWidths(2:end), sensitivities);
title('Igreja and Dias PPC model sensitivity')
xlabel('Inter-finger widths (m)');
ylabel('Sensibility (F/m)');

