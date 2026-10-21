function kTrue = getTrueCurvature(cleanFile)
D = load(cleanFile);
kTrue = D.trajectoryKappaTrue;
end