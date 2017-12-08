`.cholInvDiag_DLLInfo` <- dyn.load('\\\\SSRI-NAS-FE01.oit.duke.edu\\SSRI\\OPM\\Users\\Current\\tjb48\\Analysis\\FixedEffectsModel\\LargeFixedEffectsModel\\CholeskyDecomposition\\cholInvDiag.dll')

cholInvDiag <- Rcpp:::sourceCppFunction(function(R0) {}, FALSE, `.cholInvDiag_DLLInfo`, 'sourceCpp_1_cholInvDiag')

rm(`.cholInvDiag_DLLInfo`)
