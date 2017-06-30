`.choleskyDecomp_DLLInfo` <- dyn.load('\\\\SSRI-NAS-FE01.oit.duke.edu\\SSRI\\OPM\\Users\\Current\\tjb48\\Analysis\\FixedEffectsModel\\LargeFixedEffectsModel\\CholeskyDecomposition\\choleskyDecomp.dll')

choleskyDecomp <- Rcpp:::sourceCppFunction(function(X) {}, FALSE, `.choleskyDecomp_DLLInfo`, 'sourceCpp_3_choleskyDecomp')

rm(`.choleskyDecomp_DLLInfo`)
