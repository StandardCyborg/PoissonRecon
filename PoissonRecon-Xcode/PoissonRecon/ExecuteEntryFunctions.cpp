//
//  ExecuteEntryFunctions.cpp
//  PoissonRecon
//
//  Created by Aaron Thompson on 2/13/19.
//  Copyright © 2019 Standard Cyborg. All rights reserved.
//

#include "ExecuteEntryFunctions.hpp"
#include "PoissonReconExecute.hpp"
#include "SurfaceTrimmerExecute.hpp"

/**
 Why have these together instead of moving them into their respective *Execute source files?
 - There's a conflict between the OS-provided Point struct type and this library's Point,
   so you can't include both together.
 - If the compiler doesn't specialize both at the same time, it ends up with duplicate symbols
   when linking because it didn't know to name them differently (I think?)
 */

void PoissonReconExecute(const char *inputFilePath, const char *outputFilePath)
{
    _PoissonReconExecute<PointStreamColor<float>>(inputFilePath, outputFilePath, IsotropicUIntPack<3, FEMDegreeAndBType<1, BOUNDARY_NEUMANN>::Signature>());
}

int SurfaceTrimmerExecute(const char* inputFilePath, const char* outputFilePath)
{
    return _SurfaceTrimmerExecute<PointStreamNormal<float, 3>, PointStreamColor<float>>(inputFilePath, outputFilePath);
}
