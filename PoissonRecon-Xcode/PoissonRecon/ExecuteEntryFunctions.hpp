//
//  ExecuteEntryFunctions.hpp
//  PoissonRecon
//
//  Created by Aaron Thompson on 2/13/19.
//  Copyright © 2019 Standard Cyborg. All rights reserved.
//

#ifndef ExecuteEntryFunctions_hpp
#define ExecuteEntryFunctions_hpp

extern void PoissonReconExecute(const char *inputFilePath, const char *outputFilePath);

/** Returns 0 on success, nonzero on error */
extern int SurfaceTrimmerExecute(const char* inputFilePath, const char* outputFilePath);

#endif /* ExecuteEntryFunctions_hpp */
