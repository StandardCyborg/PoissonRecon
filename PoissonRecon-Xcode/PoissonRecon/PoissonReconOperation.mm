//
//  PoissonReconOperation.mm
//  PoissonRecon
//
//  Created by Aaron Thompson on 4/26/18.
//  Copyright © 2018 Standard Cyborg. All rights reserved.
//

#import <PoissonRecon/PoissonReconOperation.h>
#import "Parameters.hpp"
#import "ExecuteEntryFunctions.hpp"

using namespace std;

static float fclampf(float value, float min, float max) {
    return fmaxf(fminf(value, max), min);
}

static float remapAndClamp(float value, float originalMin, float originalMax, float newMin, float newMax) {
    float shiftedOriginal = value - originalMin;
    float scaled = shiftedOriginal * (newMax - newMin) / (originalMax - originalMin);
    float shiftedNew = scaled + newMin;
    
    return fclampf(shiftedNew, newMin, newMax);
}

@implementation PoissonReconOperation {
    NSString *_inputFilePath;
    NSString *_outputFilePath;
}

- (instancetype)initWithInputFilePath:(NSString *)inputPath
                       outputFilePath:(NSString *)outputPath
{
    self = [super init];
    if (self) {
        _inputFilePath = inputPath;
        _outputFilePath = outputPath;
        
        _resolution = 5;
        _smoothness = 2;
    }
    return self;
}

- (void)main
{
    const char *inputPath = [_inputFilePath UTF8String];
    NSString *tempPoissonOutputPathString = [NSTemporaryDirectory() stringByAppendingFormat:@"/poisson-%@.ply", [[NSUUID UUID] UUIDString]];
    const char *poissonOutputPath = [tempPoissonOutputPathString UTF8String];
    const char *surfaceTrimmerOutputPath = [_outputFilePath UTF8String];
    
    PoissonReconParameters poissonParams;
    poissonParams.Depth = remapAndClamp(_resolution, 1, 10, 4, 14);
    poissonParams.SamplesPerNode = remapAndClamp(_smoothness, 1, 10, 1, 15);
    
    SurfaceTrimmerParameters surfaceTrimmerParams;
    // The defaults are all fine for this one, so not exposing any knobs
    
    PoissonReconExecute(inputPath, poissonOutputPath, poissonParams);
    SurfaceTrimmerExecute(poissonOutputPath, surfaceTrimmerOutputPath, surfaceTrimmerParams);
    
    [[NSFileManager defaultManager] removeItemAtPath:tempPoissonOutputPathString error:NULL];
}

@end
