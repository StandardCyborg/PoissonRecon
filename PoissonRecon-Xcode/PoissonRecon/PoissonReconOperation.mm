//
//  PoissonReconOperation.mm
//  PoissonRecon
//
//  Created by Aaron Thompson on 4/26/18.
//  Copyright © 2018 Standard Cyborg. All rights reserved.
//

#import "PoissonReconOperation.h"
#import "ExecuteEntryFunctions.hpp"

using namespace std;

// extern int PoissonRecon_main(int argc, char *argv[]);
// extern int SurfaceTrimmer_main(int argc, const char *argv[]);

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
    }
    return self;
}

- (void)main
{
    // const char *workingPath = [[[NSBundle bundleForClass:[self class]] bundlePath] UTF8String];
    // const char *poissonOutputPath = [_outputFilePath UTF8String];
    const char *inputPath = [_inputFilePath UTF8String];
    NSString *poissonOutputPathString = [NSTemporaryDirectory() stringByAppendingPathComponent:@"poisson.ply"];
    const char *poissonOutputPath = [poissonOutputPathString UTF8String];
    const char *surfaceTrimmerOutputPath = [_outputFilePath UTF8String];
    
    PoissonReconExecute(inputPath, poissonOutputPath);
    SurfaceTrimmerExecute(poissonOutputPath, surfaceTrimmerOutputPath);
}

@end
