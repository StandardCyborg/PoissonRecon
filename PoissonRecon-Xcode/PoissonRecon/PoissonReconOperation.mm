//
//  PoissonReconOperation.mm
//  PoissonRecon
//
//  Created by Aaron Thompson on 4/26/18.
//  Copyright © 2018 Standard Cyborg. All rights reserved.
//

#import "PoissonReconOperation.h"
#import "PoissonReconExecute.hpp"
#import "SurfaceTrimmerExecute.hpp"

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
    
    // [self _runPoissonReconInWorkingPath:workingPath
    //                           inputPath:inputPath
    //                          outputPath:poissonOutputPath];
    //
    // [self _runSurfaceTrimmerInWorkingPath:workingPath
    //                             inputPath:poissonOutputPath
    //                            outputPath:surfaceTrimmerOutputPath];
    // 
    // [[NSFileManager defaultManager] removeItemAtPath:poissonOutputPathString error:NULL];
    
    PoissonReconExecute(inputPath, poissonOutputPath);
    SurfaceTrimmerExecute(poissonOutputPath, surfaceTrimmerOutputPath);
}

/*
- (void)_runPoissonReconInWorkingPath:(const char *)workingPath
                            inputPath:(const char *)inputPath
                           outputPath:(const char *)outputPath
{
    // Example usage: Bin/Linux/PoissonRecon --in /tmp/PointCloud.ply --out /tmp/PoissonRecon.ply --colors --normals --density --depth 12
    const char *argv[] = {
        workingPath,
        "--in",
        inputPath,
        "--out",
        outputPath,
        "--colors",
        "--normals",
        "--density",
        "--ascii"
    };
    int argc = sizeof(argv) / sizeof(const char *);
    
    PoissonRecon_main(argc, (char **)argv);
}
*/

/*
- (void)_runSurfaceTrimmerInWorkingPath:(const char *)workingPath
                              inputPath:(const char *)inputPath
                             outputPath:(const char *)outputPath
{
    // Example usage: Bin/Linux/SurfaceTrimmer --in /tmp/PoissonRecon.ply --out /tmp/SurfaceTrimmer.ply --trim 5
    const char *argv[] = {
        workingPath,
        "--in",
        inputPath,
        "--out",
        outputPath,
        "--trim"
    };
    int argc = 8;
    int argc = sizeof(argv) / sizeof(const char *);
 
    SurfaceTrimmer_main(argc, argv);
}
*/

@end
