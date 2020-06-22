//
//  MeshingTests.m
//  MeshingTests
//
//  Created by Aaron Thompson on 4/26/18.
//  Copyright © 2018 Standard Cyborg. All rights reserved.
//

#import <XCTest/XCTest.h>
#import <Meshing/Meshing.h>

@interface MeshingTests : XCTestCase

@end

@implementation MeshingTests

- (void)testMeshingOperation {
    // https://stackoverflow.com/questions/26811170/how-to-create-a-single-shared-framework-between-ios-and-os-x
    NSString *inputPath = [[NSBundle bundleForClass:[self class]] pathForResource:@"app_scan" ofType:@"ply"];
    NSString *outputPath = [NSTemporaryDirectory() stringByAppendingPathComponent:@"poisson_reconstructed_trimmed.ply"];
    
    NSLog(@"Running Poisson on file %@", inputPath);
    
    MeshingOperation *operation = [[MeshingOperation alloc] initWithInputFilePath:inputPath outputFilePath:outputPath];
    operation.resolution = 5;
    operation.smoothness = 1;
    operation.surfaceTrimmingAmount = 7;
    operation.closed = YES;
    [operation start];
    
    NSLog(@"Finished with output at %@", outputPath);
    XCTAssertTrue([[NSFileManager defaultManager] fileExistsAtPath:outputPath]);
}

@end
