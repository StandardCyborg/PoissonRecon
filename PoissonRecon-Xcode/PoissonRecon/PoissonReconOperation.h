//
//  PoissonReconOperation.h
//  PoissonRecon
//
//  Created by Aaron Thompson on 4/26/18.
//  Copyright © 2018 Standard Cyborg. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface PoissonReconOperation : NSOperation

- (instancetype)init NS_UNAVAILABLE;
- (instancetype)initWithInputFilePath:(NSString *)inputPath
                       outputFilePath:(NSString *)outputPath;

@end
