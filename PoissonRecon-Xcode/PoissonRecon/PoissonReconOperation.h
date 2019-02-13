//
//  PoissonReconOperation.h
//  PoissonRecon
//
//  Created by Aaron Thompson on 4/26/18.
//  Copyright © 2018 Standard Cyborg. All rights reserved.
//

#import <Foundation/Foundation.h>

/** Reconstructs a mesh from a point cloud and trims the edges on the result. */
@interface PoissonReconOperation : NSOperation

- (instancetype)init NS_UNAVAILABLE;

- (instancetype)initWithInputFilePath:(NSString *)inputPath
                       outputFilePath:(NSString *)outputPath;

/** The resolution of the reconstructed mesh vertices.
    Higher values will result in more vertices per meshes,
    and also take longer to reconstruct.
    Range is 1-10, default is 5.
 */
@property (nonatomic) int resolution;

/** The smoothness of the reconstructed mesh vertex positions.
    Range is 1-10, default is 2.
 */
@property (nonatomic) int smoothness;

@end
