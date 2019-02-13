//
//  PoissonReconExecute.cpp
//  PoissonRecon
//
//  Created by Aaron Thompson on 2/12/19.
//  Copyright © 2019 Standard Cyborg. All rights reserved.
//

/**
 Example parameters for execution of
 PoissonRecon --in Scan.ply --out PoissonRecon.ply --colors --normals --density --depth 12
 {
     ASCII: true,
     BaseDepth: 0,
     BaseVCycles: 1,
     BoundaryNames: ["free", "Dirichlet", "Neumann"],
     BType: 3,
     CGSolverAccuracy: 0.00100000005,
     Colors: true,
     Confidence: 0,
     ConfidenceBias: 0,
     DataX: 32,
     Degree: 1,
     Density: true,
     Depth 12,
     ExactInterpolation: false,
     FEMTreeRealNames: ["float", "double"],
     FullDepth: 5,
     Grid: NULL,
     In: "Scan.ply",
     InCore: false,
     Iters: 8,
     KernelDepth: 0,
     LinearFit: false,
     MaxMemoryGB: 0,
     messageWriter: {
         outputFile: NULL,
         echoSTDOUT: false,
     },
     NoComments: false,
     NonManifold: false,
     Normals: true,
     Out: "PoissonRecon.ply",
     Performance: false,
     PointWeight: 2,
     PolygonMesh: false,
     PrimalGrid: false,
     SamplesPerNode: 1.5,
     Scale: 1.10000002,
     ShowGlobalResidualNames: ["show none", "show last", "show all"],
     ShowResidual: false,
     TempDir: NULL,
     Threads: 1,
     Transform: NULL,
     Tree: NULL,
     Verbose: false,
     Width: 0
 }
 */

#import "PoissonReconExecute.hpp"
#import "MyMiscellany.h"
#import "FEMTree.h"
#import "Ply.h"
#import "PointStreamData.h"

#define DATA_DEGREE 0                            // The order of the B-Spline used to splat in data for color interpolation
#define WEIGHT_DEGREE 2                            // The order of the B-Spline used to splat in the weights for density estimation
#define NORMAL_DEGREE 2                            // The order of the B-Spline used to splat in the normals for constructing the Laplacian constraints
#define DEFAULT_FEM_DEGREE 1                    // The default finite-element degree
#define DEFAULT_FEM_BOUNDARY BOUNDARY_NEUMANN    // The default finite-element boundary type

inline XForm<float, 3 + 1> GetBoundingBoxXForm(Point<float, 3> min, Point<float, 3> max, float scaleFactor)
{
    Point<float, 3> center = (max + min) / 2;
    float scale = max[0] - min[0];
    for (int d = 1; d < 3; d++) {
        scale = std::max<float>(scale, max[d] - min[d]);
    }
    
    scale *= scaleFactor;
    
    for (int i = 0; i < 3; i++) {
        center[i] -= scale / 2;
    }
    
    auto tXForm = XForm<float, 3 + 1>::Identity();
    auto sXForm = XForm<float, 3 + 1>::Identity();
    
    for (int i = 0; i < 3; i++) {
        sXForm(i, i) = 1.0f / scale;
        tXForm(3, i) = -center[i];
    }
    
    return sXForm * tXForm;
}

inline XForm<float, 4> GetPointXForm(InputPointStream<float, 3>& stream, float scaleFactor)
{
    Point<float, 3> min, max;
    stream.boundingBox(min, max);
    return GetBoundingBoxXForm(min, max, scaleFactor);
}

struct ConstraintDual
{
    float target, weight;
    ConstraintDual(float t, float w) : target(t), weight(w) { }
    CumulativeDerivativeValues<float, 3, 0> operator()(const Point<float, 3>& p) const {
        return CumulativeDerivativeValues<float, 3, 0>(target * weight);
    };
};

struct SystemDual {
    float weight;
    SystemDual(float w)
    : weight(w)
    { }
    
    CumulativeDerivativeValues<float,  3, 0> operator()(const Point<float, 3>& p, const CumulativeDerivativeValues<float,  3, 0>& dValues) const {
        return dValues * weight;
    };
    
    CumulativeDerivativeValues<double, 3, 0> operator()(const Point<float, 3>& p, const CumulativeDerivativeValues<double, 3, 0>& dValues) const {
        return dValues * weight;
    };
};

template<typename Vertex, unsigned int ... FEMSigs, typename ... SampleData>
inline void ExtractMesh(UIntPack<FEMSigs ...>,
                        std::tuple<SampleData ...>,
                        FEMTree<sizeof ...(FEMSigs), float>& tree,
                        const DenseNodeData<float, UIntPack<FEMSigs ...>>& solution,
                        float isoValue,
                        const std::vector<typename FEMTree<sizeof ...(FEMSigs), float>::PointSample>* samples,
                        std::vector<MultiPointStreamData<float, PointStreamNormal<float, 3>, MultiPointStreamData<float, SampleData ...>>>* sampleData,
                        const typename FEMTree<sizeof ...(FEMSigs), float>::template DensityEstimator<WEIGHT_DEGREE>* density,
                        std::function<void (Vertex&, Point<float, 3>, float, MultiPointStreamData<float, PointStreamNormal<float, 3>, MultiPointStreamData<float, SampleData ...>>)> SetVertex,
                        MessageWriter& messageWriter,
                        std::vector<std::string> &comments,
                        XForm<float, sizeof...(FEMSigs) + 1> iXForm,
                        const char *Out)
{
    const bool ASCII = true;
    const int DataX = 32;
    const bool LinearFit = false;
    const bool NoComments = false;
    const bool NonManifold = false;
    const bool PolygonMesh = false;
    
    typedef UIntPack<FEMSigs ...> Sigs;
    typedef PointStreamNormal<float, 3> NormalPointSampleData;
    typedef MultiPointStreamData<float, SampleData ...> AdditionalPointSampleData;
    typedef MultiPointStreamData<float, NormalPointSampleData, AdditionalPointSampleData> TotalPointSampleData;
    static const unsigned int DataSig = FEMDegreeAndBType<DATA_DEGREE, BOUNDARY_FREE>::Signature;
    typedef typename FEMTree<3, float>::template DensityEstimator<WEIGHT_DEGREE> DensityEstimator;
    
    char tempHeader[1024];
    {
        char tempPath[1024];
        tempPath[0] = 0;
        SetTempDirectory(tempPath, sizeof(tempPath));
        
        if (strlen(tempPath) == 0) sprintf(tempPath, ".%c", FileSeparator);
        if (tempPath[strlen(tempPath) - 1] == FileSeparator) sprintf(tempHeader, "%sPR_", tempPath);
        else                                                 sprintf(tempHeader, "%s%cPR_", tempPath, FileSeparator);
    }
    CoredMeshData<Vertex> *mesh = new CoredFileMeshData<Vertex>(tempHeader);
    
    typename IsoSurfaceExtractor<3, float, Vertex>::IsoStats isoStats;
    
    if (sampleData)
    {
        SparseNodeData<ProjectiveData<TotalPointSampleData, float>, IsotropicUIntPack<3, DataSig>> _sampleData = tree.template setDataField<DataSig, false>(*samples, *sampleData, (DensityEstimator*)NULL);
        for (const RegularTreeNode<3, FEMTreeNodeData>* n = tree.tree().nextNode(); n; n = tree.tree().nextNode(n))
        {
            ProjectiveData<TotalPointSampleData, float>* clr = _sampleData(n);
            if (clr) { (*clr) *= (float)pow(DataX, tree.depth(n)); }
        }
        
        isoStats = IsoSurfaceExtractor<3, float, Vertex>::template Extract<TotalPointSampleData>(Sigs(), UIntPack<WEIGHT_DEGREE>(), UIntPack<DataSig>(), tree, density, &_sampleData, solution, isoValue, *mesh, SetVertex, !LinearFit, !NonManifold, PolygonMesh, false);
    }
#if defined(__GNUC__) && __GNUC__ < 5
// #warning "you've got me gcc version < 5"
    else isoStats = IsoSurfaceExtractor<3, float, Vertex>::template Extract<TotalPointSampleData>(Sigs(), UIntPack<WEIGHT_DEGREE>(), UIntPack<DataSig>(), tree, density, (SparseNodeData<ProjectiveData<TotalPointSampleData, float>, IsotropicUIntPack<3, DataSig>> *)NULL, solution, isoValue, *mesh, SetVertex, !LinearFit, !NonManifold, PolygonMesh, false);
#else // !__GNUC__ || __GNUC__ >= 5
    else isoStats = IsoSurfaceExtractor<3, float, Vertex>::template Extract<TotalPointSampleData>(Sigs(), UIntPack<WEIGHT_DEGREE>(), UIntPack<DataSig>(), tree, density, NULL, solution, isoValue, *mesh, SetVertex, !LinearFit, !NonManifold, PolygonMesh, false);
#endif // __GNUC__ || __GNUC__ < 4
    messageWriter("Vertices / Polygons: %d / %d\n", mesh->outOfCorePointCount() + mesh->inCorePoints.size(), mesh->polygonCount());
    std::string isoStatsString = isoStats.toString() + std::string("\n");
    messageWriter(isoStatsString.c_str());
    
    std::vector<std::string> noComments;
    if (!PlyWritePolygons<Vertex, float, 3>(Out, mesh, ASCII ? PLY_ASCII : PLY_BINARY_NATIVE, NoComments ? noComments : comments, iXForm)) {
        ERROR_OUT("Could not write mesh to: %s", Out);
    }
    
    delete mesh;
}

// Called templated as Execute<float, PointStreamColor<float>>(argc, argv, IsotropicUIntPack<3, FEMDegreeAndBType<1, BOUNDARY_NEUMANN>::Signature>())
template<typename ... SampleData, unsigned int ... FEMSigs>
inline void Execute(const char *In, const char *Out, UIntPack<FEMSigs ...>)
{
    const int BaseDepth = 0;
    const int BaseVCycles = 1;
    const float CGSolverAccuracy = 0.001;
    const int Depth = 8;
    const int FullDepth = 5;
    const int Iters = 8;
    const int PointWeight = 2;
    const float SamplesPerNode = 1.5;
    const float Scale = 1.1;
    const char *ext = "ply";
    
    typedef UIntPack<FEMSigs ...> Sigs;
    typedef UIntPack<FEMSignature<FEMSigs>::Degree ...> Degrees;
    typedef UIntPack<FEMDegreeAndBType<NORMAL_DEGREE, DerivativeBoundary<FEMSignature<FEMSigs>::BType, 1>::BType>::Signature ...> NormalSigs;
    typedef typename FEMTree<3, float>::template DensityEstimator<WEIGHT_DEGREE> DensityEstimator;
    typedef typename FEMTree<3, float>::template InterpolationInfo<float, 0> InterpolationInfo;
    typedef PointStreamNormal<float, 3> NormalPointSampleData;
    typedef MultiPointStreamData<float, SampleData ...> AdditionalPointSampleData;
    typedef MultiPointStreamData<float, NormalPointSampleData, AdditionalPointSampleData> TotalPointSampleData;
    typedef InputPointStreamWithData<float, 3, TotalPointSampleData> InputPointStream;
    typedef TransformedInputPointStreamWithData<float, 3, TotalPointSampleData> XInputPointStream;
    MessageWriter messageWriter;
    std::vector<std::string> comments;
    
    messageWriter(comments, "*************************************************************\n");
    messageWriter(comments, "** Running Screened Poisson Reconstruction (Version %s) **\n", VERSION);
    messageWriter(comments, "*************************************************************\n");
    
    XForm<float, 3 + 1> xForm = XForm<float, 3 + 1>::Identity();
    
    float isoValue = 0;
    
    FEMTree<3, float> tree(MEMORY_ALLOCATOR_BLOCK_SIZE);
    
    int pointCount;
    
    float pointWeightSum;
    std::vector<typename FEMTree<3, float>::PointSample>* samples = new std::vector<typename FEMTree<3, float>::PointSample>();
    std::vector<TotalPointSampleData>* sampleData = NULL;
    DensityEstimator* density = NULL;
    SparseNodeData<Point<float, 3>, NormalSigs>* normalInfo = NULL;
    float targetValue = (float)0.5;
    
    // Read in the samples (and color data)
    {
        InputPointStream* pointStream;
        sampleData = new std::vector<TotalPointSampleData>();
        std::vector<std::pair<Point<float, 3>, TotalPointSampleData>> inCorePoints;
        if (!strcasecmp(ext, "ply")) {
            pointStream = new PLYInputPointStreamWithData  <float, 3, TotalPointSampleData>(
                In,
                TotalPointSampleData::PlyReadProperties(),
                TotalPointSampleData::PlyReadNum,
                TotalPointSampleData::ValidPlyReadProperties);
        } else {
            pointStream = new ASCIIInputPointStreamWithData<float, 3, TotalPointSampleData>(
                In,
                TotalPointSampleData::ReadASCII);
        }
        
        typename TotalPointSampleData::Transform _xForm(xForm);
        XInputPointStream _pointStream([&](Point<float, 3>& p, TotalPointSampleData& d){ p = xForm * p, d = _xForm(d); }, *pointStream);
        xForm = Scale > 0 ? GetPointXForm(_pointStream, (float)Scale) * xForm : xForm;
        {
            typename TotalPointSampleData::Transform _xForm(xForm);
            XInputPointStream _pointStream([&](Point<float, 3>& p, TotalPointSampleData& d){ p = xForm*p, d = _xForm(d); }, *pointStream);
            auto ProcessData = [](const Point<float, 3>& p, TotalPointSampleData& d)
            {
                float l = (float)Length(std::get<0>(d.data).data);
                
                if (!l || l != l) { return -1.0f; }
                
                std::get<0>(d.data).data /= l;
                
                return 1.0f;
            };
            
            pointCount = FEMTreeInitializer<3, float>::template Initialize<TotalPointSampleData>(tree.spaceRoot(), _pointStream, Depth, *samples, *sampleData, true, tree.nodeAllocator, tree.initializer(), ProcessData);
        }
        
        delete pointStream;
        
        messageWriter("Input Points / Samples: %d / %d\n", pointCount, samples->size());
    }
    
    int kernelDepth = Depth - 2;
    
    DenseNodeData<float, Sigs> solution;
    {
        DenseNodeData<float, Sigs> constraints;
        InterpolationInfo* iInfo = NULL;
        int solveDepth = Depth;
        
        tree.resetNodeIndices();
        
        // Get the kernel density estimator
        {
            density = tree.template setDensityEstimator<WEIGHT_DEGREE>(*samples, kernelDepth, SamplesPerNode, 1);
        }
        
        // Transform the Hermite samples into a vector field
        {
            normalInfo = new SparseNodeData<Point<float, 3>, NormalSigs>();
            *normalInfo = tree.setNormalField(NormalSigs(), *samples, *sampleData, density, pointWeightSum);
#pragma omp parallel for
            for (int i = 0; i < normalInfo->size(); i++) { (*normalInfo)[i] *= -1.0f; }
            messageWriter("Point weight / Estimated Area: %g / %g\n", pointWeightSum, pointCount * pointWeightSum);
        }
        
        // Trim the tree and prepare for multigrid
        {
            constexpr int MAX_DEGREE = NORMAL_DEGREE> Degrees::Max() ? NORMAL_DEGREE : Degrees::Max();
            tree.template finalizeForMultigrid<MAX_DEGREE>(FullDepth, typename FEMTree<3, float>::template HasNormalDataFunctor<NormalSigs>(*normalInfo), normalInfo, density);
        }
        
        // Add the FEM constraints
        {
            constraints = tree.initDenseNodeData(Sigs());
            typename FEMIntegrator::template Constraint<Sigs, IsotropicUIntPack<3, 1>, NormalSigs, IsotropicUIntPack<3, 0>, 3> F;
            unsigned int derivatives2[3];
            for (int d = 0; d < 3; d++) { derivatives2[d] = 0; }
            
            typedef IsotropicUIntPack<3, 1> Derivatives1;
            typedef IsotropicUIntPack<3, 0> Derivatives2;
            
            for (int d = 0; d < 3; d++)
            {
                unsigned int derivatives1[3];
                for (int dd = 0; dd < 3; dd++) { derivatives1[dd] = dd == d ? 1 : 0; }
                
                F.weights[d][TensorDerivatives<Derivatives1>::Index(derivatives1)][TensorDerivatives<Derivatives2>::Index(derivatives2)] = 1;
            }
            
            tree.addFEMConstraints(F, *normalInfo, constraints, solveDepth);
        }
        
        // Free up the normal info
        delete normalInfo, normalInfo = NULL;
        
        // Add the interpolation constraints
        if (PointWeight > 0)
        {
            iInfo = FEMTree<3, float>::template InitializeApproximatePointInterpolationInfo<float, 0>(
                         tree,
                         *samples,
                         ConstraintDual(targetValue, (float)PointWeight * pointWeightSum),
                         SystemDual((float)PointWeight * pointWeightSum),
                         true,
                         1);
            tree.addInterpolationConstraints(constraints, solveDepth, *iInfo);
        }
        
        messageWriter("Leaf Nodes / Active Nodes / Ghost Nodes: %d / %d / %d\n", (int)tree.leaves(), (int)tree.nodes(), (int)tree.ghostNodes());
        
        // Solve the linear system
        {
            typename FEMTree<3, float>::SolverInfo sInfo;
            sInfo.cgDepth = 0;
            sInfo.cascadic = true;
            sInfo.vCycles = 1;
            sInfo.iters = Iters;
            sInfo.cgAccuracy = CGSolverAccuracy;
            sInfo.verbose = false;
            sInfo.showResidual = false;
            sInfo.showGlobalResidual = SHOW_GLOBAL_RESIDUAL_NONE;
            sInfo.sliceBlockSize = 1;
            sInfo.baseDepth = BaseDepth;
            sInfo.baseVCycles = BaseVCycles;
            
            typename FEMIntegrator::template System<Sigs, IsotropicUIntPack<3, 1>> F({ 0.0, 1.0 });
            solution = tree.solveSystem(Sigs(), F, constraints, solveDepth, sInfo, iInfo);
            if (iInfo) delete iInfo, iInfo = NULL;
        }
    }
    
    {
        double valueSum = 0, weightSum = 0;
        typename FEMTree<3, float>::template MultiThreadedEvaluator<Sigs, 0> evaluator(&tree, solution);
#pragma omp parallel for reduction(+ : valueSum, weightSum)
        for (int j = 0; j < samples->size(); j++)
        {
            ProjectiveData<Point<float, 3>, float>& sample = (*samples)[j].sample;
            float w = sample.weight;
            if (w > 0) {
                weightSum += w;
                valueSum += evaluator.values(sample.data / sample.weight, omp_get_thread_num(), (*samples)[j].node)[0] * w;
            }
        }
        
        isoValue = (float)(valueSum / weightSum);
        messageWriter("Iso-Value: %e = %g / %g\n", isoValue, valueSum, weightSum);
    }
    
    typedef PlyVertexWithData<float, 3, MultiPointStreamData<float, PointStreamNormal<float, 3>, PointStreamValue<float>, AdditionalPointSampleData>> Vertex;
    std::function<void (Vertex&, Point<float, 3>, float, TotalPointSampleData)> SetVertex = [](Vertex& v, Point<float, 3> p, float w, TotalPointSampleData d) {
        v.point = p;
        std::get<0>(v.data.data) = std::get<0>(d.data);
        std::get<1>(v.data.data).data = w;
        std::get<2>(v.data.data) = std::get<1>(d.data);
    };
    
    ExtractMesh<Vertex>(UIntPack<FEMSigs ...>(),
                        std::tuple<SampleData ...>(),
                        tree,
                        solution,
                        isoValue,
                        samples,
                        sampleData,
                        density,
                        SetVertex,
                        messageWriter,
                        comments,
                        xForm.inverse(),
                        Out);
    
    if (sampleData) { delete sampleData; sampleData = NULL; }
    if (density) { delete density, density = NULL; }
}

void PoissonReconExecute(const char *inputFilePath, const char *outputFilePath)
{
    Execute<PointStreamColor<float>>(inputFilePath, outputFilePath, IsotropicUIntPack<3, FEMDegreeAndBType<1, BOUNDARY_NEUMANN>::Signature>());
}























/**
 Command-line parameters:
 {
 "In": "/path/to/input.ply",
 "Out": "/path/to/output.ply",
 "Smooth": 5,
 "Trim": 5,
 "IslandAreaRatio": 0.001,
 "PolygonMesh": false,
 "Verbose": false
 }
 */

static long long EdgeKey(int key1, int key2)
{
    if (key1 < key2) {
        return (((long long)key1) << 32) | ((long long)key2);
    } else {
        return (((long long)key2) << 32) | ((long long)key1);
    }
}

template <typename ... VertexData>
static PlyVertexWithData<float, 3, MultiPointStreamData<float, PointStreamValue<float>, VertexData...>>
InterpolateVertices(const PlyVertexWithData<float, 3, MultiPointStreamData<float, PointStreamValue<float>, VertexData...>>& v1,
                    const PlyVertexWithData<float, 3, MultiPointStreamData<float, PointStreamValue<float>, VertexData...>>& v2,
                    float value)
{
    if (std::get<0>(v1.data.data).data == std::get<0>(v2.data.data).data) {
        return (v1 + v2) / float(2.0);
    }
    
    float dx = (std::get<0>(v1.data.data).data - value) / (std::get<0>(v1.data.data).data - std::get<0>(v2.data.data).data);
    
    return v1 * (1.f - dx) + v2 * dx;
}

template <typename ... VertexData>
static void SmoothValues(std::vector<PlyVertexWithData<float, 3, MultiPointStreamData<float, PointStreamValue<float>, VertexData...>>>& vertices,
                         const std::vector<std::vector<int>>& polygons)
{
    std::vector<int> count(vertices.size());
    std::vector<float> sums(vertices.size(), 0);
    
    for (size_t i = 0; i < polygons.size(); i++) {
        int sz = int(polygons[i].size());
        
        for (int j = 0; j < sz; j++) {
            int j1 = j;
            int j2 = (j + 1) % sz;
            int v1 = polygons[i][j1];
            int v2 = polygons[i][j2];
            count[v1]++;
            count[v2]++;
            sums[v1] += std::get<0>(vertices[v2].data.data).data;
            sums[v2] += std::get<0>(vertices[v1].data.data).data;
        }
    }
    
    for (size_t i = 0; i < vertices.size(); i++) {
        std::get<0>(vertices[i].data.data).data = (sums[i] + std::get<0>(vertices[i].data.data).data) / (count[i] + 1);
    }
}

template <typename ... VertexData>
static void SplitPolygon(
                         const std::vector<int>& polygon,
                         std::vector<PlyVertexWithData<float, 3, MultiPointStreamData<float, PointStreamValue<float>, VertexData...>>>& vertices,
                         std::vector<std::vector<int>>* ltPolygons, std::vector<std::vector<int>>* gtPolygons,
                         std::vector<bool>* ltFlags, std::vector<bool>* gtFlags,
                         std::unordered_map<long long, int>& vertexTable,
                         float trimValue)
{
    int sz = int(polygon.size());
    std::vector<bool> gt(sz);
    int gtCount = 0;
    
    for (int j = 0; j < sz; j++) {
        gt[j] = (std::get<0>(vertices[polygon[j]].data.data).data > trimValue);
        
        if (gt[j]) {
            gtCount++;
        }
    }
    
    if (gtCount == sz) {
        if (gtPolygons) {
            gtPolygons->push_back(polygon);
        }
        if (gtFlags) {
            gtFlags->push_back(false);
        }
    }
    else if (gtCount == 0) {
        if (ltPolygons) {
            ltPolygons->push_back(polygon);
        }
        if (ltFlags) {
            ltFlags->push_back(false);
        }
    }
    else {
        int start;
        for (start = 0; start < sz; start++) {
            if (gt[start] && !gt[(start + sz - 1) % sz]) {
                break;
            }
        }
        
        bool gtFlag = true;
        std::vector<int> poly;
        
        // Add the initial vertex
        {
            int j1 = (start + int(sz) - 1) % sz, j2 = start;
            int v1 = polygon[j1], v2 = polygon[j2];
            int vIdx;
            std::unordered_map<long long, int>::iterator iter = vertexTable.find(EdgeKey(v1, v2));
            
            if (iter == vertexTable.end()) {
                vertexTable[EdgeKey(v1, v2)] = vIdx = int(vertices.size());
                vertices.push_back(InterpolateVertices(vertices[v1], vertices[v2], trimValue));
            } else {
                vIdx = iter->second;
            }
            
            poly.push_back(vIdx);
        }
        
        for (int _j = 0; _j <= sz; _j++) {
            int j1 = (_j + start + sz - 1) % sz, j2 = (_j + start) % sz;
            int v1 = polygon[j1], v2 = polygon[j2];
            
            if (gt[j2] == gtFlag) {
                poly.push_back(v2);
            } else {
                int vIdx;
                std::unordered_map<long long, int>::iterator iter = vertexTable.find(EdgeKey(v1, v2));
                
                if (iter == vertexTable.end()) {
                    vertexTable[EdgeKey(v1, v2)] = vIdx = int(vertices.size());
                    vertices.push_back(InterpolateVertices(vertices[v1], vertices[v2], trimValue));
                } else {
                    vIdx = iter->second;
                }
                
                poly.push_back(vIdx);
                
                if (gtFlag)
                {
                    if (gtPolygons) {
                        gtPolygons->push_back(poly);
                    }
                    if (ltFlags) {
                        ltFlags->push_back(true);
                    }
                }
                else
                {
                    if (ltPolygons) {
                        ltPolygons->push_back(poly);
                    }
                    if (gtFlags) {
                        gtFlags->push_back(true);
                    }
                }
                poly.clear(), poly.push_back(vIdx), poly.push_back(v2);
                gtFlag = !gtFlag;
            }
        }
    }
}

template <class Vertex>
static void Triangulate(const std::vector<Vertex>& vertices, const std::vector<std::vector<int>>& polygons, std::vector<std::vector<int>>& triangles)
{
    triangles.clear();
    
    for (size_t i = 0; i < polygons.size(); i++) {
        if (polygons.size() > 3) {
            std::vector<Point<float, 3>> _vertices(polygons[i].size());
            for (int j = 0; j < int(polygons[i].size()); j++) {
                _vertices[j] = vertices[polygons[i][j]].point;
            }
            
            std::vector<TriangleIndex> _triangles = MinimalAreaTriangulation<float, 3>((ConstPointer(Point<float, 3>))GetPointer(_vertices), _vertices.size());
            
            // Add the triangles to the mesh
            size_t idx = triangles.size();
            triangles.resize(idx + _triangles.size());
            
            for (int j = 0; j < int(_triangles.size()); j++) {
                triangles[idx + j].resize(3);
                for (int k = 0; k < 3; k++) {
                    triangles[idx + j][k] = polygons[i][_triangles[j].idx[k]];
                }
            }
        }
        else if (polygons[i].size() == 3) {
            triangles.push_back(polygons[i]);
        }
    }
}

template <class Vertex>
static double PolygonArea(const std::vector<Vertex>& vertices, const std::vector<int>& polygon)
{
    if (polygon.size() < 3) {
        return 0.0;
    }
    else if (polygon.size() == 3) {
        return Area(vertices[polygon[0]].point, vertices[polygon[1]].point, vertices[polygon[2]].point);
    }
    else {
        Point<float, 3> center;
        for (size_t i = 0; i < polygon.size(); i++) {
            center += vertices[polygon[i]].point;
        }
        
        center /= float(polygon.size());
        double area = 0;
        
        for (size_t i = 0; i < polygon.size(); i++) {
            area += Area(center, vertices[polygon[i]].point, vertices[polygon[(i + 1) % polygon.size()]].point);
        }
        
        return area;
    }
}

template <class Vertex>
static void RemoveHangingVertices(std::vector<Vertex>& vertices, std::vector<std::vector<int>>& polygons)
{
    std::unordered_map<int, int> vMap;
    std::vector<bool> vertexFlags(vertices.size(), false);
    
    for (size_t i = 0; i < polygons.size(); i++) {
        for (size_t j = 0; j < polygons[i].size(); j++) {
            vertexFlags[polygons[i][j]] = true;
        }
    }
    
    int vCount = 0;
    for (int i = 0; i < int(vertices.size()); i++) {
        if (vertexFlags[i]) {
            vMap[i] = vCount++;
        }
    }
    
    for (size_t i = 0; i < polygons.size(); i++) {
        for (size_t j = 0; j < polygons[i].size(); j++) {
            polygons[i][j] = vMap[polygons[i][j]];
        }
    }
    
    std::vector<Vertex> _vertices(vCount);
    for (int i = 0; i < int(vertices.size()); i++) {
        if (vertexFlags[i]) {
            _vertices[vMap[i]] = vertices[i];
        }
    }
    
    vertices = _vertices;
}

static void
SetConnectedComponents(const std::vector<std::vector<int>>& polygons,
                       std::vector<std::vector<int>>& components)
{
    std::vector<int> polygonRoots(polygons.size());
    for (size_t i = 0; i < polygons.size(); i++) {
        polygonRoots[i] = int(i);
    }
    
    std::unordered_map<long long, int> edgeTable;
    for (size_t i = 0; i < polygons.size(); i++) {
        int sz = int(polygons[i].size());
        
        for (int j = 0; j < sz; j++) {
            int j1 = j, j2 = (j + 1) % sz;
            int v1 = polygons[i][j1], v2 = polygons[i][j2];
            long long eKey = EdgeKey(v1, v2);
            std::unordered_map<long long, int>::iterator iter = edgeTable.find(eKey);
            
            if (iter == edgeTable.end()) {
                edgeTable[eKey] = int(i);
            } else {
                int p = iter->second;
                while (polygonRoots[p] != p) {
                    int temp = polygonRoots[p];
                    polygonRoots[p] = int(i);
                    p = temp;
                }
                polygonRoots[p] = int(i);
            }
        }
    }
    
    for (size_t i = 0; i < polygonRoots.size(); i++) {
        int p = int(i);
        while (polygonRoots[p] != p) {
            p = polygonRoots[p];
        }
        
        int root = p;
        p = int(i);
        
        while (polygonRoots[p] != p) {
            int temp = polygonRoots[p];
            polygonRoots[p] = root;
            p = temp;
        }
    }
    
    int cCount = 0;
    std::unordered_map<int, int> vMap;
    for (int i = 0; i < int(polygonRoots.size()); i++) {
        if (polygonRoots[i] == i) {
            vMap[i] = cCount++;
        }
    }
    
    components.resize(cCount);
    for (int i = 0; i < int(polygonRoots.size()); i++) {
        components[vMap[polygonRoots[i]]].push_back(i);
    }
}

template <typename ... VertexData>
static inline int _Execute(const char *In, const char *Out)
{
    const int Smooth = 5;
    const int Trim = 5;
    const float IslandAreaRatio = 0.001;
    const bool PolygonMesh = false;
    
    MessageWriter messageWriter;
    typedef PlyVertexWithData<float, 3, MultiPointStreamData<float, PointStreamValue<float>, VertexData ...>> Vertex;
    float min, max;
    std::vector<Vertex> vertices;
    std::vector<std::vector<int>> polygons;
    
    int ft;
    std::vector<std::string> comments;
    PlyReadPolygons<Vertex>(In, vertices, polygons, Vertex::PlyReadProperties(), Vertex::PlyReadNum, ft, comments);
    
    for (int i = 0; i < Smooth; i++) {
        SmoothValues(vertices, polygons);
    }
    
    min = max = std::get<0>(vertices[0].data.data).data;
    for (size_t i = 0; i < vertices.size(); i++) {
        min = std::min<float>(min, std::get<0>(vertices[i].data.data).data), max = std::max<float>(max, std::get<0>(vertices[i].data.data).data);
    }
    
    std::unordered_map<long long, int> vertexTable;
    std::vector<std::vector<int>> ltPolygons, gtPolygons;
    std::vector<bool> ltFlags, gtFlags;
    
    messageWriter(comments, "*********************************************\n");
    messageWriter(comments, "** Running Surface Trimmer **\n");
    messageWriter(comments, "*********************************************\n");
    
    for (size_t i = 0; i < polygons.size(); i++) {
        SplitPolygon(polygons[i], vertices, &ltPolygons, &gtPolygons, &ltFlags, &gtFlags, vertexTable, Trim);
    }
    
    if (IslandAreaRatio > 0) {
        std::vector<std::vector<int>> _ltPolygons, _gtPolygons;
        std::vector<std::vector<int>> ltComponents, gtComponents;
        SetConnectedComponents(ltPolygons, ltComponents);
        SetConnectedComponents(gtPolygons, gtComponents);
        std::vector<double> ltAreas(ltComponents.size(), 0.0), gtAreas(gtComponents.size(), 0.0);
        std::vector<bool> ltComponentFlags(ltComponents.size(), false), gtComponentFlags(gtComponents.size(), false);
        double area = 0.0;
        
        for (size_t i = 0; i < ltComponents.size(); i++) {
            for (size_t j = 0; j < ltComponents[i].size(); j++) {
                ltAreas[i] += PolygonArea<Vertex>(vertices, ltPolygons[ltComponents[i][j]]);
                ltComponentFlags[i] = (ltComponentFlags[i] || ltFlags[ltComponents[i][j]]);
            }
            area += ltAreas[i];
        }
        
        for (size_t i = 0; i < gtComponents.size(); i++) {
            for (size_t j = 0; j < gtComponents[i].size(); j++) {
                gtAreas[i] += PolygonArea<Vertex>(vertices, gtPolygons[gtComponents[i][j]]);
                gtComponentFlags[i] = (gtComponentFlags[i] || gtFlags[gtComponents[i][j]]);
            }
            area += gtAreas[i];
        }
        
        for (size_t i = 0; i < ltComponents.size(); i++) {
            if (ltAreas[i] < area * IslandAreaRatio && ltComponentFlags[i]) {
                for (size_t j = 0; j < ltComponents[i].size(); j++) {
                    _gtPolygons.push_back(ltPolygons[ltComponents[i][j]]);
                }
            }
            else {
                for (size_t j = 0; j < ltComponents[i].size(); j++) {
                    _ltPolygons.push_back(ltPolygons[ltComponents[i][j]]);
                }
            }
        }
        
        for (size_t i = 0; i < gtComponents.size(); i++) {
            if (gtAreas[i] < area * IslandAreaRatio && gtComponentFlags[i]) {
                for (size_t j = 0; j < gtComponents[i].size(); j++) {
                    _ltPolygons.push_back(gtPolygons[gtComponents[i][j]]);
                }
            }
            else {
                for (size_t j = 0; j < gtComponents[i].size(); j++) {
                    _gtPolygons.push_back(gtPolygons[gtComponents[i][j]]);
                }
            }
        }
        
        ltPolygons = _ltPolygons;
        gtPolygons = _gtPolygons;
    }
    
    if (!PolygonMesh) {
        {
            std::vector<std::vector<int>> polys = ltPolygons;
            Triangulate<Vertex>(vertices, ltPolygons, polys);
            ltPolygons = polys;
        }
        {
            std::vector<std::vector<int>> polys = gtPolygons;
            Triangulate<Vertex>(vertices, gtPolygons, polys);
            gtPolygons = polys;
        }
    }
    
    RemoveHangingVertices(vertices, gtPolygons);
    
    if (!PlyWritePolygons<Vertex>(Out, vertices, gtPolygons, Vertex::PlyWriteProperties(), Vertex::PlyWriteNum, ft, comments)) {
        ERROR_OUT("Could not write mesh to: %s", Out);
    }
    
    return EXIT_SUCCESS;
}

int SurfaceTrimmerExecute(const char* inputFilePath, const char* outputFilePath)
{
    typedef MultiPointStreamData<float, PointStreamValue<float>, PointStreamNormal<float, 3>, PointStreamColor<float>> VertexData;
    typedef PlyVertexWithData<float, 3, VertexData> Vertex;
    bool readFlags[Vertex::PlyReadNum];
    
    if (!PlyReadHeader((char *)inputFilePath, Vertex::PlyReadProperties(), Vertex::PlyReadNum, readFlags)) {
        ERROR_OUT("Failed to read ply header: %s", inputFilePath);
    }
    
    if (!VertexData::ValidPlyReadProperties<0>(readFlags + 3)) {
        ERROR_OUT("Ply file does not contain values");
    }
    
    return _Execute<PointStreamNormal<float, 3>, PointStreamColor<float>>(inputFilePath, outputFilePath);
}
