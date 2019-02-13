//
//  PoissonReconExecute.mm
//  PoissonRecon
//
//  Created by Aaron Thompson on 5/8/18.
//  Copyright © 2018 Standard Cyborg. All rights reserved.
//


#undef ARRAY_DEBUG                              // If enabled, array access is tested for validity
#define DATA_DEGREE 0                           // The order of the B-Spline used to splat in data for color interpolation
// This can be changed to zero if more interpolatory performance is desired.
#define WEIGHT_DEGREE 2                         // The order of the B-Spline used to splat in the weights for density estimation
#define NORMAL_DEGREE 2                         // The order of the B-Spline used to splat in the normals for constructing the Laplacian constraints
#define DEFAULT_FEM_DEGREE 1                    // The default finite-element degree
#define DEFAULT_FEM_BOUNDARY BOUNDARY_NEUMANN   // The default finite-element boundary type
#define DIMENSION 3                             // The dimension of the system

#import "PoissonReconExecute.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "MyMiscellany.h"
#include "CmdLineParser.h"
#include "PPolynomial.h"
#include "FEMTree.h"
#include "Ply.h"
#include "PointStreamData.h"
#include "Image.h"

MessageWriter messageWriter;

const float DefaultPointWeightMultiplier = 2.f;

cmdLineParameter<char*>
    InputPoints("in"),
    OutputMesh("out"),
    TempDir("tempDir"),
    Grid("grid"),
    Tree("tree"),
    Transform("xForm");

cmdLineReadable
    Performance("performance"),
    ShowResidual("showResidual"),
    NoComments("noComments"),
    PolygonMesh("polygonMesh"),
    NonManifold("nonManifold"),
    ASCII("ascii"),
    Density("density"),
    LinearFit("linearFit"),
    PrimalGrid("primalGrid"),
    ExactInterpolation("exact"),
    Normals("normals"),
    Colors("colors"),
    InCore("inCore"),
    Verbose("verbose");

cmdLineParameter<int>
    BSplineDegree("degree", DEFAULT_FEM_DEGREE),
    MaxDepth("depth", 8),
    KernelDepth("kernelDepth"),
    Iterations("iters", 8),
    FullDepth("fullDepth", 5),
    BaseMGSolverDepth("baseDepth", 0),
    BaseMGSolverVCycles("baseVCycles", 1),
    BoundaryType("bType", DEFAULT_FEM_BOUNDARY+1),
    MaxMemoryGB("maxMemory", 0),
    Threads("threads", omp_get_num_procs());

cmdLineParameter<float>
    DataX("data", 32.f),
    MinSamplesPerNode("samplesPerNode", 1.5f),
    Scale("scale", 1.1f),
    GridWidth("width", 0.f),
    Confidence("confidence", 0.f),
    ConfidenceBias("confidenceBias", 0.f),
    CGSolverAccuracy("cgAccuracy", 1e-3f),
    PointWeight("pointWeight");

cmdLineReadable* params[] =
{
    &BSplineDegree, &BoundaryType,
    &InputPoints, &MaxDepth, &OutputMesh, &Transform,
    &GridWidth,
    &Scale, &Verbose, &CGSolverAccuracy, &NoComments,
    &KernelDepth, &MinSamplesPerNode, &Confidence, &NonManifold, &PolygonMesh, &ASCII, &ShowResidual,
    &ConfidenceBias,
    &BaseMGSolverDepth, &BaseMGSolverVCycles,
    &PointWeight,
    &Grid, &Threads,
    &Tree,
    &Density,
    &FullDepth,
    &Iterations,
    &DataX,
    &Colors,
    &Normals,
    &LinearFit,
    &PrimalGrid,
    &TempDir,
    &ExactInterpolation,
    &Performance,
    &MaxMemoryGB,
    &InCore,
    NULL
};

void ShowUsage(char* ex)
{
    printf("Usage: %s\n", ex);
    printf("\t --%s <input points>\n", InputPoints.name);
    printf("\t[--%s <ouput triangle mesh>]\n", OutputMesh.name);
    printf("\t[--%s <ouput grid>]\n", Grid.name);
    printf("\t[--%s <ouput fem tree>]\n", Tree.name);
    printf("\t[--%s <b-spline degree> = %d]\n", BSplineDegree.name, BSplineDegree.value);
    printf("\t[--%s <boundary type> = %d]\n", BoundaryType.name, BoundaryType.value);
    for (int i = 0; i < BOUNDARY_COUNT; i++) printf("\t\t%d] %s\n", i+1, BoundaryNames[i]);
    printf("\t[--%s <maximum reconstruction depth> = %d]\n", MaxDepth.name, MaxDepth.value);
    printf("\t[--%s <grid width>]\n", GridWidth.name);
    printf("\t[--%s <full depth> = %d]\n", FullDepth.name, FullDepth.value);
    printf("\t[--%s <coarse MG solver depth> = %d]\n", BaseMGSolverDepth.name, BaseMGSolverDepth.value);
    printf("\t[--%s <coarse MG solver v-cycles> = %d]\n", BaseMGSolverVCycles.name, BaseMGSolverVCycles.value);
    printf("\t[--%s <scale factor> = %f]\n", Scale.name, Scale.value);
    printf("\t[--%s <minimum number of samples per node> = %f]\n", MinSamplesPerNode.name, MinSamplesPerNode.value);
    printf("\t[--%s <interpolation weight> = %.3e * <b-spline degree>]\n", PointWeight.name, DefaultPointWeightMultiplier);
    printf("\t[--%s <iterations> = %d]\n", Iterations.name, Iterations.value);
    printf("\t[--%s]\n", ExactInterpolation.name);
    printf("\t[--%s <pull factor> = %f]\n", DataX.name, DataX.value);
    printf("\t[--%s]\n", Colors.name);
    printf("\t[--%s]\n", Normals.name);
#ifdef _OPENMP
    printf("\t[--%s <num threads> = %d]\n", Threads.name, Threads.value);
#endif // _OPENMP
    printf("\t[--%s <normal confidence exponent> = %f]\n", Confidence.name, Confidence.value);
    printf("\t[--%s <normal confidence bias exponent> = %f]\n", ConfidenceBias.name, ConfidenceBias.value);
    printf("\t[--%s]\n", NonManifold.name);
    printf("\t[--%s]\n", PolygonMesh.name);
    printf("\t[--%s <cg solver accuracy> = %g]\n", CGSolverAccuracy.name, CGSolverAccuracy.value);
    printf("\t[--%s <maximum memory (in GB)> = %d]\n", MaxMemoryGB.name, MaxMemoryGB.value);
    printf("\t[--%s]\n", Performance.name);
    printf("\t[--%s]\n", Density.name);
    printf("\t[--%s]\n", LinearFit.name);
    printf("\t[--%s]\n", PrimalGrid.name);
    printf("\t[--%s]\n", ASCII.name);
    printf("\t[--%s]\n", NoComments.name);
    printf("\t[--%s]\n", TempDir.name);
    printf("\t[--%s]\n", InCore.name);
    printf("\t[--%s]\n", Verbose.name);
}

double Weight(double v, double start, double end)
{
    v = (v - start) / (end - start);
    if (v < 0) return 1.;
    else if (v > 1) return 0.;
    else
    {
        // P(x) = a x^3 + b x^2 + c x + d
        //        P (0) = 1, P (1) = 0, P'(0) = 0, P'(1) = 0
        // =>     d = 1, a + b + c + d = 0, c = 0, 3a + 2b + c = 0
        // =>     c = 0, d = 1, a + b = -1, 3a + 2b = 0
        // =>     a = 2, b = -3, c = 0, d = 1
        // =>     P(x) = 2 x^3 - 3 x^2 + 1
        return 2. * v * v * v - 3. * v * v + 1.;
    }
}

template<unsigned int Dim>
struct FEMTreeProfiler
{
    FEMTree<Dim, float>& tree;
    double t;
    
    FEMTreeProfiler(FEMTree<Dim, float>& t) : tree(t) {; }
    void start(void){ t = Time(), FEMTree<Dim, float>::ResetLocalMemoryUsage(); }
    void print(const char* header) const
    {
        FEMTree<Dim, float>::MemoryUsage();
        if (header) printf("%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n", header, Time()-t, FEMTree<Dim, float>::LocalMemoryUsage(), FEMTree<Dim, float>::MaxMemoryUsage(), MemoryInfo::PeakMemoryUsageMB());
        else         printf("%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n",          Time()-t, FEMTree<Dim, float>::LocalMemoryUsage(), FEMTree<Dim, float>::MaxMemoryUsage(), MemoryInfo::PeakMemoryUsageMB());
    }
    void dumpOutput(const char* header) const
    {
        FEMTree<Dim, float>::MemoryUsage();
        if (header) messageWriter("%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n", header, Time()-t, FEMTree<Dim, float>::LocalMemoryUsage(), FEMTree<Dim, float>::MaxMemoryUsage(), MemoryInfo::PeakMemoryUsageMB());
        else         messageWriter("%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n",          Time()-t, FEMTree<Dim, float>::LocalMemoryUsage(), FEMTree<Dim, float>::MaxMemoryUsage(), MemoryInfo::PeakMemoryUsageMB());
    }
    void dumpOutput2(std::vector<std::string>& comments, const char* header) const
    {
        FEMTree<Dim, float>::MemoryUsage();
        if (header) messageWriter(comments, "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n", header, Time()-t, FEMTree<Dim, float>::LocalMemoryUsage(), FEMTree<Dim, float>::MaxMemoryUsage(), MemoryInfo::PeakMemoryUsageMB());
        else         messageWriter(comments,    "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n",          Time()-t, FEMTree<Dim, float>::LocalMemoryUsage(), FEMTree<Dim, float>::MaxMemoryUsage(), MemoryInfo::PeakMemoryUsageMB());
    }
};

template<unsigned int Dim>
XForm<float, Dim + 1> GetBoundingBoxXForm(Point<float, Dim> min, Point<float, Dim> max, float scaleFactor)
{
    Point<float, Dim> center = (max + min) / 2;
    float scale = max[0] - min[0];
    for (int d = 1; d < Dim; d++) { scale = std::max<float>(scale, max[d] - min[d]); }
    scale *= scaleFactor;
    for (int i = 0; i < Dim; i++) { center[i] -= scale / 2; }
    XForm<float, Dim + 1> tXForm = XForm<float, Dim + 1>::Identity(), sXForm = XForm<float, Dim + 1>::Identity();
    for (int i = 0; i < Dim; i++) { sXForm(i,i) = (float)(1./scale), tXForm(Dim, i) = -center[i]; }
    return sXForm * tXForm;
}

template<unsigned int Dim>
XForm<float, Dim + 1> GetBoundingBoxXForm(Point<float, Dim> min, Point<float, Dim> max, float width, float scaleFactor, int& depth)
{
    // Get the target resolution (along the largest dimension)
    float resolution = (max[0]-min[0]) / width;
    for (int d = 1; d<Dim; d++) resolution = std::max<float>(resolution, (max[d]-min[d]) / width);
    resolution *= scaleFactor;
    depth = 0;
    while((1<<depth)<resolution) depth++;
    
    Point<float, Dim> center = (max + min) / 2;
    float scale = (1<<depth) * width;
    
    for (int i = 0; i < Dim; i++) center[i] -= scale/2;
    XForm<float, Dim + 1> tXForm = XForm<float, Dim + 1>::Identity(), sXForm = XForm<float, Dim + 1>::Identity();
    for (int i = 0; i < Dim; i++) sXForm(i,i) = (float)(1./scale), tXForm(Dim,i) = -center[i];
    return sXForm * tXForm;
}

template<unsigned int Dim>
XForm<float, Dim + 1> GetPointXForm(InputPointStream<float, Dim>& stream, float width, float scaleFactor, int& depth)
{
    Point<float, Dim> min, max;
    stream.boundingBox(min, max);
    return GetBoundingBoxXForm(min, max, width, scaleFactor, depth);
}
template<unsigned int Dim>
XForm<float, Dim + 1> GetPointXForm(InputPointStream<float, Dim>& stream, float scaleFactor)
{
    Point<float, Dim> min, max;
    stream.boundingBox(min, max);
    return GetBoundingBoxXForm(min, max, scaleFactor);
}

template<unsigned int Dim>
struct ConstraintDual
{
    float target, weight;
    ConstraintDual(float t, float w) : target(t), weight(w){ }
    CumulativeDerivativeValues<float, Dim, 0> operator()(const Point<float, Dim>& p) const { return CumulativeDerivativeValues<float, Dim, 0>(target*weight); };
};

template<unsigned int Dim>
struct SystemDual
{
    float weight;
    SystemDual(float w) : weight(w){ }
    CumulativeDerivativeValues<float, Dim, 0> operator()(const Point<float, Dim>& p, const CumulativeDerivativeValues<float, Dim, 0>& dValues) const { return dValues * weight; };
};

template<typename Vertex, unsigned int ... FEMSigs, typename ... SampleData>
void ExtractMesh(UIntPack<FEMSigs ...>, std::tuple<SampleData ...>, FEMTree<sizeof ... (FEMSigs), float>& tree, const DenseNodeData<float, UIntPack<FEMSigs ...>>& solution, float isoValue, const std::vector<typename FEMTree<sizeof ... (FEMSigs), float>::PointSample>* samples, std::vector<MultiPointStreamData<float, PointStreamNormal<float, DIMENSION>, MultiPointStreamData<float, SampleData ...>>>* sampleData, const typename FEMTree<sizeof ... (FEMSigs), float>::template DensityEstimator<WEIGHT_DEGREE>* density, std::function<void (Vertex&, Point<float, DIMENSION>, float, MultiPointStreamData<float, PointStreamNormal<float, DIMENSION>, MultiPointStreamData<float, SampleData ...>>)> SetVertex, std::vector<std::string> &comments, XForm<float, sizeof...(FEMSigs)+1> iXForm)
{
    static const int Dim = sizeof ... (FEMSigs);
    typedef UIntPack<FEMSigs ...> Sigs;
    typedef PointStreamNormal<float, Dim> NormalPointSampleData;
    typedef MultiPointStreamData<float, SampleData ...> AdditionalPointSampleData;
    typedef MultiPointStreamData<float, NormalPointSampleData, AdditionalPointSampleData> TotalPointSampleData;
    static const unsigned int DataSig = FEMDegreeAndBType<DATA_DEGREE, BOUNDARY_FREE>::Signature;
    typedef typename FEMTree<Dim, float>::template DensityEstimator<WEIGHT_DEGREE> DensityEstimator;
    
    FEMTreeProfiler<Dim> profiler(tree);
    
    char tempHeader[1024];
    {
        char tempPath[1024];
        tempPath[0] = 0;
        if (TempDir.set) strcpy(tempPath, TempDir.value);
        else SetTempDirectory(tempPath, sizeof(tempPath));
        if (strlen(tempPath) == 0) sprintf(tempPath, ".%c", FileSeparator);
        if (tempPath[ strlen(tempPath)-1 ] == FileSeparator) sprintf(tempHeader, "%sPR_", tempPath);
        else                                                  sprintf(tempHeader, "%s%cPR_", tempPath, FileSeparator);
    }
    CoredMeshData<Vertex> *mesh;
    if (InCore.set) mesh = new CoredVectorMeshData<Vertex>();
    else            mesh = new CoredFileMeshData<Vertex>(tempHeader);
    
    profiler.start();
    typename IsoSurfaceExtractor<Dim, float, Vertex>::IsoStats isoStats;
    if (sampleData)
    {
        SparseNodeData<ProjectiveData<TotalPointSampleData, float>, IsotropicUIntPack<Dim, DataSig>> _sampleData = tree.template setDataField<DataSig, false>(*samples, *sampleData, (DensityEstimator*)NULL);
        for (const RegularTreeNode<Dim, FEMTreeNodeData>* n = tree.tree().nextNode(); n; n = tree.tree().nextNode(n))
        {
            ProjectiveData<TotalPointSampleData, float>* clr = _sampleData(n);
            if (clr) (*clr) *= (float)pow(DataX.value, tree.depth(n));
        }
        isoStats = IsoSurfaceExtractor<Dim, float, Vertex>::template Extract<TotalPointSampleData>(Sigs(), UIntPack<WEIGHT_DEGREE>(), UIntPack<DataSig>(), tree, density, &_sampleData, solution, isoValue, *mesh, SetVertex, !LinearFit.set, !NonManifold.set, PolygonMesh.set, false);
    }
#if defined(__GNUC__) && __GNUC__ <5
#warning "you've got me gcc version<5"
    else isoStats = IsoSurfaceExtractor<Dim, float, Vertex>::template Extract<TotalPointSampleData>(Sigs(), UIntPack<WEIGHT_DEGREE>(), UIntPack<DataSig>(), tree, density, (SparseNodeData<ProjectiveData<TotalPointSampleData, float>, IsotropicUIntPack<Dim, DataSig>> *)NULL, solution, isoValue, *mesh, SetVertex, !LinearFit.set, !NonManifold.set, PolygonMesh.set, false);
#else // !__GNUC__ || __GNUC__ >= 5
    else isoStats = IsoSurfaceExtractor<Dim, float, Vertex>::template Extract<TotalPointSampleData>(Sigs(), UIntPack<WEIGHT_DEGREE>(), UIntPack<DataSig>(), tree, density, NULL, solution, isoValue, *mesh, SetVertex, !LinearFit.set, !NonManifold.set, PolygonMesh.set, false);
#endif // __GNUC__ || __GNUC__ < 4
    messageWriter("Vertices / Polygons: %d / %d\n", mesh->outOfCorePointCount()+mesh->inCorePoints.size(), mesh->polygonCount());
    std::string isoStatsString = isoStats.toString() + std::string("\n");
    messageWriter(isoStatsString.c_str());
    if (PolygonMesh.set) profiler.dumpOutput2(comments, "#         Got polygons:");
    else                  profiler.dumpOutput2(comments, "#        Got triangles:");
    
    std::vector<std::string> noComments;
    if (!PlyWritePolygons<Vertex, float, Dim>(OutputMesh.value, mesh, ASCII.set ? PLY_ASCII : PLY_BINARY_NATIVE, NoComments.set ? noComments : comments, iXForm))
        ERROR_OUT("Could not write mesh to: %s", OutputMesh.value);
    
    delete mesh;
}

template<unsigned int Dim>
void WriteGrid(ConstPointer(float) values, int res, const char *fileName)
{
    int resolution = 1;
    for (int d = 0; d<Dim; d++) resolution *= res;
    
    char *ext = GetFileExtension(fileName);
    
    if (Dim == 2 && ImageWriter::ValidExtension(ext))
    {
        float avg = 0;
#pragma omp parallel for reduction(+ : avg)
        for (int i = 0; i < resolution; i++) { avg += values[i]; }
        avg /= (float)resolution;
        
        float std = 0;
#pragma omp parallel for reduction(+ : std)
        for (int i = 0; i < resolution; i++) std += (values[i] - avg) * (values[i] - avg);
        std = (float)sqrt(std / resolution);
        
        if (Verbose.set) printf("Grid to image: [%.2f,%.2f] -> [0,255]\n", avg - 2*std, avg + 2*std);
        
        unsigned char *pixels = new unsigned char[resolution * 3];
#pragma omp parallel for
        for (int i = 0; i < resolution; i++)
        {
            float v = (float)std::min<float>((float)1., std::max<float>((float)-1., (values[i] - avg) / (2 * std)));
            v = (float)((v + 1.) / 2. * 256.);
            unsigned char color = (unsigned char)std::min<float>((float)255., std::max<float>((float)0., v));
            for (int c = 0; c < 3; c++) { pixels[i * 3 + c] = color; }
        }
        ImageWriter::Write(fileName, pixels, res, res, 3);
        delete[] pixels;
    }
    else
    {
        
        FILE *fp = fopen(fileName, "wb");
        if (!fp) ERROR_OUT("Failed to open grid file for writing: %s", fileName);
        else
        {
            fwrite(&res, sizeof(int), 1, fp);
            if (typeid(float) == typeid(float)) fwrite(values, sizeof(float), resolution, fp);
            else
            {
                float *fValues = new float[resolution];
                for (int i = 0; i < resolution; i++) fValues[i] = float(values[i]);
                fwrite(fValues, sizeof(float), resolution, fp);
                delete[] fValues;
            }
            fclose(fp);
        }
    }
    delete[] ext;
}


template<typename ... SampleData, unsigned int ... FEMSigs>
void Execute(int argc, char* argv[], UIntPack<FEMSigs ...>)
{
    static const int Dim = sizeof ... (FEMSigs);
    typedef UIntPack<FEMSigs ...> Sigs;
    typedef UIntPack<FEMSignature<FEMSigs>::Degree ...> Degrees;
    typedef UIntPack<FEMDegreeAndBType<NORMAL_DEGREE, DerivativeBoundary<FEMSignature<FEMSigs>::BType, 1>::BType>::Signature ...> NormalSigs;
    static const unsigned int DataSig = FEMDegreeAndBType<DATA_DEGREE, BOUNDARY_FREE>::Signature;
    typedef typename FEMTree<Dim, float>::template DensityEstimator<WEIGHT_DEGREE> DensityEstimator;
    typedef typename FEMTree<Dim, float>::template InterpolationInfo<float, 0> InterpolationInfo;
    typedef PointStreamNormal<float, Dim> NormalPointSampleData;
    typedef MultiPointStreamData<float, SampleData ...> AdditionalPointSampleData;
    typedef MultiPointStreamData<float, NormalPointSampleData, AdditionalPointSampleData> TotalPointSampleData;
    typedef InputPointStreamWithData<float, Dim, TotalPointSampleData> InputPointStream;
    typedef TransformedInputPointStreamWithData<float, Dim, TotalPointSampleData> XInputPointStream;
    std::vector<std::string> comments;
    messageWriter(comments, "*************************************************************\n");
    messageWriter(comments, "*************************************************************\n");
    messageWriter(comments, "** Running Screened Poisson Reconstruction (Version %s) **\n", VERSION);
    messageWriter(comments, "*************************************************************\n");
    messageWriter(comments, "*************************************************************\n");
    
    XForm<float, Dim + 1> xForm, iXForm;
    if (Transform.set)
    {
        FILE* fp = fopen(Transform.value, "r");
        if (!fp)
        {
            WARN("Could not read x-form from: %s", Transform.value);
            xForm = XForm<float, Dim + 1>::Identity();
        }
        else
        {
            for (int i = 0; i < Dim + 1; i++) for (int j = 0; j < Dim + 1; j++)
            {
                float f;
                if (fscanf(fp, " %f ", &f) != 1) ERROR_OUT("Failed to read xform");
                xForm(i,j) = (float)f;
            }
            fclose(fp);
        }
    }
    else xForm = XForm<float, Dim + 1>::Identity();
    
    char str[1024];
    for (int i = 0; params[i]; i++)
        if (params[i]->set)
        {
            params[i]->writeValue(str);
            if (strlen(str)) messageWriter(comments, "\t--%s %s\n", params[i]->name, str);
            else                messageWriter(comments, "\t--%s\n", params[i]->name);
        }
    
    double startTime = Time();
    float isoValue = 0;
    
    FEMTree<Dim, float> tree(MEMORY_ALLOCATOR_BLOCK_SIZE);
    FEMTreeProfiler<Dim> profiler(tree);
    
    if (MaxDepth.set && GridWidth.value > 0)
    {
        WARN("Both --%s and --%s set, ignoring --%s", MaxDepth.name, GridWidth.name, GridWidth.name);
        GridWidth.value = 0;
    }
    
    int pointCount;
    
    float pointWeightSum;
    std::vector<typename FEMTree<Dim, float>::PointSample>* samples = new std::vector<typename FEMTree<Dim, float>::PointSample>();
    std::vector<TotalPointSampleData>* sampleData = NULL;
    DensityEstimator* density = NULL;
    SparseNodeData<Point<float, Dim>, NormalSigs>* normalInfo = NULL;
    float targetValue = (float)0.5;
    
    // Read in the samples (and color data)
    {
        profiler.start();
        InputPointStream* pointStream;
        char* ext = GetFileExtension(InputPoints.value);
        sampleData = new std::vector<TotalPointSampleData>();
        std::vector<std::pair<Point<float, Dim>, TotalPointSampleData>> inCorePoints;
        if (InCore.set)
        {
            InputPointStream *_pointStream;
            if     (!strcasecmp(ext, "bnpts")) _pointStream = new BinaryInputPointStreamWithData<float, Dim, TotalPointSampleData>(InputPoints.value, TotalPointSampleData::ReadBinary);
            else if (!strcasecmp(ext, "ply"))  _pointStream = new PLYInputPointStreamWithData<float, Dim, TotalPointSampleData>(InputPoints.value, TotalPointSampleData::PlyReadProperties(), TotalPointSampleData::PlyReadNum, TotalPointSampleData::ValidPlyReadProperties);
            else                               _pointStream = new ASCIIInputPointStreamWithData<float, Dim, TotalPointSampleData>(InputPoints.value, TotalPointSampleData::ReadASCII);
            Point<float, Dim> p;
            TotalPointSampleData d;
            while(_pointStream->nextPoint(p, d)) inCorePoints.push_back(std::pair<Point<float, Dim>, TotalPointSampleData>(p, d));
            delete _pointStream;
            
            pointStream = new MemoryInputPointStreamWithData<float, Dim, TotalPointSampleData>(inCorePoints.size(), &inCorePoints[0]);
        }
        else
        {
            if     (!strcasecmp(ext, "bnpts")) pointStream = new BinaryInputPointStreamWithData<float, Dim, TotalPointSampleData>(InputPoints.value, TotalPointSampleData::ReadBinary);
            else if (!strcasecmp(ext, "ply"))  pointStream = new PLYInputPointStreamWithData<float, Dim, TotalPointSampleData>(InputPoints.value, TotalPointSampleData::PlyReadProperties(), TotalPointSampleData::PlyReadNum, TotalPointSampleData::ValidPlyReadProperties);
            else                               pointStream = new ASCIIInputPointStreamWithData<float, Dim, TotalPointSampleData>(InputPoints.value, TotalPointSampleData::ReadASCII);
        }
        delete[] ext;
        typename TotalPointSampleData::Transform _xForm(xForm);
        XInputPointStream _pointStream([&](Point<float, Dim>& p, TotalPointSampleData& d){ p = xForm*p, d = _xForm(d); }, *pointStream);
        if (GridWidth.value > 0) xForm = GetPointXForm<Dim>(_pointStream, GridWidth.value, (float)(Scale.value > 0 ? Scale.value : 1.), MaxDepth.value) * xForm;
        else                     xForm = Scale.value > 0 ? GetPointXForm<Dim>(_pointStream, (float)Scale.value) * xForm : xForm;
        {
            typename TotalPointSampleData::Transform _xForm(xForm);
            XInputPointStream _pointStream([&](Point<float, Dim>& p, TotalPointSampleData& d){ p = xForm*p, d = _xForm(d); }, *pointStream);
            auto ProcessDataWithConfidence = [&](const Point<float, Dim>& p, TotalPointSampleData& d)
            {
                float l = (float)Length(std::get<0>(d.data).data);
                if (!l || l != l) return (float)-1.;
                return (float)pow(l, Confidence.value);
            };
            auto ProcessData = [](const Point<float, Dim>& p, TotalPointSampleData& d)
            {
                float l = (float)Length(std::get<0>(d.data).data);
                if (!l || l != l) return (float)-1.;
                std::get<0>(d.data).data /= l;
                return (float)1.;
            };
            if (Confidence.value > 0) pointCount = FEMTreeInitializer<Dim, float>::template Initialize<TotalPointSampleData>(tree.spaceRoot(), _pointStream, MaxDepth.value, *samples, *sampleData, true, tree.nodeAllocator, tree.initializer(), ProcessDataWithConfidence);
            else                      pointCount = FEMTreeInitializer<Dim, float>::template Initialize<TotalPointSampleData>(tree.spaceRoot(), _pointStream, MaxDepth.value, *samples, *sampleData, true, tree.nodeAllocator, tree.initializer(), ProcessData);
        }
        iXForm = xForm.inverse();
        delete pointStream;
        
        messageWriter("Input Points / Samples: %d / %d\n", pointCount, samples->size());
        profiler.dumpOutput2(comments, "# Read input into tree:");
    }
    int kernelDepth = KernelDepth.set ? KernelDepth.value : MaxDepth.value-2;
    if (kernelDepth > MaxDepth.value)
    {
        WARN("%s can't be greater than %s: %d <= %d", KernelDepth.name, MaxDepth.name, KernelDepth.value, MaxDepth.value);
        kernelDepth = MaxDepth.value;
    }
    
    DenseNodeData<float, Sigs> solution;
    {
        DenseNodeData<float, Sigs> constraints;
        InterpolationInfo* iInfo = NULL;
        int solveDepth = MaxDepth.value;
        
        tree.resetNodeIndices();
        
        // Get the kernel density estimator
        {
            profiler.start();
            density = tree.template setDensityEstimator<WEIGHT_DEGREE>(*samples, kernelDepth, MinSamplesPerNode.value, 1);
            profiler.dumpOutput2(comments, "#   Got kernel density:");
        }
        
        // Transform the Hermite samples into a vector field
        {
            profiler.start();
            normalInfo = new SparseNodeData<Point<float, Dim>, NormalSigs>();
            if (ConfidenceBias.value > 0) *normalInfo = tree.setNormalField(NormalSigs(), *samples, *sampleData, density, pointWeightSum, [&](float conf){ return (float)(log(conf) * ConfidenceBias.value / log(1<<(Dim-1))); });
            else                          *normalInfo = tree.setNormalField(NormalSigs(), *samples, *sampleData, density, pointWeightSum);
#pragma omp parallel for
            for (int i = 0; i < normalInfo->size(); i++) {(*normalInfo)[i] *= (float)-1.; }
            profiler.dumpOutput2(comments, "#     Got normal field:");
            messageWriter("Point weight / Estimated Area: %g / %g\n", pointWeightSum, pointCount*pointWeightSum);
        }
        
        if (!Density.set) delete density, density = NULL;
        if (DataX.value <= 0 || (!Colors.set && !Normals.set)) delete sampleData, sampleData = NULL;
        
        // Trim the tree and prepare for multigrid
        {
            profiler.start();
            constexpr int MAX_DEGREE = NORMAL_DEGREE> Degrees::Max() ? NORMAL_DEGREE : Degrees::Max();
            tree.template finalizeForMultigrid<MAX_DEGREE>(FullDepth.value, typename FEMTree<Dim, float>::template HasNormalDataFunctor<NormalSigs>(*normalInfo), normalInfo, density);
            profiler.dumpOutput2(comments, "#       Finalized tree:");
        }
        // Add the FEM constraints
        {
            profiler.start();
            constraints = tree.initDenseNodeData(Sigs());
            typename FEMIntegrator::template Constraint<Sigs, IsotropicUIntPack<Dim, 1>, NormalSigs, IsotropicUIntPack<Dim, 0>, Dim> F;
            unsigned int derivatives2[Dim];
            for (int d = 0; d<Dim; d++) derivatives2[d] = 0;
            typedef IsotropicUIntPack<Dim, 1> Derivatives1;
            typedef IsotropicUIntPack<Dim, 0> Derivatives2;
            for (int d = 0; d<Dim; d++)
            {
                unsigned int derivatives1[Dim];
                for (int dd = 0; dd<Dim; dd++) derivatives1[dd] = dd == d ?  1 : 0;
                F.weights[d][ TensorDerivatives<Derivatives1>::Index(derivatives1) ][ TensorDerivatives<Derivatives2>::Index(derivatives2) ] = 1;
            }
            tree.addFEMConstraints(F, *normalInfo, constraints, solveDepth);
            profiler.dumpOutput2(comments, "#  Set FEM constraints:");
        }
        
        // Free up the normal info
        delete normalInfo, normalInfo = NULL;
        
        // Add the interpolation constraints
        if (PointWeight.value > 0)
        {
            profiler.start();
            if (ExactInterpolation.set) iInfo = FEMTree<Dim, float>::template       InitializeExactPointInterpolationInfo<float, 0> (tree, *samples, ConstraintDual<Dim>(targetValue, (float)PointWeight.value * pointWeightSum), SystemDual<Dim>((float)PointWeight.value * pointWeightSum), true, false);
            else                        iInfo = FEMTree<Dim, float>::template InitializeApproximatePointInterpolationInfo<float, 0> (tree, *samples, ConstraintDual<Dim>(targetValue, (float)PointWeight.value * pointWeightSum), SystemDual<Dim>((float)PointWeight.value * pointWeightSum), true, 1);
            tree.addInterpolationConstraints(constraints, solveDepth, *iInfo);
            profiler.dumpOutput2(comments, "#Set point constraints:");
        }
        
        messageWriter("Leaf Nodes / Active Nodes / Ghost Nodes: %d / %d / %d\n", (int)tree.leaves(), (int)tree.nodes(), (int)tree.ghostNodes());
        messageWriter("Memory Usage: %.3f MB\n", float(MemoryInfo::Usage())/(1<<20));
        
        // Solve the linear system
        {
            profiler.start();
            typename FEMTree<Dim, float>::SolverInfo sInfo;
            sInfo.cgDepth = 0, sInfo.cascadic = true, sInfo.vCycles = 1, sInfo.iters = Iterations.value, sInfo.cgAccuracy = CGSolverAccuracy.value, sInfo.verbose = Verbose.set, sInfo.showResidual = ShowResidual.set, sInfo.showGlobalResidual = SHOW_GLOBAL_RESIDUAL_NONE, sInfo.sliceBlockSize = 1;
            sInfo.baseDepth = BaseMGSolverDepth.value, sInfo.baseVCycles = BaseMGSolverVCycles.value;
            typename FEMIntegrator::template System<Sigs, IsotropicUIntPack<Dim, 1>> F({ 0., 1. });
            solution = tree.solveSystem(Sigs(), F, constraints, solveDepth, sInfo, iInfo);
            profiler.dumpOutput2(comments, "# Linear system solved:");
            if (iInfo) delete iInfo, iInfo = NULL;
        }
    }
    
    {
        profiler.start();
        double valueSum = 0, weightSum = 0;
        typename FEMTree<Dim, float>::template MultiThreadedEvaluator<Sigs, 0> evaluator(&tree, solution);
#pragma omp parallel for reduction(+ : valueSum, weightSum)
        for (int j = 0; j < samples->size(); j++)
        {
            ProjectiveData<Point<float, Dim>, float>& sample = (*samples)[j].sample;
            float w = sample.weight;
            if (w > 0) weightSum += w, valueSum += evaluator.values(sample.data / sample.weight, omp_get_thread_num(), (*samples)[j].node)[0] * w;
        }
        isoValue = (float)(valueSum / weightSum);
        if (DataX.value <= 0 || (!Colors.set && !Normals.set)) delete samples, samples = NULL;
        profiler.dumpOutput("Got average:");
        messageWriter("Iso-Value: %e = %g / %g\n", isoValue, valueSum, weightSum);
    }
    if (Tree.set)
    {
        FILE* fp = fopen(Tree.value, "wb");
        if (!fp) ERROR_OUT("Failed to open file for writing: %s", Tree.value);
        FEMTree<Dim, float>::WriteParameter(fp);
        DenseNodeData<float, Sigs>::WriteSignatures(fp);
        tree.write(fp);
        solution.write(fp);
        fclose(fp);
    }
    
    if (Grid.set)
    {
        int res = 0;
        profiler.start();
        Pointer(float) values = tree.template regularGridEvaluate<true>(solution, res, -1, PrimalGrid.set);
        int resolution = 1;
        for (int d = 0; d < Dim; d++) { resolution *= res; }
#pragma omp parallel for
        for (int i = 0; i < resolution; i++) { values[i] -= isoValue; }
        profiler.dumpOutput("Got grid:");
        WriteGrid<DIMENSION>(values, res, Grid.value);
        DeletePointer(values);
        if (Verbose.set)
        {
            printf("Transform:\n");
            for (int i = 0; i < Dim + 1; i++)
            {
                printf("\t");
                for (int j = 0; j < Dim + 1; j++) printf(" %f", iXForm(j,i));
                printf("\n");
            }
        }
    }
    
    if (OutputMesh.set)
    {
        if (Normals.set)
        {
            if (Density.set)
            {
                typedef PlyVertexWithData<float, Dim, MultiPointStreamData<float, PointStreamNormal<float, Dim>, PointStreamValue<float>, AdditionalPointSampleData>> Vertex;
                std::function<void (Vertex&, Point<float, Dim>, float, TotalPointSampleData)> SetVertex = [](Vertex& v, Point<float, Dim> p, float w, TotalPointSampleData d){ v.point = p, std::get<0>(v.data.data) = std::get<0>(d.data), std::get<1>(v.data.data).data = w, std::get<2>(v.data.data) = std::get<1>(d.data); };
                ExtractMesh<Vertex>(UIntPack<FEMSigs ...>(), std::tuple<SampleData ...>(), tree, solution, isoValue, samples, sampleData, density, SetVertex, comments, iXForm);
            }
            else
            {
                typedef PlyVertexWithData<float, Dim, MultiPointStreamData<float, PointStreamNormal<float, Dim>, AdditionalPointSampleData>> Vertex;
                std::function<void (Vertex&, Point<float, Dim>, float, TotalPointSampleData)> SetVertex = [](Vertex& v, Point<float, Dim> p, float w, TotalPointSampleData d){ v.point = p, std::get<0>(v.data.data) = std::get<0>(d.data), std::get<1>(v.data.data) = std::get<1>(d.data); };
                ExtractMesh<Vertex>(UIntPack<FEMSigs ...>(), std::tuple<SampleData ...>(), tree, solution, isoValue, samples, sampleData, density, SetVertex, comments, iXForm);
            }
        }
        else
        {
            if (Density.set)
            {
                typedef PlyVertexWithData<float, Dim, MultiPointStreamData<float, PointStreamValue<float>, AdditionalPointSampleData>> Vertex;
                std::function<void (Vertex&, Point<float, Dim>, float, TotalPointSampleData)> SetVertex = [](Vertex& v, Point<float, Dim> p, float w, TotalPointSampleData d){ v.point = p, std::get<0>(v.data.data).data = w, std::get<1>(v.data.data) = std::get<1>(d.data); };
                ExtractMesh<Vertex>(UIntPack<FEMSigs ...>(), std::tuple<SampleData ...>(), tree, solution, isoValue, samples, sampleData, density, SetVertex, comments, iXForm);
            }
            else
            {
                typedef PlyVertexWithData<float, Dim, MultiPointStreamData<float, AdditionalPointSampleData>> Vertex;
                std::function<void (Vertex&, Point<float, Dim>, float, TotalPointSampleData)> SetVertex = [](Vertex& v, Point<float, Dim> p, float w, TotalPointSampleData d){ v.point = p, std::get<0>(v.data.data) = std::get<1>(d.data); };
                ExtractMesh<Vertex>(UIntPack<FEMSigs ...>(), std::tuple<SampleData ...>(), tree, solution, isoValue, samples, sampleData, density, SetVertex, comments, iXForm);
            }
        }
        if (sampleData){ delete sampleData; sampleData = NULL; }
    }
    if (density) delete density, density = NULL;
    messageWriter(comments, "#          Total Solve: %9.1f (s), %9.1f (MB)\n", Time()-startTime, FEMTree<Dim, float>::MaxMemoryUsage());
}

template<unsigned int Dim, typename ... SampleData>
void Execute(int argc, char* argv[])
{
    switch (BoundaryType.value)
    {
        case BOUNDARY_FREE+1:
        {
            switch(BSplineDegree.value)
            {
                case 1: return Execute<SampleData ...>(argc, argv, IsotropicUIntPack<Dim, FEMDegreeAndBType<1, BOUNDARY_FREE>::Signature>());
                case 2: return Execute<SampleData ...>(argc, argv, IsotropicUIntPack<Dim, FEMDegreeAndBType<2, BOUNDARY_FREE>::Signature>());
                    //                case 3: return Execute<float, SampleData ...>(argc, argv, IsotropicUIntPack<Dim, FEMDegreeAndBType<3, BOUNDARY_FREE>::Signature>());
                    //                case 4: return Execute<float, SampleData ...>(argc, argv, IsotropicUIntPack<Dim, FEMDegreeAndBType<4, BOUNDARY_FREE>::Signature>());
                default: ERROR_OUT("Only B-Splines of degree 1 - 2 are supported");
            }
        }
        case BOUNDARY_NEUMANN+1:
        {
            switch(BSplineDegree.value)
            {
                case 1: return Execute<SampleData ...>(argc, argv, IsotropicUIntPack<Dim, FEMDegreeAndBType<1, BOUNDARY_NEUMANN>::Signature>());
                case 2: return Execute<SampleData ...>(argc, argv, IsotropicUIntPack<Dim, FEMDegreeAndBType<2, BOUNDARY_NEUMANN>::Signature>());
                    //                case 3: return Execute<float, SampleData ...>(argc, argv, IsotropicUIntPack<Dim, FEMDegreeAndBType<3, BOUNDARY_NEUMANN>::Signature>());
                    //                case 4: return Execute<float, SampleData ...>(argc, argv, IsotropicUIntPack<Dim, FEMDegreeAndBType<4, BOUNDARY_NEUMANN>::Signature>());
                default: ERROR_OUT("Only B-Splines of degree 1 - 2 are supported");
            }
        }
        case BOUNDARY_DIRICHLET+1:
        {
            switch(BSplineDegree.value)
            {
                case 1: return Execute<SampleData ...>(argc, argv, IsotropicUIntPack<Dim, FEMDegreeAndBType<1, BOUNDARY_DIRICHLET>::Signature>());
                case 2: return Execute<SampleData ...>(argc, argv, IsotropicUIntPack<Dim, FEMDegreeAndBType<2, BOUNDARY_DIRICHLET>::Signature>());
                    //            case 3: return Execute<float, SampleData ...>(argc, argv, IsotropicUIntPack<Dim, FEMDegreeAndBType<3, BOUNDARY_DIRICHLET>::Signature>());
                    //            case 4: return Execute<float, SampleData ...>(argc, argv, IsotropicUIntPack<Dim, FEMDegreeAndBType<4, BOUNDARY_DIRICHLET>::Signature>());
                default: ERROR_OUT("Only B-Splines of degree 1 - 2 are supported");
            }
        }
        default: ERROR_OUT("Not a valid boundary type: %d", BoundaryType.value);
    }
}

int main(int argc, char* argv[])
{
    Timer timer;
#ifdef ARRAY_DEBUG
    WARN("Array debugging enabled");
#endif // ARRAY_DEBUG
    
    cmdLineParse(argc-1, &argv[1], params);
    if (MaxMemoryGB.value > 0) SetPeakMemoryMB(MaxMemoryGB.value<<10);
    omp_set_num_threads(Threads.value> 1 ? Threads.value : 1);
    messageWriter.echoSTDOUT = Verbose.set;
    
    if (!InputPoints.set)
    {
        ShowUsage(argv[0]);
        return 0;
    }
    if (DataX.value <= 0) Normals.set = Colors.set = false;
    if (BaseMGSolverDepth.value>FullDepth.value)
    {
        if (BaseMGSolverDepth.set) WARN("Base depth must be smaller than full depth: %d <= %d", BaseMGSolverDepth.value, FullDepth.value);
        BaseMGSolverDepth.value = FullDepth.value;
    }
    
    if (!PointWeight.set) PointWeight.value = DefaultPointWeightMultiplier*BSplineDegree.value;
    if (Colors.set) Execute<DIMENSION, PointStreamColor<float>>(argc, argv);
    else            Execute<DIMENSION>(argc, argv);
    
    if (Performance.set)
    {
        printf("Time (Wall/CPU): %.2f / %.2f\n", timer.wallTime(), timer.cpuTime());
        printf("Peak Memory (MB): %d\n", MemoryInfo::PeakMemoryUsageMB());
    }
    return EXIT_SUCCESS;
}
