#pragma once

#include <map>
// #include "CellHelpers.h"
// #include "SPHKernelFuncs.h"

namespace zeno{
struct PBF {
//params
public:    
    int numSubsteps = 5;
    float dt= 1.0 / 20.0;
    float pRadius = 3.0;
    vec3f bounds_max{40.0, 40.0, 40.0};
    vec3f bounds_min{0,0,0};
    vec3f gravity{0, -10.0, 0};

    float mass = 1.0;
    float rho0 = 1.0;
    float h = 1.1;
    float neighborSearchRadius = h * 1.05;
    float coeffDq = 0.3;
    float coeffK = 0.001;
    float lambdaEpsilon = 100.0; // to prevent the singularity

private:
    void preSolve();
    void solve();
    void postSolve();

    void computeLambda();
    void computeDpos();


    //physical fields
    int numParticles;
    std::vector<vec3f> pos;
    std::vector<vec3f> oldPos;
    std::vector<vec3f> vel;
    std::vector<float> lambda;
    std::vector<vec3f> dpos;

    // std::shared_ptr<zeno::PrimitiveObject> prim;

    //helpers
    void boundaryHandling(vec3f &p);
    inline float computeScorr(const vec3f& distVec, float coeffDq, float coeffK, float h);

    //data for cells
    inline vec3i getCellXYZ(const vec3f& p);
    inline int getCellID(const vec3f& p);
    inline int getCellHash(int i, int j, int k);
    inline bool isInBound(const vec3i& cell);
    inline int cellXYZ2ID(const vec3i& xyz);
    inline vec3i cellID2XYZ(int i);
    std::array<int, 3> numCellXYZ;
    int numCell;
    float dx; //cell size
    float dxInv; 
    void initCellData();
    void initNeighborList();
    struct Cell
    {
        int x,y,z;
        std::vector<int> parInCell; 
    };
    std::map<int, Cell>  cell;

    //neighborList
    std::vector<std::vector<int>> neighborList;
    void neighborSearch();

public:
    void setParams()
    {
        // //用户自定义参数
        // dt = get_input<zeno::NumericObject>("dt")->get<float>();
        // pRadius = get_input<zeno::NumericObject>("particle_radius")->get<float>();
        // bounds_min = get_input<zeno::NumericObject>("bounds_min")->get<vec3f>();
        // bounds_max = get_input<zeno::NumericObject>("bounds_max")->get<vec3f>();
        // gravity = get_input<zeno::NumericObject>("gravity")->get<vec3f>();
        // rho0 = get_input<zeno::NumericObject>("rho0")->get<float>();
        // lambdaEpsilon = get_input<zeno::NumericObject>("lambdaEpsilon")->get<float>();
        // coeffDq = get_input<zeno::NumericObject>("coeffDq")->get<float>();
        // coeffK = get_input<zeno::NumericObject>("coeffK")->get<float>();

        // dx = get_input<zeno::NumericObject>("dx")->get<float>();
        
        // //可以推导出来的参数
        // // auto diam = pRadius*2;
        // // mass = 0.8 * diam*diam*diam * rho0;
        // // h = 4* pRadius;
        // neighborSearchRadius = h;
    }

    void initCube()
    {
        int numParticles = 10000;
        vec3f initPos{10.0,10.0,10.0};
        int cubeSize = 20;
        float spacing = 1;

        // auto &pos = prim->verts;

        std::cout<<"type:"<<typeid(pos).name()<<"\n";
        pos.resize(numParticles);
        int num_per_row = (int) (cubeSize / spacing) + 1; //21
        int num_per_floor = num_per_row * num_per_row; //21 * 21 =441
        for (size_t i = 0; i < numParticles; i++)
        {
            int floor = i / (num_per_floor);
            int row = (i % num_per_floor) / num_per_row ;
            int col = (i % num_per_floor) % num_per_row ;
            pos[i] = vec3f(col*spacing, floor*spacing, row*spacing) + initPos;
        }
    }

    // virtual void apply() override{
    void apply(){
        // auto &pos = prim->verts;

        static bool firstTime = true;
        if(firstTime == true)
        {
            firstTime = false;
            setParams();
            // numParticles = prim->verts.size();
            numParticles = pos.size();

            //fields
            oldPos.resize(numParticles);
            vel.resize(numParticles);
            lambda.resize(numParticles);
            dpos.resize(numParticles);

            initCellData();
            initNeighborList(); 
        }

        preSolve();
        neighborSearch();//grid-baed neighborSearch
        for (size_t i = 0; i < numSubsteps; i++)
            solve(); 
        postSolve();  

    }
};







inline float kernelPoly6(float dist, float h=0.1)
{
    float coeff = 315.0 / 64.0 / 3.14159265358979323846;
    float res = 0.0;
    if(dist > 0 && dist < h)
    {
        float x = (h * h - dist * dist) / (h * h * h);
        res = coeff * x * x * x;
    }
    return res;
}

//SPH kernel gradient
inline zeno::vec3f kernelSpikyGradient(const zeno::vec3f& r, float h=0.1)
{
    float coeff = -45.0 / 3.14159265358979323846;
    zeno::vec3f res{0.0, 0.0, 0.0};
    float dist = length(r);
    if (dist > 0 && dist < h)
    {
        float x = (h - dist) / (h * h * h);
        float factor = coeff * x * x;
        res = r * factor / dist;
    }
    return res;
}








void PBF::preSolve()
{
    // auto &pos = prim->verts;
    for (int i = 0; i < numParticles; i++)
        oldPos[i] = pos[i];

    //update the pos
    for (int i = 0; i < numParticles; i++)
    {
        vec3f tempVel = vel[i];
        tempVel += gravity * dt;
        pos[i] += tempVel * dt;
        boundaryHandling(pos[i]);
    }
}


void PBF::boundaryHandling(vec3f & p)
{
    vec3f bmin = bounds_min + pRadius;
    vec3f bmax = bounds_max - pRadius;

    for (size_t dim = 0; dim < 3; dim++)
    {
        float r = ((float) rand() / (RAND_MAX));
        if (p[dim] <= bmin[dim])
            p[dim] = bmin[dim] + 1e-5 * r;
        else if (p[dim]>= bmax[dim])
            p[dim] = bmax[dim] - 1e-5 * r;
    }
}

void PBF::solve()
{
    // auto &pos = prim->verts;

    computeLambda();

    computeDpos();

    //apply the dpos to the pos
    for (size_t i = 0; i < numParticles; i++)
        pos[i] += dpos[i];
}

void PBF::computeLambda()
{
    lambda.clear();
    lambda.resize(numParticles);
    // auto &pos = prim->verts;

    for (size_t i = 0; i < numParticles; i++)
    {
        vec3f gradI{0.0, 0.0, 0.0};
        float sumSqr = 0.0;
        float densityCons = 0.0;

        for (size_t j = 0; j < neighborList[i].size(); j++)
        {
            int pj = neighborList[i][j];
            vec3f distVec = pos[i] - pos[pj];
            vec3f gradJ = kernelSpikyGradient(distVec, h);
            gradI += gradJ;
            sumSqr += dot(gradJ, gradJ);
            densityCons += kernelPoly6(length(distVec), h);
        }
        densityCons = (mass * densityCons / rho0) - 1.0;

        //compute lambda
        sumSqr += dot(gradI, gradI);
        lambda[i] = (-densityCons) / (sumSqr + lambdaEpsilon);
    }
}

void PBF::computeDpos()
{
    dpos.clear();
    dpos.resize(numParticles);
    // auto &pos = prim->verts;

    for (size_t i = 0; i < numParticles; i++)
    {
        vec3f dposI{0.0, 0.0, 0.0};
        for (size_t j = 0; j < neighborList[i].size(); j++)
        {
            int pj = neighborList[i][j];
            vec3f distVec = pos[i] - pos[pj];

            float sCorr = computeScorr(distVec, coeffDq, coeffK, h);
            dposI += (lambda[i] + lambda[pj] + sCorr) * kernelSpikyGradient(distVec, h);
        }
        dposI /= rho0;
        dpos[i] = dposI;
    }
}

//helper for computeDpos()
inline float PBF::computeScorr(const vec3f& distVec, float coeffDq, float coeffK, float h)
{
    float x = kernelPoly6(length(distVec), h) / kernelPoly6(coeffDq * h, h);
    x = x * x;
    x = x * x;
    return (-coeffK) * x;
}


void PBF::postSolve()
{
    // auto &pos = prim->verts;

    for (size_t i = 0; i < numParticles; i++)
        vel[i] = (pos[i] - oldPos[i]) / dt;
}

















//helpers for neighborSearch
inline bool PBF::isInBound(const vec3i& cellXYZ)
{
    return cellXYZ[0] >= 0 && cellXYZ[0] < numCellXYZ[0] &&
           cellXYZ[1] >= 0 && cellXYZ[1] < numCellXYZ[1] &&
           cellXYZ[2] >= 0 && cellXYZ[2] < numCellXYZ[2];
}

inline int PBF::getCellID(const vec3f& p)
{
    vec3i xyz = p*dxInv;
    int numPerRow = numCellXYZ[0];
    int numPerFloor = numCellXYZ[0] * numCellXYZ[1];
    int res = numPerFloor * xyz[2] + numPerRow * xyz[1] + xyz[0];
    return res;
}

inline int PBF::cellXYZ2ID(const vec3i& xyz)
{
    int numPerRow = numCellXYZ[0];
    int numPerFloor = numCellXYZ[0] * numCellXYZ[1];
    int res = numPerFloor * xyz[2] + numPerRow * xyz[1] + xyz[0];
    return res;
}

inline vec3i PBF::cellID2XYZ(int i)
{
    //to calculate the x y z coord of cell
    int numPerRow = numCellXYZ[0];
    int numPerFloor = numCellXYZ[0] * numCellXYZ[1];

    int floor = (i / numPerFloor);
    int row = (i % numPerFloor) / numPerRow;
    int col = (i % numPerFloor) % numPerRow; 
    
    vec3i res{col,row,floor};
    return res;
}

inline vec3i PBF::getCellXYZ(const vec3f& p)
{
    vec3i res{p*dxInv};
    return res;
}

inline int PBF::getCellHash(int i, int j, int k)
{
    int res = ( (73856093 * i) ^ 
                (19349663 * j) ^ 
                (83492791 * k) ) 
                % (2 * numParticles);
    return res;
}

inline void PBF::initNeighborList()
{     
    //prepare neighbor list 
    neighborList.resize(numParticles);
    for (size_t i = 0; i < numParticles; i++)
        neighborList[i].reserve(10);
}

inline void PBF::initCellData()
{
    //calc cell size
    dxInv = 1.0/dx;
    int numX = int((bounds_max[0]-bounds_min[0]) *dxInv) + 1;
    int numY = int((bounds_max[1]-bounds_min[1]) *dxInv) + 1;
    int numZ = int((bounds_max[2]-bounds_min[2]) *dxInv) + 1;
    // numCellXYZ.resize(3);
    numCellXYZ[0] = numX;
    numCellXYZ[1] = numY;
    numCellXYZ[2] = numZ;
    numCell = numX * numY * numZ;


    //prepare cell data
    for (size_t i = 0; i < numCell; i++)
    {
        // //to calculate the x y z coord of cell
        vec3i xyz = cellID2XYZ(i);
        int hash = getCellHash(xyz[0],xyz[1],xyz[2]);

        cell[hash].x = xyz[0];
        cell[hash].y = xyz[1];
        cell[hash].z = xyz[2];
        cell[hash].parInCell.reserve(10); //pre-allocate memory to speed up
    }
}



// This neighborSearch algorithm uses grid-based searching,
// which is simple but can be improved
void PBF::neighborSearch()
{
    // auto &pos = prim->verts;
    //clear parInCell and neighborList
    for (size_t i = 0; i < numCell; i++)
        cell[i].parInCell.clear();
    for (size_t i = 0; i < numParticles; i++)
        neighborList[i].clear();
    
    //update the parInCell list
    vec3i cellXYZ;
    int cellID;
    for (size_t i = 0; i < numParticles; i++) // i is the particle ID
    {
        cellID = getCellID(pos[i]);
        // int hash = getCellHash(cellXYZ[0], cellXYZ[1], cellXYZ[2]);
        cell[cellID].parInCell.push_back(i);
    }

    //update the neighborList
    for (size_t i = 0; i < numParticles; i++)
    {
        cellXYZ = getCellXYZ(pos[i]);
        for (int off_x = -1; off_x < 2; off_x++)
            for (int off_y = -1; off_y < 2; off_y++)
                for (int off_z = -1; off_z < 2; off_z++)
                {
                    vec3i off{off_x, off_y, off_z};
                    vec3i toCheckXYZ = cellXYZ + off;
                    int toCheck = cellXYZ2ID(toCheckXYZ);
                    if (isInBound(toCheckXYZ))
                    {
                        Cell theCell = cell[toCheck];
                        std::vector<int> parInTheCell = theCell.parInCell;

                        for (int j = 0; j < parInTheCell.size(); j++)
                        {
                            int p = parInTheCell[j];
                            if(p!=i && length(pos[i] - pos[p]) < neighborSearchRadius)
                            {
                                neighborList[i].push_back(p);
                            }
                        }
                    }
                }  
    }
}












}