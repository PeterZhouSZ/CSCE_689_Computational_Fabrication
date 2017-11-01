//Computational Fabrication Assignment #1
// Created by David Levin 2014
// Modified by Shinjiro Sueda 2016

#include <iostream>
#include <vector>
#include <omp.h>

#include "../include/CompFab.h"
#include "../include/Mesh.h"
using namespace std;

//Triangle list (global)
typedef std::vector<CompFab::Triangle> TriangleList;

TriangleList g_triangleList;
CompFab::VoxelGrid *g_voxelGrid;

// This test function is adapted from Moller-Trumbore intersection algorithm. 
// See https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm 
bool rayTriangleIntersects(CompFab::Triangle &triangle, CompFab::Vec3 dir, CompFab::Vec3 pos){

    CompFab::Vec3 e1 = triangle.m_v2 - triangle.m_v1;
    CompFab::Vec3 e2 = triangle.m_v3 - triangle.m_v1;

    // Calculate planes normal vector
    //cross product
    CompFab::Vec3 pvec = dir % e2; 
    
    //dot product
    float det = e1 * pvec; 

    // Ray is parallel to plane
    if (det <1e-8 && det > -1e-8){
        return false;
    }

    float inv_det = 1/det;

    // Distance from v1 to ray pos
    CompFab::Vec3 tvec = pos - triangle.m_v1;
    float u = tvec * pvec * inv_det;
    if (u < 0 || u > 1){
        return false;
    }

    CompFab::Vec3 qvec = tvec % e1;
    float v = dir * qvec * inv_det;
    if (v<0 || u+v>1){
        return false;
    }

    float t = (e2 * qvec) * inv_det;
    if (t > 1e-8) return true;
    return false;
}



bool loadMesh(char *filename, unsigned int num)
{
    g_triangleList.clear();
    
    Mesh *tempMesh = new Mesh(filename, true);
    
    CompFab::Vec3 v1, v2, v3;
    
    //copy triangles to global list
    for(unsigned int tri =0; tri<tempMesh->t.size(); ++tri)
    {
        v1 = tempMesh->v[tempMesh->t[tri][0]];
        v2 = tempMesh->v[tempMesh->t[tri][1]];
        v3 = tempMesh->v[tempMesh->t[tri][2]];
        g_triangleList.push_back(CompFab::Triangle(v1,v2,v3));
    }
    
    //Create Voxel Grid
    CompFab::Vec3 bbMax, bbMin;
    BBox(*tempMesh, bbMin, bbMax);
    
    //Build Voxel Grid
    double bbX = bbMax[0] - bbMin[0];
    double bbY = bbMax[1] - bbMin[1];
    double bbZ = bbMax[2] - bbMin[2];
    double spacing;
    
    if(bbX > bbY && bbX > bbZ)
    {
        spacing = bbX/(num-1);
    } else if(bbY > bbX && bbY > bbZ) {
        spacing = bbY/(num-1);
    } else {
        spacing = bbZ/(num-1);
    }
    
    g_voxelGrid = new CompFab::VoxelGrid(bbMin, num, num, num, spacing);
    
    delete tempMesh;
    
    return true;
    
}

void saveVoxelsToObj(const char * outfile)
{
    Mesh box;
    Mesh mout;
    int nx = g_voxelGrid->m_numX;
    int ny = g_voxelGrid->m_numY;
    int nz = g_voxelGrid->m_numZ;
    double spacing = g_voxelGrid->m_spacing;
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                if(!g_voxelGrid->isInside(ii,jj,kk)){
                    continue;
                }
                CompFab::Vec3 coord(0.5f + ((double)ii)*spacing, 0.5f + ((double)jj)*spacing, 0.5f+((double)kk)*spacing);
                CompFab::Vec3 box0 = coord - hspacing;
                CompFab::Vec3 box1 = coord + hspacing;
                makeCube(box, box0, box1);
                mout.append(box);
            }
        }
    }
    
    mout.save_obj(outfile);
}


int main(int argc, char **argv)
{
    unsigned int num = 16; //number of voxels (e.g. 16x16x16)
    
    if(argc < 5)
    {
        std::cout<<"Usage: voxelizer input.obj output.obj num isRobust(0 or 1)\n";
        exit(0);
    }
    num = atoi(argv[3]);
    
    // The loadMesh() function loads the mesh and then creates:
    // - g_triangleList: The list of triangles in the input OBJ mesh
    // - g_voxelGrid: The VoxelGrid object, with all voxels marked as unoccupied
    std::cout<<"Load Mesh : "<<argv[1]<<"\n";
    loadMesh(argv[1], num);
    
    std::cout<<"Voxelizing into "<<num<<"x"<<num<<"x"<<num<<"\n";

    // Whether to use multiple rays to intersect: 0 for 1 ray; 1 for 6 rays
    int multiRays = atoi(argv[4]);

    
    // Below, write a triple for-loop to go through all the voxels in X, Y, Z.
    //   g_voxelGrid->m_numX is the number of voxels in the X direction.
    //   g_voxelGrid->m_numY is the number of voxels in the Y direction.
    //   g_voxelGrid->m_numZ is the number of voxels in the Z direction.
    //   g_voxelGrid->m_spacing is the size of each voxel.
    //
    // Inside the triple for-loop, check if the voxel is inside or outside the
    // mesh. Use the g_voxelGrid->isInside(...) method to set whether that voxel
    // is inside or outside. E.g.,
    //   g_voxelGrid->isInside(0,0,0) = false;

    // <TRIPLE FOR-LOOP HERE>

    double spacing = g_voxelGrid->m_spacing;
    CompFab::Vec3 cornerPos = g_voxelGrid->m_corner;

    // a vector for multiple directions of 6 rays    
    std::vector<CompFab::Vec3> dirVec;

    // Specify a direction or multiple directions

    if(!multiRays){
        dirVec.push_back(CompFab::Vec3(1.0, 0.0, 0.0));
    }else{
        dirVec.push_back(CompFab::Vec3(1.0, 0.0, 0.0));
        dirVec.push_back(CompFab::Vec3(-1.0, 0.0, 0.0));
        dirVec.push_back(CompFab::Vec3(0.0, 1.0, 0.0));
        dirVec.push_back(CompFab::Vec3(0.0, -1.0, 0.0));
        dirVec.push_back(CompFab::Vec3(0.0, 0.0, 1.0));
        dirVec.push_back(CompFab::Vec3(0.0, 0.0, -1.0));
    }

    // Triple for loop

    // parallize each voxel 

    #pragma omp parallel for collapse(3)
    for(int i=0; i < g_voxelGrid->m_numX; i++){
        for(int j=0; j < g_voxelGrid->m_numY; j++){
            for(int k=0; k < g_voxelGrid->m_numZ; k++){

                int nIntersections = 0;
                
                //Position of the current voxel
                CompFab::Vec3 voxPos(cornerPos.m_x + i * spacing, 
                                     cornerPos.m_y + j * spacing, 
                                     cornerPos.m_z + k * spacing);

            

                if (!multiRays){
                    // loop through all the triangles 
                    for (int m = 0; m < g_triangleList.size(); m++){
                        if (rayTriangleIntersects(g_triangleList[m], dirVec[0], voxPos)){
                            nIntersections += 1;
                        }
                    }

                    if (nIntersections % 2!= 0){
                        g_voxelGrid->isInside(i,j,k) = true;
                    }else{
                        g_voxelGrid->isInside(i,j,k) = false;
                    }
                    
                }else{
                    // multirays
                    int vote = 0;

                    for (int it = 0; it < 6; it ++){
                        int nIntersections = 0;

                        for(int m = 0; m < g_triangleList.size(); m++){
                            if(rayTriangleIntersects(g_triangleList[m], dirVec[it], voxPos)){
                                nIntersections += 1;
                            }
                        }

                        if(nIntersections % 2 != 0){
                            vote += 1;
                        }
                    }

                    // Mark a voxel as being inside if 4 or more tests return inside
                    if (vote >= 4){
                        g_voxelGrid->isInside(i,j,k) = true;
                        
                    } else{
                        g_voxelGrid->isInside(i,j,k) = false;
                    }
                }

            }
        }
    }

    //Write out voxel data as obj
    saveVoxelsToObj(argv[2]);
    
    delete g_voxelGrid;
}