#ifndef BRANCHTREE_H
#define BRANCHTREE_H

#include "v3d_interface.h"

using namespace std;

struct Branch{
    V3DLONG headPointIndex,endPointIndex;
    Branch* parent;
    int level;
    float length, distance;

    Branch():level(-1),headPointIndex(-1),endPointIndex(-1)
    {
        parent=0;
    }
    ~Branch(){}

    bool get_pointsIndex_of_branch(vector<V3DLONG> &points,NeuronTree &nt);
    bool get_r_pointsIndex_of_branch(vector<V3DLONG> &r_points, NeuronTree &nt);
    void getBranchFeature(unsigned char *pdata, V3DLONG* sz, NeuronTree& nt);
};

struct BranchTree
{
    vector<Branch> branchs;
    NeuronTree nt;
    BranchTree():nt()
    {
        branchs.clear();
    }
    ~BranchTree(){
        branchs.clear();
    }

    bool initialize(NeuronTree t);
    bool get_level_index(vector<int> &level_index,int level);
    int get_max_level();

    void getBranchFeature(unsigned char *pdata, V3DLONG* sz);
};

#endif // BRANCHTREE_H
