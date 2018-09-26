#include <cassert>
#include <iostream>

#include "KDTree.h"


struct Point
{
    float x, y, z;
};

class PointWrapper
{
    const Point * data;
public:

    PointWrapper(const void * data) :
        data(reinterpret_cast<const Point *>(data))
    {
    }

    float Get(int dim)
    {
        switch (dim)
        {
            case 0:
                return data->x;
                break;
            case 1:
                return data->y;
                break;
            case 2:
                return data->z;
                break;
            default:
                assert(false && "Error accessing point coordinates");
        }

        return 0.f;
    }
};

int main()
{
    std::vector<Point> pointVector;
    pointVector.resize(1000);

    for (unsigned int i = 0; i < pointVector.size(); ++i)
    {
        pointVector[i].x = (float)i;
        pointVector[i].y = (float)i;
        pointVector[i].z = (float)i;
    }

    gbi::KDTree<PointWrapper, 3> kdTree;

    int rebalanceCount = 0;

    for (auto & point : pointVector)
    {
        kdTree.Insert(&point);

        bool rebalance = true;
        while (rebalance)
        {
            rebalance = kdTree.RebalanceIteration();
            ++rebalanceCount;
        }
    }

    std::cout << "Rebalance count: " << rebalanceCount << std::endl;

    return 0;
}
