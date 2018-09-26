#include <array>
#include <cstdint>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

#ifdef ENABLE_GBI_ASSERTS
    #include <cassert>
    #define gbiAssert assert
#else
    #define gbiAssert(x)
#endif

namespace gbi
{

typedef uint32_t UInt;
typedef int32_t Int;
typedef void* PointData;
typedef std::vector<PointData>::size_type SizeT;

template<typename PointWrapper, UInt Dimension>
class KDTree
{
    std::vector<PointData> pointDataVector;
    std::unordered_map<PointData, SizeT> indexedPointData;

    std::unordered_map<SizeT, UInt> nodeFloor;

    std::unordered_map<SizeT, SizeT> lowerIndex;
    std::unordered_map<SizeT, SizeT> upperIndex;
    std::unordered_map<SizeT, SizeT> parentIndex;

    std::unordered_map<PointData, std::array<PointData, Dimension * 2>> boundaries;

    std::unordered_map<SizeT, Int> balanceByIndex;
    std::map<Int, std::map<UInt, std::set<SizeT>>> indexByBalanceAndFloor;

    // Data
    PointData origin;

public:

    // Constructor
    KDTree() :
        origin(nullptr)
    {
    }

private:

    // Methods
    SizeT InsertInternal(PointData data)
    {
        gbiAssert(data != nullptr && "Cannot insert nullptr point");

        pointDataVector.push_back(data);
        SizeT index = pointDataVector.size() - 1;
        indexedPointData[data] = index;

        return index;
    }

    PointData GetLower(PointData pointData)
    {
        gbiAssert(indexedPointData.find(pointData) != indexedPointData.end() && "Trying to find lower of unregistered point");

        SizeT dataIndex = indexedPointData[pointData];

        if (lowerIndex.find(dataIndex) != lowerIndex.end())
            return pointDataVector.at(lowerIndex[dataIndex]);

        return nullptr;
    }

    PointData GetUpper(PointData pointData)
    {
        gbiAssert(indexedPointData.find(pointData) != indexedPointData.end() && "Trying to find upper of unregistered point");

        SizeT dataIndex = indexedPointData[pointData];

        if (upperIndex.find(dataIndex) != upperIndex.end())
            return pointDataVector.at(upperIndex[dataIndex]);

        return nullptr;
    }

    PointData GetParent(PointData pointData)
    {
        gbiAssert(indexedPointData.find(pointData) != indexedPointData.end() && "Trying to find parent of unregistered point");

        PointData parent = nullptr;

        if (parentIndex.find(indexedPointData[pointData]) != parentIndex.end())
        {
            parent = pointDataVector.at(parentIndex[indexedPointData[pointData]]);
        }
        else
        {
            gbiAssert(pointData == origin && "KDTree corruption: point orphan node found, only accepted for origin");
        }

        return parent;
    }

    void SetLower(PointData parent, PointData child)
    {
        gbiAssert(indexedPointData.find(parent) != indexedPointData.end());
        gbiAssert(indexedPointData.find(child) != indexedPointData.end());
        gbiAssert(nodeFloor.find(indexedPointData[parent]) != nodeFloor.end());

        lowerIndex[indexedPointData[parent]] = indexedPointData[child];
        parentIndex[indexedPointData[child]] = indexedPointData[parent];
        nodeFloor[indexedPointData[child]] = nodeFloor[indexedPointData[parent]] + 1;
    }

    void SetUpper(PointData parent, PointData child)
    {
        gbiAssert(indexedPointData.find(parent) != indexedPointData.end());
        gbiAssert(indexedPointData.find(child) != indexedPointData.end());
        gbiAssert(nodeFloor.find(indexedPointData[parent]) != nodeFloor.end());

        upperIndex[indexedPointData[parent]] = indexedPointData[child];
        parentIndex[indexedPointData[child]] = indexedPointData[parent];
        nodeFloor[indexedPointData[child]] = nodeFloor[indexedPointData[parent]] + 1;
    }

    void UpdateBoundaries(const std::vector<PointData> & updateList, PointData point)
    {
        gbiAssert(indexedPointData.find(point) != indexedPointData.end() && "Boundary not found in KDTree");

        PointWrapper boundary(point);
        for (auto pointToUpdateData : updateList)
        {
            gbiAssert(indexedPointData.find(pointToUpdateData) != indexedPointData.end() && "Point to update not found in KDTree");

            PointWrapper pointToUpdate(pointToUpdateData);

            if (boundaries.find(pointToUpdateData) == boundaries.end())
            {
                for (size_t i = 0; i < Dimension; ++i)
                {
                    boundaries[pointToUpdateData][2 * i] = point;
                    boundaries[pointToUpdateData][2 * i + 1] = point;
                }
            }
            else
            {
                for (size_t i = 0; i < Dimension; ++i)
                {
                    if (PointWrapper(point).Get(i) < PointWrapper(boundaries[pointToUpdateData][2 * i]).Get(i))
                    {
                        boundaries[pointToUpdateData][2 * i] = point;
                    }

                    if (PointWrapper(boundaries[pointToUpdateData][2 * i + 1]).Get(i) < PointWrapper(point).Get(i))
                    {
                        boundaries[pointToUpdateData][2 * i + 1] = point;
                    }
                }
            }
        }
    }

    void InternalRemoveBalancePriority(SizeT pointIndex)
    {
        if (balanceByIndex.find(pointIndex) != balanceByIndex.end())
        {
            Int pointBalance = balanceByIndex[pointIndex];

            gbiAssert(nodeFloor.find(pointIndex) != nodeFloor.end());
            gbiAssert(
                indexByBalanceAndFloor.find(pointBalance) != indexByBalanceAndFloor.end() &&
                indexByBalanceAndFloor[pointBalance].find(nodeFloor[pointIndex]) != indexByBalanceAndFloor[pointBalance].end() &&
                indexByBalanceAndFloor[pointBalance][nodeFloor[pointIndex]].find(pointIndex) != indexByBalanceAndFloor[pointBalance][nodeFloor[pointIndex]].end()
            );

            indexByBalanceAndFloor[pointBalance][nodeFloor[pointIndex]].erase(pointIndex);
            if (indexByBalanceAndFloor[pointBalance][nodeFloor[pointIndex]].size() == 0)
            {
                indexByBalanceAndFloor[pointBalance].erase(nodeFloor[pointIndex]);
                if (indexByBalanceAndFloor[pointBalance].size() == 0)
                {
                    indexByBalanceAndFloor.erase(pointBalance);
                }
            }
        }
    }

    void UpdateBalance(PointData pointToUpdate, bool increment)
    {
        gbiAssert(indexedPointData.find(pointToUpdate) != indexedPointData.end() && "Cannot find point when updating balance");
        gbiAssert(nodeFloor.find(indexedPointData[pointToUpdate]) != nodeFloor.end());

        SizeT updateIndex = indexedPointData[pointToUpdate];

        InternalRemoveBalancePriority(updateIndex);

        if (increment)
        {
            ++balanceByIndex[updateIndex];
        }
        else
        {
            --balanceByIndex[updateIndex];
        }

        indexByBalanceAndFloor[balanceByIndex[updateIndex]][nodeFloor[updateIndex]].insert(updateIndex);
    }

    bool IsLeaf(SizeT pointIndex) const
    {
        return (lowerIndex.find(pointIndex) == lowerIndex.end()) && (upperIndex.find(pointIndex) == upperIndex.end());
    }

    void UpdateBoundariesWithPoint(PointData other, std::array<PointData, Dimension * 2> & boundariesSlot, bool & insertedSomething)
    {
        if (other != nullptr)
        {
            gbiAssert(indexedPointData.find(other) != indexedPointData.end());

            for (size_t i = 0; i < Dimension; ++i)
            {
                boundariesSlot[2 * i] = other;
                boundariesSlot[2 * i + 1] = other;
            }

            if (boundaries.find(other) != boundaries.end())
            {
                for (size_t i = 0; i < Dimension; ++i)
                {
                    boundariesSlot[2 * i] =
                        PointWrapper(boundariesSlot[2 * i]).Get(i) < PointWrapper(boundaries[other][2 * i]).Get(i) ?
                        boundariesSlot[2 * i] :
                        boundaries[other][2 * i];

                    boundariesSlot[2 * i + 1] =
                        PointWrapper(boundariesSlot[2 * i + 1]).Get(i) < PointWrapper(boundaries[other][2 * i + 1]).Get(i) ?
                        boundaries[other][2 * i + 1] :
                        boundariesSlot[2 * i + 1];
                }
            }

            insertedSomething = true;
        }
    }

    void UpdateBoundaries(SizeT pointIndex)
    {
        if (boundaries.find(pointDataVector.at(pointIndex)) != boundaries.end())
            boundaries.erase(pointDataVector.at(pointIndex));

        PointData lowerPoint = lowerIndex.find(pointIndex) != lowerIndex.end() ? pointDataVector.at(lowerIndex[pointIndex]) : nullptr;
        PointData upperPoint = upperIndex.find(pointIndex) != upperIndex.end() ? pointDataVector.at(upperIndex[pointIndex]) : nullptr;

        std::array<PointData, Dimension * 2> boundariesSlot;
        bool insertedSomething = false;

        UpdateBoundariesWithPoint(lowerPoint, boundariesSlot, insertedSomething);
        UpdateBoundariesWithPoint(upperPoint, boundariesSlot, insertedSomething);

        if (insertedSomething)
            boundaries[pointDataVector.at(pointIndex)] = boundariesSlot;
    }

    void MoveLastElementTo(SizeT itemIndex)
    {
        SizeT lastIndex = pointDataVector.size() - 1;
        PointData previousDataAtLocation = pointDataVector.at(itemIndex);

        gbiAssert(pointDataVector.size() > itemIndex);
        gbiAssert(indexedPointData.find(pointDataVector.at(itemIndex)) != indexedPointData.end());
        gbiAssert(indexedPointData.find(pointDataVector.at(lastIndex)) != indexedPointData.end());

        if (itemIndex != lastIndex)
        {
            if (parentIndex.find(lastIndex) != parentIndex.end())
            {
                if (lowerIndex.find(parentIndex[lastIndex]) != lowerIndex.end() && lowerIndex[parentIndex[lastIndex]] == lastIndex)
                    lowerIndex[parentIndex[lastIndex]] = itemIndex;
                else if (upperIndex.find(parentIndex[lastIndex]) != upperIndex.end() && upperIndex[parentIndex[lastIndex]] == lastIndex)
                    upperIndex[parentIndex[lastIndex]] = itemIndex;
                else
                    gbiAssert(false && "Has parent that doesn't refer item as child");
            }

            if (balanceByIndex.find(lastIndex) != balanceByIndex.end())
            {
                balanceByIndex[itemIndex] = balanceByIndex[lastIndex];
                indexByBalanceAndFloor[balanceByIndex[lastIndex]][nodeFloor[itemIndex]].insert(itemIndex);
            }
            else
            {
                balanceByIndex.erase(itemIndex);
            }

            if (lowerIndex.find(lastIndex) != lowerIndex.end())
            {
                lowerIndex[itemIndex] = lowerIndex[lastIndex];
            }
            else
            {
                lowerIndex.erase(itemIndex);
            }

            if (upperIndex.find(lastIndex) != upperIndex.end())
            {
                upperIndex[itemIndex] = upperIndex[lastIndex];
            }
            else
            {
                upperIndex.erase(itemIndex);
            }

            if (parentIndex.find(lastIndex) != parentIndex.end())
            {
                parentIndex[itemIndex] = parentIndex[lastIndex];
            }
            else
            {
                parentIndex.erase(itemIndex);
            }

            //if (boundaries.find(pointDataVector.at(lastIndex)) != boundaries.end())
            //{
            //    boundaries[pointDataVector.at(itemIndex)] = boundaries[pointDataVector.at(lastIndex)];
            //}
            //else
            //{
            boundaries.erase(pointDataVector.at(itemIndex));
            //}

            gbiAssert(nodeFloor.find(lastIndex) != nodeFloor.end());

            nodeFloor[itemIndex] = nodeFloor[lastIndex];
            pointDataVector[itemIndex] = pointDataVector[lastIndex];
            indexedPointData[pointDataVector.at(lastIndex)] = itemIndex;
        }

        indexByBalanceAndFloor[balanceByIndex[lastIndex]][nodeFloor[lastIndex]].erase(lastIndex);
        if (indexByBalanceAndFloor[balanceByIndex[lastIndex]][nodeFloor[lastIndex]].size() == 0)
        {
            indexByBalanceAndFloor[balanceByIndex[lastIndex]].erase(nodeFloor[lastIndex]);
            if (indexByBalanceAndFloor[balanceByIndex[lastIndex]].size() == 0)
            {
                indexByBalanceAndFloor.erase(balanceByIndex[lastIndex]);
            }
        }
        balanceByIndex.erase(lastIndex);

        pointDataVector.resize(pointDataVector.size() - 1);
        indexedPointData.erase(previousDataAtLocation);
        nodeFloor.erase(lastIndex);
        lowerIndex.erase(lastIndex);
        upperIndex.erase(lastIndex);
        parentIndex.erase(lastIndex);
    }

    PointData RemoveLeaf(SizeT pointIndex, bool & incrementNewLeafBalance)
    {
        PointData parent = nullptr;

        gbiAssert(lowerIndex.find(pointIndex) == lowerIndex.end());
        gbiAssert(upperIndex.find(pointIndex) == upperIndex.end());
        gbiAssert(boundaries.find(pointDataVector.at(pointIndex)) == boundaries.end());

        // Remove association in parent
        if (parentIndex.find(pointIndex) != parentIndex.end())
        {
            gbiAssert((
                (lowerIndex.find(parentIndex[pointIndex]) != lowerIndex.end() && lowerIndex[parentIndex[pointIndex]] == pointIndex) ||
                (upperIndex.find(parentIndex[pointIndex]) != upperIndex.end() && upperIndex[parentIndex[pointIndex]] == pointIndex)) &&
                "Inconsistent parent-child reference"
            );

            if (lowerIndex.find(parentIndex[pointIndex]) != lowerIndex.end() && lowerIndex[parentIndex[pointIndex]] == pointIndex)
            {
                lowerIndex.erase(parentIndex[pointIndex]);
                incrementNewLeafBalance = true;
            }
            else
            {
                upperIndex.erase(parentIndex[pointIndex]);
                incrementNewLeafBalance = false;
            }

            parent = pointDataVector[parentIndex[pointIndex]];
        }

        if (balanceByIndex.find(pointIndex) != balanceByIndex.end())
        {
            gbiAssert(indexByBalanceAndFloor.find(balanceByIndex[pointIndex]) != indexByBalanceAndFloor.end());
            gbiAssert(nodeFloor.find(pointIndex) != nodeFloor.end());
            gbiAssert(indexByBalanceAndFloor[balanceByIndex[pointIndex]].find(nodeFloor[pointIndex]) != indexByBalanceAndFloor[balanceByIndex[pointIndex]].end());
            gbiAssert(indexByBalanceAndFloor[balanceByIndex[pointIndex]][nodeFloor[pointIndex]].find(pointIndex) != indexByBalanceAndFloor[balanceByIndex[pointIndex]][nodeFloor[pointIndex]].end());

            indexByBalanceAndFloor[balanceByIndex[pointIndex]][nodeFloor[pointIndex]].erase(pointIndex);
            if (indexByBalanceAndFloor[balanceByIndex[pointIndex]][nodeFloor[pointIndex]].size() == 0)
            {
                indexByBalanceAndFloor[balanceByIndex[pointIndex]].erase(nodeFloor[pointIndex]);
                if (indexByBalanceAndFloor[balanceByIndex[pointIndex]].size() == 0)
                {
                    indexByBalanceAndFloor.erase(balanceByIndex[pointIndex]);
                }
            }
            balanceByIndex.erase(pointIndex);
        }

        nodeFloor.erase(pointIndex);
        parentIndex.erase(pointIndex);

        MoveLastElementTo(pointIndex);

        return parent;
    }

    // Same with different input variables
    bool HasLowerCoordinate(UInt dim, PointData var1, SizeT var2)
    {
        gbiAssert(indexedPointData.find(var1) != indexedPointData.end());

        return HasLowerCoordinate(dim, indexedPointData[var1], var2);
    }

    bool HasLowerCoordinate(UInt dim, SizeT var1, SizeT var2)
    {
        PointWrapper p1(pointDataVector.at(var1));
        PointWrapper p2(pointDataVector.at(var2));

        return p1.Get(dim) < p2.Get(dim);
    }

    bool HasHigherCoordinate(UInt dim, SizeT var1, SizeT var2)
    {
        return !HasLowerCoordinate(dim, var1, var2);
    }

    // Same� with different input variables
    bool HasHigherCoordinate(UInt dim, PointData var1, SizeT var2)
    {
        gbiAssert(indexedPointData.find(var1) != indexedPointData.end());

        return HasHigherCoordinate(dim, indexedPointData[var1], var2);
    }

    SizeT GetBestReplacementCandidate(SizeT index)
    {
        SizeT result;

        gbiAssert(!IsLeaf(index));
        gbiAssert(boundaries.find(pointDataVector.at(index)) != boundaries.end());
        gbiAssert(balanceByIndex.find(index) != balanceByIndex.end());
        gbiAssert(nodeFloor.find(index) != nodeFloor.end());

        UInt dim = nodeFloor[index] % Dimension;

        if (balanceByIndex[index] < 0)
        {
            gbiAssert(lowerIndex.find(index) != lowerIndex.end());

            result = lowerIndex[index];

            if (boundaries.find(pointDataVector.at(lowerIndex[index])) != boundaries.end() && HasHigherCoordinate(dim, boundaries[pointDataVector.at(lowerIndex[index])][dim * 2], result))
            {
                gbiAssert(indexedPointData.find(boundaries[pointDataVector.at(lowerIndex[index])][dim * 2]) != indexedPointData.end());

                result = indexedPointData[boundaries[pointDataVector.at(lowerIndex[index])][dim * 2]];
            }
        }
        else
        {
            gbiAssert(upperIndex.find(index) != upperIndex.end());

            result = upperIndex[index];

            if (boundaries.find(pointDataVector.at(upperIndex[index])) != boundaries.end() && HasLowerCoordinate(dim, boundaries[pointDataVector.at(upperIndex[index])][dim * 2], result))
            {
                gbiAssert(indexedPointData.find(boundaries[pointDataVector.at(upperIndex[index])][dim * 2]) != indexedPointData.end());

                result = indexedPointData[boundaries[pointDataVector.at(upperIndex[index])][dim * 2]];
            }
        }

        return result;
    }

    void SwapNodes(SizeT dst, SizeT src)
    {
        gbiAssert(pointDataVector.size() > dst);
        gbiAssert(pointDataVector.size() > src);

        if (pointDataVector[dst] == origin)
            origin = pointDataVector[src];

        boundaries.erase(pointDataVector[dst]);
        boundaries.erase(pointDataVector[src]);

        std::swap(pointDataVector[dst], pointDataVector[src]);
        std::swap(indexedPointData[pointDataVector[dst]], indexedPointData[pointDataVector[src]]);
    }

public:

    void Insert(PointData point)
    {
        if (point != nullptr)
        {
            PointData parent = nullptr;
            PointData current = origin;
            UInt dim = 0;

            std::vector<PointData> updateList;

            bool upper = false;
            while (current != nullptr)
            {
                updateList.push_back(current);
                parent = current;

                PointWrapper insertedPoint(point);
                PointWrapper currentPoint(current);

                if (insertedPoint.Get(dim) < currentPoint.Get(dim) || (insertedPoint.Get(dim) == currentPoint.Get(dim) && balanceByIndex[indexedPointData[current]] > 0))
                {
                    upper = false;
                    current = GetLower(current);
                }
                else
                {
                    upper = true;
                    current = GetUpper(current);
                }

                UpdateBalance(parent, upper);

                dim = (dim + 1) % Dimension;
            }

            InsertInternal(point);
            UpdateBoundaries(updateList, point);

            if (parent != nullptr)
            {
                if (upper)
                {
                    SetUpper(parent, point);
                }
                else
                {
                    SetLower(parent, point);
                }
            }
            else
            {
                origin = point;

                gbiAssert(indexedPointData.find(origin) != indexedPointData.end() && "Forgot to insert origin before assigning floor");
                nodeFloor[indexedPointData[origin]] = 0;
            }
        }
    }

    void Erase(PointData point)
    {
        PointData newLeaf = nullptr;
        bool incrementNewLeafBalance = false;

        if (indexedPointData.find(point) != indexedPointData.end())
        {
            SizeT pointIndex = indexedPointData[point];

            std::vector<SizeT> swapChain;
            SizeT current = pointIndex;

            while (!IsLeaf(current))
            {
                swapChain.push_back(current);
                current = GetBestReplacementCandidate(current);
            }
            swapChain.push_back(current);

            gbiAssert(swapChain.size() > 0);

            auto it1 = swapChain.begin();
            auto it2 = it1;
            ++it2;
            for (; it2 != swapChain.end(); ++it1, ++it2)
            {
                SwapNodes(*it1, *it2);
            }

            newLeaf = RemoveLeaf(*swapChain.rbegin(), incrementNewLeafBalance);
        }

        if (newLeaf != nullptr)
        {
            gbiAssert(indexedPointData.find(newLeaf) != indexedPointData.end());
            gbiAssert(balanceByIndex.find(indexedPointData[newLeaf]) != balanceByIndex.end());

            PointData current = newLeaf;

            UpdateBalance(newLeaf, incrementNewLeafBalance);

            while (current != nullptr)
            {
                PointData parent = GetParent(current);

                UpdateBoundaries(indexedPointData[current]);

                if (parent != nullptr)
                {
                    gbiAssert(indexedPointData.find(parent) != indexedPointData.end());

                    if (lowerIndex.find(indexedPointData[parent]) != lowerIndex.end() && pointDataVector.at(lowerIndex[indexedPointData[parent]]) == current)
                    {
                        UpdateBalance(parent, true);
                    }
                    else if (upperIndex.find(indexedPointData[parent]) != upperIndex.end() && pointDataVector.at(upperIndex[indexedPointData[parent]]) == current)
                    {
                        UpdateBalance(parent, false);
                    }
                    else
                    {
                        gbiAssert(false && "Corrupted KDTree: Parent child inconsistency");
                    }
                }

                current = parent;
            }
        }
        else
        {
            gbiAssert(pointDataVector.size() == 0);
            gbiAssert(indexedPointData.size() == 0);
            gbiAssert(nodeFloor.size() == 0);
            gbiAssert(lowerIndex.size() == 0);
            gbiAssert(upperIndex.size() == 0);
            gbiAssert(parentIndex.size() == 0);
            gbiAssert(boundaries.size() == 0);
            gbiAssert(balanceByIndex.size() == 0);
            gbiAssert(indexByBalanceAndFloor.size() == 0);

            origin = nullptr;
        }
    }

    bool RebalanceIteration()
    {
        bool didRebalance = false;

        if (indexByBalanceAndFloor.size() > 0)
        {
            auto itLowestBalance = indexByBalanceAndFloor.begin();
            auto itHighestBalance = indexByBalanceAndFloor.rbegin();

            bool rebalanceHighFound = true;
            SizeT toRebalanceHigh = 0;
            bool rebalanceLowFound = true;
            SizeT toRebalanceLow = 0;
            if (itLowestBalance != indexByBalanceAndFloor.end() && itLowestBalance->first < -1 && itLowestBalance->second.size() > 0)
            {
                auto itUpperNode = itLowestBalance->second.begin();
                gbiAssert(itUpperNode != itLowestBalance->second.end());
                gbiAssert(itUpperNode->second.begin() != itUpperNode->second.end());
                toRebalanceLow = *itUpperNode->second.begin();
            }
            else
            {
                rebalanceLowFound = false;
            }

            if (itHighestBalance != indexByBalanceAndFloor.rend() && itHighestBalance->first > 1 && itHighestBalance->second.size() > 0)
            {
                auto itUpperNode = itHighestBalance->second.begin();
                gbiAssert(itUpperNode != itHighestBalance->second.end());
                gbiAssert(itUpperNode->second.begin() != itUpperNode->second.end());
                toRebalanceHigh = *itUpperNode->second.begin();
            }
            else
            {
                rebalanceHighFound = false;
            }

            if (rebalanceLowFound || rebalanceHighFound)
            {
                PointData toRebalancePtr = nullptr;

                if (rebalanceLowFound && rebalanceHighFound)
                {
                    if (std::abs(balanceByIndex[toRebalanceLow]) > std::abs(balanceByIndex[toRebalanceHigh]))
                    {
                        toRebalancePtr = pointDataVector.at(toRebalanceLow);
                    }
                    else
                    {
                        toRebalancePtr = pointDataVector.at(toRebalanceHigh);
                    }
                }
                else if (rebalanceLowFound)
                {
                    toRebalancePtr = pointDataVector.at(toRebalanceLow);
                }
                else if (rebalanceHighFound)
                {
                    toRebalancePtr = pointDataVector.at(toRebalanceHigh);
                }

                if (toRebalancePtr != nullptr)
                {
                    Erase(toRebalancePtr);
                    Insert(toRebalancePtr);

                    didRebalance = true;
                }
            }
        }

        return didRebalance;
    }
};

}