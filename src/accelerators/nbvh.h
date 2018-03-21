
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_ACCELERATORS_NBVH_H
#define PBRT_ACCELERATORS_NBVH_H

// accelerators/bvh.h*
#include "pbrt.h"
#include "primitive.h"
#include <atomic>

namespace pbrt {
struct BVHBuildNode;

// NBVHAccel Forward Declarations
struct BVHPrimitiveInfo;
struct MortonPrimitive;
struct LinearBVHNode;
struct PairHeapCompare;

// NBVHAccel Declarations
class NBVHAccel : public Aggregate {
  public:
    // NBVHAccel Public Types
    enum class SplitMethod { SAH, HLBVH, Middle, EqualCounts };

    // NBVHAccel Public Methods
    NBVHAccel(const std::vector<std::shared_ptr<Primitive>> &p,
             int maxPrimsInNode = 1,
             SplitMethod splitMethod = SplitMethod::SAH,
			 bool logging = false);
    Bounds3f WorldBound() const;
    ~NBVHAccel();
    bool Intersect(const Ray &ray, SurfaceInteraction *isect) const;
    bool IntersectP(const Ray &ray) const;
	void initSAHCost() const;
	int ComputeSAHCost(std::pair<double, double> & sah_sa) const;
	void Optimize();

	// make hits as float (for direct fraction use during computation) is 
	// no good (incementing +1 have no effect for number >= 2^24)
	std::vector<int64_t> hitsNodeS;		// node shadow ray hits 
	std::vector<int64_t> hitsPrimitiveS;
	std::vector<int64_t> hitsNodeR;		// node regular ray hits

  private:
    // NBVHAccel Private Methods
    BVHBuildNode *recursiveBuild(
        MemoryArena &arena, std::vector<BVHPrimitiveInfo> &primitiveInfo,
        int start, int end, int *totalNodes,
        std::vector<std::shared_ptr<Primitive>> &orderedPrims);
    BVHBuildNode *HLBVHBuild(
        MemoryArena &arena, const std::vector<BVHPrimitiveInfo> &primitiveInfo,
        int *totalNodes,
        std::vector<std::shared_ptr<Primitive>> &orderedPrims) const;
    BVHBuildNode *emitLBVH(
        BVHBuildNode *&buildNodes,
        const std::vector<BVHPrimitiveInfo> &primitiveInfo,
        MortonPrimitive *mortonPrims, int nPrimitives, int *totalNodes,
        std::vector<std::shared_ptr<Primitive>> &orderedPrims,
        std::atomic<int> *orderedPrimsOffset, int bitIndex) const;
    BVHBuildNode *buildUpperSAH(MemoryArena &arena,
                                std::vector<BVHBuildNode *> &treeletRoots,
                                int start, int end, int *totalNodes) const;
    void flattenNBVHTree(BVHBuildNode *node, int nodeOffset, int *offset);

	typedef bool (NBVHAccel::*ShadowIntersection)(const NBVHAccel* nbvh, const Ray &ray) const;
	typedef bool (NBVHAccel::*RegularIntersection)(const NBVHAccel* nbvh, const Ray &ray, SurfaceInteraction *isect) const;
	ShadowIntersection usedIntersectionP;
	RegularIntersection usedIntersection;
	bool IntersectPRegular(const NBVHAccel* nbvh, const Ray &ray) const;
	bool IntersectPLogging(const NBVHAccel* nbvh, const Ray &ray) const;
	bool IntersectRegular(const NBVHAccel* nbvh, const Ray &ray, SurfaceInteraction *isect) const;
	bool IntersectLogging(const NBVHAccel* nbvh, const Ray &ray, SurfaceInteraction *isect) const;
	void AddChildrenContract(
		std::vector<std::tuple<int, float, unsigned int> > &childrenOffsets, 
		int buildNodeOffset, unsigned int replaceOffset);
	int64_t Contract(int buildNodeOffset, int nodeOffset, int *offset, float &cost);

    // NBVHAccel Private Data
	int totalNodes;
	int totalNodesOptimized;
    const int maxPrimsInNode;
    const SplitMethod splitMethod;
    std::vector<std::shared_ptr<Primitive>> primitives;
    LinearBVHNode *nodes = nullptr;
	LinearBVHNode *nodesBuild = nullptr;
};

std::shared_ptr<NBVHAccel> CreateNBVHAccelerator(
    const std::vector<std::shared_ptr<Primitive>> &prims, const ParamSet &ps);

}  // namespace pbrt

#endif  // PBRT_ACCELERATORS_BVH_H
