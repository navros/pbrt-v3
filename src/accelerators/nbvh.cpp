
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


// accelerators/bvh.cpp*
#include "accelerators/nbvh.h"
#include "interaction.h"
#include "paramset.h"
#include "stats.h"
#include "parallel.h"
#include <algorithm>

namespace pbrt {

STAT_MEMORY_COUNTER("Memory/BVH tree", treeBytes);
STAT_RATIO("BVH/Primitives per leaf node", totalPrimitives, totalLeafNodes);
STAT_COUNTER("BVH/Interior nodes", interiorNodes);
STAT_COUNTER("BVH/Leaf nodes", leafNodes);

STAT_COUNTER("BVH_cout/Traversation steps", traversationSteps);
STAT_COUNTER("BVH_cout/Intersection tests", intersectTests);
STAT_COUNTER("BVH_cout/Traversation steps shadow", traversationStepsShadow);
STAT_COUNTER("BVH_cout/Intersection tests shadow", intersectTestsShadow);
STAT_COUNTER("BVH_cout/Traversation steps regular", traversationStepsRegular);
STAT_COUNTER("BVH_cout/Intersection tests regular", intersectTestsRegular);

STAT_COUNTER("BVH/Build time", buildTime);
STAT_SAH_COST("SAH/SAH cost", surfaceRoot, innerSAHArea, leafSAHArea);
STAT_COUNTER("SAH/primitives - Aggregate", primitives_aggregate);
STAT_COUNTER("SAH/nodes - Interior", sah_interiorNodes);
STAT_COUNTER("SAH/nodes - Leaf", sah_leafNodes);
STAT_COUNTER("SAH/nodes - Leaf-hybrid", sah_leafHybridNodes);
STAT_COUNTER("SAH/nodes - total", bvh_nodesTotal);

STAT_COUNTER("BVH_logging/Otptimization time", optimizeTime);
STAT_COUNTER("BVH_logging/traversations count", logging_traversations);
STAT_COUNTER("BVH_logging/traversations count shadow", logging_traversationsShadow);
STAT_COUNTER("BVH_logging/traversations count regular", logging_traversationsRegular);
STAT_COUNTER("BVH_logging/primitive intersect tests", logging_primitivesTests);
STAT_COUNTER("BVH_logging/contracted nodes", logging_contractedNodes);
STAT_COUNTER("BVH_logging/swapped nodes", logging_swappedNodes);
STAT_COUNTER("BVH_logging/swapped primitives", logging_swappedPrimitives);

#define VISIT_THRESHOLD 20 //temporary.. TODO
#define NODE_N 4 //n-ary BVH nodes

/*
const unsigned char layoutLTU[64] = {
	// y=layout x=axes_direction
	// 000 001 010 011 100 101 110 111
	// order from lower significant bits to higher
	0x6C, 0xB1, 0x72, 0xC9, 0x63, 0x8D, 0x4E, 0x39,
	0x9C, 0xC9, 0x6C, 0xC6, 0x93, 0x39, 0x63, 0x36,
	0xD8, 0x8D, 0xD2, 0x2D, 0x78, 0x87, 0x72, 0x27,
	0xB4, 0x2D, 0xD8, 0x36, 0x9C, 0x27, 0x78, 0x1E,
	0x78, 0xE1, 0x4E, 0x39, 0x6C, 0xB1, 0x4B, 0x2D,
	0xE4, 0x39, 0x78, 0x1E, 0xB4, 0x2D, 0x6C, 0x1B,
	0x18, 0x21, 0x12, 0x09, 0x18, 0x21, 0x12, 0x09, // 3 child node
	0x24, 0x09, 0x18, 0x06, 0x24, 0x09, 0x18, 0x06  // 3 child node
};
*/
const unsigned char layoutLUT[256] = {
	// left to right order
	0,3,2,1,  1,0,3,2,  2,0,3,1,  1,2,0,3,  3,0,2,1,  1,3,0,2,  2,3,0,1,  1,2,3,0,
	0,3,1,2,  1,2,0,3,  0,3,2,1,  2,1,0,3,  3,0,1,2,  1,2,3,0,  3,0,2,1,  2,1,3,0,
	0,2,1,3,  1,3,0,2,  2,0,1,3,  1,3,2,0,  0,2,3,1,  3,1,0,2,  2,0,3,1,  3,1,2,0,
	0,1,3,2,  1,3,2,0,  0,2,1,3,  2,1,3,0,  0,3,1,2,  3,1,2,0,  0,2,3,1,  2,3,1,0,
	0,2,3,1,  1,0,2,3,  2,3,0,1,  1,2,3,0,  0,3,2,1,  1,0,3,2,  3,2,0,1,  1,3,2,0,
	0,1,2,3,  1,2,3,0,  0,2,3,1,  2,3,1,0,  0,1,3,2,  1,3,2,0,  0,3,2,1,  3,2,1,0,
	0,2,1,0,  1,0,2,0,  2,0,1,0,  1,2,0,0,  0,2,1,0,  1,0,2,0,  2,0,1,0,  1,2,0,0, // 3 child node
	0,1,2,0,  1,2,0,0,  0,2,1,0,  2,1,0,0,  0,1,2,0,  1,2,0,0,  0,2,1,0,  2,1,0,0  // 3 child node
};


// NBVHAccel Local Declarations
struct BVHPrimitiveInfo {
    BVHPrimitiveInfo() {}
    BVHPrimitiveInfo(size_t primitiveNumber, const Bounds3f &bounds)
        : primitiveNumber(primitiveNumber),
          bounds(bounds),
          centroid(.5f * bounds.pMin + .5f * bounds.pMax) {}
    size_t primitiveNumber;
    Bounds3f bounds;
    Point3f centroid;
};

struct BVHBuildNode {
    // BVHBuildNode Public Methods
    void InitLeaf(int first, int n, const Bounds3f &b) {
        firstPrimOffset = first;
        nPrimitives = n;
        bounds = b;
        children[0] = children[1] = nullptr;
        ++leafNodes;
        ++totalLeafNodes;
        totalPrimitives += n;
    }
    void InitInterior(int axis, BVHBuildNode *c0, BVHBuildNode *c1) {
        children[0] = c0;
        children[1] = c1;
        bounds = Union(c0->bounds, c1->bounds);
        splitAxis = axis;
        nPrimitives = 0;
        ++interiorNodes;
    }
    Bounds3f bounds;
    BVHBuildNode *children[2];
    int splitAxis, firstPrimOffset, nPrimitives;
};

struct MortonPrimitive {
    int primitiveIndex;
    uint32_t mortonCode;
};

struct LBVHTreelet {
    int startIndex, nPrimitives;
    BVHBuildNode *buildNodes;
};

struct LinearBVHNode {
    Bounds3f bounds;
    union {
        int primitivesOffset;	// leaf
        int childrenOffset;		// interior
    };
	union {
		uint16_t nPrimitives;	// leaf
		struct {				// interior
			uint8_t nodeLayout;	// 3bit layout
			uint8_t nChildren;	// amount of children
		};
	};
    uint8_t axis;				// interior node: 3*2b{xyz}, based on contraction
    uint8_t childOrder;			// leaf flag 0 or children ordering offsets (4*2b)

	unsigned char operator[](int index) const { // Get offset by child index
		return (childOrder >> (2 * index)) & 0x3; }
};

struct PairHeapCompare {
	bool  operator()(std::tuple<int, float, unsigned int>&a, std::tuple<int, float, unsigned int>&b) const {
		return std::get<1>(a) < std::get<1>(b);
	}
};

// NBVHAccel Utility Functions
inline uint32_t LeftShift3(uint32_t x) {
    CHECK_LE(x, (1 << 10));
    if (x == (1 << 10)) --x;
#ifdef PBRT_HAVE_BINARY_CONSTANTS
    x = (x | (x << 16)) & 0b00000011000000000000000011111111;
    // x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x | (x << 8)) & 0b00000011000000001111000000001111;
    // x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x | (x << 4)) & 0b00000011000011000011000011000011;
    // x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x = (x | (x << 2)) & 0b00001001001001001001001001001001;
    // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
#else
    x = (x | (x << 16)) & 0x30000ff;
    // x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x | (x << 8)) & 0x300f00f;
    // x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x | (x << 4)) & 0x30c30c3;
    // x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x = (x | (x << 2)) & 0x9249249;
    // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
#endif // PBRT_HAVE_BINARY_CONSTANTS
    return x;
}

inline uint32_t EncodeMorton3(const Vector3f &v) {
    CHECK_GE(v.x, 0);
    CHECK_GE(v.y, 0);
    CHECK_GE(v.z, 0);
    return (LeftShift3(v.z) << 2) | (LeftShift3(v.y) << 1) | LeftShift3(v.x);
}

static void RadixSort(std::vector<MortonPrimitive> *v) {
    std::vector<MortonPrimitive> tempVector(v->size());
    PBRT_CONSTEXPR int bitsPerPass = 6;
    PBRT_CONSTEXPR int nBits = 30;
    static_assert((nBits % bitsPerPass) == 0,
                  "Radix sort bitsPerPass must evenly divide nBits");
    PBRT_CONSTEXPR int nPasses = nBits / bitsPerPass;

    for (int pass = 0; pass < nPasses; ++pass) {
        // Perform one pass of radix sort, sorting _bitsPerPass_ bits
        int lowBit = pass * bitsPerPass;

        // Set in and out vector pointers for radix sort pass
        std::vector<MortonPrimitive> &in = (pass & 1) ? tempVector : *v;
        std::vector<MortonPrimitive> &out = (pass & 1) ? *v : tempVector;

        // Count number of zero bits in array for current radix sort bit
        PBRT_CONSTEXPR int nBuckets = 1 << bitsPerPass;
        int bucketCount[nBuckets] = {0};
        PBRT_CONSTEXPR int bitMask = (1 << bitsPerPass) - 1;
        for (const MortonPrimitive &mp : in) {
            int bucket = (mp.mortonCode >> lowBit) & bitMask;
            CHECK_GE(bucket, 0);
            CHECK_LT(bucket, nBuckets);
            ++bucketCount[bucket];
        }

        // Compute starting index in output array for each bucket
        int outIndex[nBuckets];
        outIndex[0] = 0;
        for (int i = 1; i < nBuckets; ++i)
            outIndex[i] = outIndex[i - 1] + bucketCount[i - 1];

        // Store sorted values in output array
        for (const MortonPrimitive &mp : in) {
            int bucket = (mp.mortonCode >> lowBit) & bitMask;
            out[outIndex[bucket]++] = mp;
        }
    }
    // Copy final result from _tempVector_, if needed
    if (nPasses & 1) std::swap(*v, tempVector);
}

// NBVHAccel Method Definitions
NBVHAccel::NBVHAccel(const std::vector<std::shared_ptr<Primitive>> &p,
                   int maxPrimsInNode, SplitMethod splitMethod, 
	               bool logging)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)),
      splitMethod(splitMethod),
      primitives(p){
    ProfilePhase _(Prof::AccelConstruction);
    if (primitives.empty()) return;
	std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
    // Build BVH from _primitives_

    // Initialize _primitiveInfo_ array for primitives
    std::vector<BVHPrimitiveInfo> primitiveInfo(primitives.size());
    for (size_t i = 0; i < primitives.size(); ++i)
        primitiveInfo[i] = {i, primitives[i]->WorldBound()};

    // Build BVH tree for primitives using _primitiveInfo_
    MemoryArena arena(1024 * 1024);
    totalNodes = 0;
    std::vector<std::shared_ptr<Primitive>> orderedPrims;
    orderedPrims.reserve(primitives.size());
    BVHBuildNode *root;
    if (splitMethod == SplitMethod::HLBVH)
        root = HLBVHBuild(arena, primitiveInfo, &totalNodes, orderedPrims);
    else
        root = recursiveBuild(arena, primitiveInfo, 0, primitives.size(),
                              &totalNodes, orderedPrims);
    primitives.swap(orderedPrims);
    LOG(INFO) << StringPrintf("BVH created with %d nodes for %d "
                              "primitives (%.2f MB)", totalNodes,
                              (int)primitives.size(),
                              float(totalNodes * sizeof(LinearBVHNode)) /
                              (1024.f * 1024.f));

    // Compute representation of depth-first traversal of BVH tree
    treeBytes += totalNodes * sizeof(LinearBVHNode) + sizeof(*this) +
                 primitives.size() * sizeof(primitives[0]);
    nodes = AllocAligned<LinearBVHNode>(totalNodes);
    int offset = 0;
	flattenNBVHTree(root, 0, &offset);
    CHECK_EQ(totalNodes, offset+1);

	// logging arrays
	totalNodesOptimized = totalNodes;
	if (logging || true) { // ..TODO hitArrays independent Optimize if no logging
		hitsNodeS.resize(totalNodes, 0);
		hitsNodeR.resize(totalNodes, 0);
		hitsPrimitiveS.resize(primitives.size(), 0);
		usedIntersectionP = &NBVHAccel::IntersectPLogging;
		usedIntersection = &NBVHAccel::IntersectLogging;
	}
	else {
		usedIntersectionP = &NBVHAccel::IntersectPRegular;
		usedIntersection = &NBVHAccel::IntersectRegular;
	}

	optimized = -1; //no phase done yet
	logging_contractedNodes = 0;
	std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
	buildTime += std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
}

void NBVHAccel::initSAHCost() const {
	// SAH cost
	bvh_nodesTotal = 1;
	std::pair<double, double> sah_sa = std::pair<double, double>(0, 0);
	ComputeSAHCost(sah_sa);
	surfaceRoot = WorldBound().SurfaceArea();
	innerSAHArea = sah_sa.first;
	leafSAHArea = sah_sa.second;
}

int NBVHAccel::ComputeSAHCost(std::pair<double, double> & sah_sa) const {
	primitives_aggregate++;
	bvh_nodesTotal += totalNodesOptimized - 1;

	for (size_t i = 0; i < totalNodesOptimized; i++) {
		if (nodes[i].childOrder > 0) {
			// interior node
			sah_interiorNodes++;
			sah_sa.first += nodes[i].bounds.SurfaceArea();
		}
		else {
			// leaf
			int truePrimitives = 0;
			Bounds3f leafBounds;
			for (size_t j = nodes[i].primitivesOffset; j < nodes[i].primitivesOffset + nodes[i].nPrimitives; j++){
				if (primitives[j]->ComputeSAHCost(sah_sa) > 0) {
					if (truePrimitives == 0)
						leafBounds = primitives[j]->WorldBound();
					else
						leafBounds = Union(leafBounds, primitives[j]->WorldBound());
					truePrimitives++;
				}
			}
			if (truePrimitives > 0)
				sah_sa.second += truePrimitives * leafBounds.SurfaceArea();
			if (truePrimitives != nodes[i].nPrimitives)
				sah_leafHybridNodes++;
			else
				sah_leafNodes++;
		}
	}
	return 0;
}

Bounds3f NBVHAccel::WorldBound() const {
    return nodes ? nodes[0].bounds : Bounds3f();
}

void NBVHAccel::Optimize(int phase) {
	// compactacion and child ordering
	if (optimized == phase)
		return;
	optimized = phase;

	std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
	logging_traversations = traversationSteps;
	logging_traversationsShadow = traversationStepsShadow;
	logging_traversationsRegular = traversationStepsRegular;
	logging_primitivesTests = intersectTests;
	
	
	int offset = 0;
	float cost = 0.0f;

	if (static_cast<OptimizePhase>(phase) == OptimizePhase::All) {
		// 1PHASE optimization
		nodesBuild = nodes;
		nodes = AllocAligned<LinearBVHNode>(totalNodes);

		int64_t primitives_hits = Contract(0, 0, &offset, cost);
		totalNodesOptimized = offset + 1;
		//TODO.. decrease size by contracted nodes
		std::cout << "root cost " << cost << " primitives hits total = " << primitives_hits << std::endl;//TODO remove it

		hitsPrimitiveS.clear();
		hitsNodeS.clear();
		hitsNodeR.clear();
		usedIntersection = &NBVHAccel::IntersectRegular;
		usedIntersectionP = &NBVHAccel::IntersectPRegular;
	}
	else {
		// 2PHASE
		int64_t primitives_hits = 0;
		if (static_cast<OptimizePhase>(optimized) == OptimizePhase::Contract) {
			nodesBuild = nodes;
			nodes = AllocAligned<LinearBVHNode>(totalNodes);

			ContractOnly(0, 0, &offset);
			totalNodesOptimized = offset + 1;
			std::fill(hitsPrimitiveS.begin(), hitsPrimitiveS.end(), 0);
			std::fill(hitsNodeS.begin(), hitsNodeS.end(), 0);
			hitsNodeR.clear();
			usedIntersection = &NBVHAccel::IntersectRegular;
			usedIntersectionP = &NBVHAccel::IntersectPLoggingMulti;
		}
		else if (static_cast<OptimizePhase>(optimized) == OptimizePhase::Reorder) {
			//swap child nodes and primitives
			primitives_hits = SortChildren(0, cost);
			hitsPrimitiveS.clear();
			hitsNodeS.clear();
			usedIntersectionP = &NBVHAccel::IntersectPRegular;
		}
		std::cout << "root cost " << cost << " primitives hits total = " << primitives_hits << std::endl;//TODO remove it
	}

	std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
	optimizeTime += std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();

	// optimize instanced BVHs in primitives
	/*
	for (int i = 0; i < primitives.size(); i++){
		primitives[i]->Optimize();
	}
	*/
}

void NBVHAccel::AddChildrenContract(
	std::vector<std::tuple<int, float, unsigned int> > &childrenOffsets, 
	int buildNodeOffset, 
	unsigned int replaceOffset) {
	// insert children offsets with calculated hit probability (both shadow and regular rays statistic)
	LinearBVHNode* node = &nodesBuild[buildNodeOffset];
	unsigned int nodeOffset;
	for (int i = 0; i < node->nChildren; i++) {
		nodeOffset = i == 0 ? replaceOffset : childrenOffsets.size();
		if (nodesBuild[node->childrenOffset + i].childOrder == 0) {
			// lowest priority in heap for leaves => no contraction
			childrenOffsets.push_back(std::tuple<int, float, unsigned int>(
				node->childrenOffset + i, 0.0f, nodeOffset));
		}
		else {
			if (hitsNodeS[node->childrenOffset + i] + hitsNodeR[node->childrenOffset + i] != 0)//TODO ?: > VISIT_THRESHOLD
				childrenOffsets.push_back(std::tuple<int, float, unsigned int>(
					node->childrenOffset + i,
					(float)(hitsNodeS[node->childrenOffset + i] + hitsNodeR[node->childrenOffset + i]) 
					/ (float)(hitsNodeS[buildNodeOffset] + hitsNodeR[buildNodeOffset]), 
					nodeOffset));
			else {
				// too few visits hence probability set by surface area
				// or for render with no logging phase
				childrenOffsets.push_back(std::tuple<int, float, unsigned int>(
					node->childrenOffset + i,
					nodesBuild[node->childrenOffset + i].bounds.SurfaceArea() / nodesBuild[buildNodeOffset].bounds.SurfaceArea(), 
					nodeOffset));
			}
		}
		std::push_heap(childrenOffsets.begin(), childrenOffsets.end(), PairHeapCompare()); // heapify
	}
}

int64_t NBVHAccel::Contract(int buildNodeOffset, int nodeOffset, int *offset, float &cost) {
	cost = 1.0f; // Ct
	int64_t primitives_hits = 0;
	LinearBVHNode *buildNode = &nodesBuild[buildNodeOffset];
	LinearBVHNode *linearNode = &nodes[nodeOffset];
	linearNode->bounds = buildNode->bounds;

	if (buildNode->childOrder == 0) {
		// leaf node

		// change the primitives order by hit probability
		int end_offset = buildNode->primitivesOffset + buildNode->nPrimitives;
		for (int i = buildNode->primitivesOffset; i < end_offset; i++) {
			for (int j = i; j > buildNode->primitivesOffset; --j) {
				if (hitsPrimitiveS[j-1] < hitsPrimitiveS[j]) {
					++logging_swappedPrimitives;
					std::swap(hitsPrimitiveS[j-1], hitsPrimitiveS[j]);
					std::swap(primitives[j-1], primitives[j]);
				}
			}
		}

		// new leaf information
		linearNode->childOrder = 0;
		linearNode->nPrimitives = buildNode->nPrimitives;
		linearNode->primitivesOffset = buildNode->primitivesOffset;
		for (int i = 0; i < buildNode->nPrimitives; i++)
			primitives_hits += hitsPrimitiveS[buildNode->primitivesOffset + i];
		cost = buildNode->nPrimitives;
	}
	else {
		// interior node - contract and order children
		linearNode->nodeLayout = 0;
		linearNode->axis = buildNode->axis & 0x3;
		std::vector<std::tuple<int, float, unsigned int> > childrenOffsets; // tuple of node offset, cost value, child index
		AddChildrenContract(childrenOffsets, buildNodeOffset, 0);

		// interior nodes contraction
		while (std::get<1>(childrenOffsets.front()) > 1.0f - 1.0f / childrenOffsets.size()
			&& childrenOffsets.size() != NODE_N) { // while not all leaves or full
			++logging_contractedNodes;

			// select next node for contraction
			int contractSelected = std::get<0>(childrenOffsets.front());
			unsigned int contractedChild = std::get<2>(childrenOffsets.front());

			// remove selected node from set
			std::pop_heap(childrenOffsets.begin(), childrenOffsets.end(), PairHeapCompare());
			childrenOffsets.pop_back();

			// modify current node and add children of contracted node
			linearNode->nodeLayout += childrenOffsets.size() == 1 ? contractedChild : (contractedChild << 1); // expecting max 2 contraction
			linearNode->axis += (nodesBuild[contractSelected].axis & 0x3) << 2 * childrenOffsets.size();
			AddChildrenContract(childrenOffsets, contractSelected, contractedChild);
		}
		
		// Tree recursion and child nodes ordering by hit probability for shadow rays
		int myOffset = *offset + 1;
		(*offset) += childrenOffsets.size();
		std::vector<float> H, I, C;
		H.reserve(childrenOffsets.size());
		I.reserve(childrenOffsets.size());
		C.reserve(childrenOffsets.size());
		float parentHits = (float)hitsNodeS[buildNodeOffset];
		for (int i = 0; i < childrenOffsets.size(); i++) {
			// recurse
			int64_t primitiv_hits_i = Contract(std::get<0>(childrenOffsets[i]), myOffset + i, offset, C[i]);
			if (hitsNodeS[buildNodeOffset] > 0) {
				H[i] = (float)primitiv_hits_i / parentHits;
				I[i] = (float)hitsNodeS[std::get<0>(childrenOffsets[i])] / parentHits;
				std::get<1>(childrenOffsets[i]) = I[i] == 0.0f ? 0.0f : H[i] / (I[i]*C[i]);
			}
			else {
				H[i] = 0.0f;
				I[i] = 0.0f;
				std::get<1>(childrenOffsets[i]) = 0.0f;
			}

			// update current node values
			primitives_hits += primitiv_hits_i;
			std::get<0>(childrenOffsets[i]) = myOffset + i; // change to new array offset for sorting
		}
		
		// change the childern order by overall score
		for (int i = 0; i < childrenOffsets.size(); i++) {
			for (int j = i; j > 0; --j) {
				if (std::get<1>(childrenOffsets[j - 1]) < std::get<1>(childrenOffsets[j])) {
					++logging_swappedNodes;
					std::swap(std::get<1>(childrenOffsets[j - 1]), std::get<1>(childrenOffsets[j]));
					std::swap(std::get<2>(childrenOffsets[j - 1]), std::get<2>(childrenOffsets[j]));
					std::swap(nodes[std::get<0>(childrenOffsets[j-1])], nodes[std::get<0>(childrenOffsets[j])]);
					std::swap(H[j - 1], H[j]);
					std::swap(I[j - 1], I[j]);
					std::swap(C[j - 1], C[j]);
				}
			}
		}

		// child ordering after swapping
		linearNode->childOrder = 0;
		for (int i = 0; i < childrenOffsets.size(); i++)
			linearNode->childOrder += i << 2 * std::get<2>(childrenOffsets[i]);

		if (childrenOffsets.size() == 2) {
			// exception for binary node - fake 000 layout with just 2 possible orderings
			linearNode->childOrder += linearNode->childOrder << 4; // 4 upper bits repeat
			linearNode->axis = buildNode->axis; // 3 times repeated
		}
		if (childrenOffsets.size() == 3)
			// redefine layout LTU index for 3-children node
			linearNode->nodeLayout = linearNode->nodeLayout == 0 ? 0x6 : 0x7;

		// final cost
		float product_H = 1.0f;
		for (int i = 0; i < childrenOffsets.size() - 1; i++) {
			cost += I[i] * C[i] * product_H;
			product_H *= (1.0f - H[i]);
		}
		
		// nodes stats
		linearNode->childrenOffset = myOffset;
		linearNode->nChildren = childrenOffsets.size();
	}
	return primitives_hits;
}


void NBVHAccel::ContractOnly(int buildNodeOffset, int nodeOffset, int *offset) {
	LinearBVHNode *buildNode = &nodesBuild[buildNodeOffset];
	LinearBVHNode *linearNode = &nodes[nodeOffset];
	linearNode->bounds = buildNode->bounds;

	if (buildNode->childOrder == 0) {
		// leaf node
		linearNode->childOrder = 0;
		linearNode->nPrimitives = buildNode->nPrimitives;
		linearNode->primitivesOffset = buildNode->primitivesOffset;
	}
	else {
		// interior node - contract and order children
		linearNode->nodeLayout = 0;
		linearNode->axis = buildNode->axis & 0x3;
		std::vector<std::tuple<int, float, unsigned int> > childrenOffsets; // tuple of node offset, cost value, child index
		AddChildrenContract(childrenOffsets, buildNodeOffset, 0);

		// interior nodes contraction
		while (std::get<1>(childrenOffsets.front()) > 1.0f - 1.0f / childrenOffsets.size()
			&& childrenOffsets.size() != NODE_N) { // while not all leaves or full
			++logging_contractedNodes;

			// select next node for contraction
			int contractSelected = std::get<0>(childrenOffsets.front());
			unsigned int contractedChild = std::get<2>(childrenOffsets.front());

			// remove selected node from set
			std::pop_heap(childrenOffsets.begin(), childrenOffsets.end(), PairHeapCompare());
			childrenOffsets.pop_back();

			// modify current node and add children of contracted node
			linearNode->nodeLayout += childrenOffsets.size() == 1 ? contractedChild : (contractedChild << 1); // expecting max 2 contraction
			linearNode->axis += (nodesBuild[contractSelected].axis & 0x3) << 2 * childrenOffsets.size();
			AddChildrenContract(childrenOffsets, contractSelected, contractedChild);
		}

		// Tree recursion and child nodes ordering by hit probability for shadow rays
		int myOffset = *offset + 1;
		(*offset) += childrenOffsets.size();
		for (int i = 0; i < childrenOffsets.size(); i++) {
			// recurse
			ContractOnly(std::get<0>(childrenOffsets[i]), myOffset + i, offset);
		}

		// child ordering value
		linearNode->childOrder = 0;
		for (int i = 0; i < childrenOffsets.size(); i++)
			linearNode->childOrder += i << 2 * std::get<2>(childrenOffsets[i]);

		if (childrenOffsets.size() == 2) {
			// exception for binary node - fake 000 layout with just 2 possible orderings
			linearNode->childOrder += linearNode->childOrder << 4; // 4 upper bits repeat
			linearNode->axis = buildNode->axis; // 3 times repeated
		}
		if (childrenOffsets.size() == 3)
			// redefine layout LTU index for 3-children node
			linearNode->nodeLayout = linearNode->nodeLayout == 0 ? 0x6 : 0x7;

		// nodes stats
		linearNode->childrenOffset = myOffset;
		linearNode->nChildren = childrenOffsets.size();
	}
}

int64_t NBVHAccel::SortChildren(int nodeOffset, float &cost) {
	cost = 1.0f; // Ct
	int64_t primitives_hits = 0;
	LinearBVHNode *linearNode = &nodes[nodeOffset];
	
	if (linearNode->childOrder == 0) {
		// leaf node

		// change the primitives order by hit probability
		int end_offset = linearNode->primitivesOffset + linearNode->nPrimitives;
		for (int i = linearNode->primitivesOffset; i < end_offset; i++) {
			for (int j = i; j > linearNode->primitivesOffset; --j) {
				if (hitsPrimitiveS[j - 1] < hitsPrimitiveS[j]) {
					++logging_swappedPrimitives;
					std::swap(hitsPrimitiveS[j - 1], hitsPrimitiveS[j]);
					std::swap(primitives[j - 1], primitives[j]);
				}
			}
		}
		for (int i = 0; i < linearNode->nPrimitives; i++)
			primitives_hits += hitsPrimitiveS[linearNode->primitivesOffset + i];
		cost = linearNode->nPrimitives;
	}
	else {
		// interior node
		
		std::vector<std::tuple<int, float, unsigned int> > childrenOffsets; // tuple of node offset, cost value, child index
		childrenOffsets.resize(linearNode->nChildren);
		for (unsigned int i = 0; i < linearNode->nChildren; i++) {
			int nodeOffset = (linearNode->childOrder >> 2 * i) & 0x3;
			childrenOffsets[nodeOffset] = std::tuple<int, float, unsigned int>(
				linearNode->childrenOffset + nodeOffset, 0.0f, i);
		}
		/* THE SAME
		for (unsigned int offset = 0; offset < linearNode->nChildren; offset++) {
			unsigned int childIndex = 0;
			for (; childIndex < linearNode->nChildren; childIndex++)
				if (offset == ((linearNode->childOrder >> 2 * childIndex) & 0x3))
					break;
			childrenOffsets.push_back(std::tuple<int, float, unsigned int>(
				linearNode->childrenOffset + offset, 0.0f, childIndex));
		}
		*/

		// Tree recursion and child nodes ordering by hit probability for shadow rays
		std::vector<float> H, I, C;
		H.reserve(childrenOffsets.size());
		I.reserve(childrenOffsets.size());
		C.reserve(childrenOffsets.size());
		float parentHits = (float)hitsNodeS[nodeOffset];
		for (int i = 0; i < childrenOffsets.size(); i++) {
			// recurse
			int64_t primitiv_hits_i = SortChildren(std::get<0>(childrenOffsets[i]), C[i]);
			if (parentHits > 0) {
				H[i] = (float)primitiv_hits_i / parentHits;
				I[i] = (float)hitsNodeS[std::get<0>(childrenOffsets[i])] / parentHits;
				std::get<1>(childrenOffsets[i]) = I[i] == 0.0f ? 0.0f : H[i] / (I[i] * C[i]);
			}
			else {
				H[i] = 0.0f;
				I[i] = 0.0f;
				std::get<1>(childrenOffsets[i]) = 0.0f;
			}

			// update current node values
			primitives_hits += primitiv_hits_i;
		}

		// change the childern order by overall score
		for (int i = 0; i < childrenOffsets.size(); i++) {
			for (int j = i; j > 0; --j) {
				if (std::get<1>(childrenOffsets[j - 1]) < std::get<1>(childrenOffsets[j])) {
					++logging_swappedNodes;
					std::swap(std::get<1>(childrenOffsets[j - 1]), std::get<1>(childrenOffsets[j]));
					std::swap(std::get<2>(childrenOffsets[j - 1]), std::get<2>(childrenOffsets[j]));
					std::swap(nodes[std::get<0>(childrenOffsets[j - 1])], nodes[std::get<0>(childrenOffsets[j])]);
					std::swap(H[j - 1], H[j]);
					std::swap(I[j - 1], I[j]);
					std::swap(C[j - 1], C[j]);
				}
			}
		}

		// child ordering after swapping
		linearNode->childOrder = 0;
		for (int i = 0; i < childrenOffsets.size(); i++)
			linearNode->childOrder += i << 2 * std::get<2>(childrenOffsets[i]);

		if (childrenOffsets.size() == 2) {
			// exception for binary node - fake 000 layout with just 2 possible orderings
			linearNode->childOrder += linearNode->childOrder << 4; // 4 upper bits repeat
		}

		// final cost
		float product_H = 1.0f;
		for (int i = 0; i < childrenOffsets.size() - 1; i++) {
			cost += I[i] * C[i] * product_H;
			product_H *= (1.0f - H[i]);
		}
	}
	return primitives_hits;
}


struct BucketInfo {
    int count = 0;
    Bounds3f bounds;
};

BVHBuildNode *NBVHAccel::recursiveBuild(
    MemoryArena &arena, std::vector<BVHPrimitiveInfo> &primitiveInfo, int start,
    int end, int *totalNodes,
    std::vector<std::shared_ptr<Primitive>> &orderedPrims) {
    CHECK_NE(start, end);
    BVHBuildNode *node = arena.Alloc<BVHBuildNode>();
    (*totalNodes)++;
    // Compute bounds of all primitives in BVH node
    Bounds3f bounds;
    for (int i = start; i < end; ++i)
        bounds = Union(bounds, primitiveInfo[i].bounds);
    int nPrimitives = end - start;
    if (nPrimitives == 1) {
        // Create leaf _BVHBuildNode_
        int firstPrimOffset = orderedPrims.size();
        for (int i = start; i < end; ++i) {
            int primNum = primitiveInfo[i].primitiveNumber;
            orderedPrims.push_back(primitives[primNum]);
        }
        node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
        return node;
    } else {
        // Compute bound of primitive centroids, choose split dimension _dim_
        Bounds3f centroidBounds;
        for (int i = start; i < end; ++i)
            centroidBounds = Union(centroidBounds, primitiveInfo[i].centroid);
        int dim = centroidBounds.MaximumExtent();

        // Partition primitives into two sets and build children
        int mid = (start + end) / 2;
        if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
            // Create leaf _BVHBuildNode_
            int firstPrimOffset = orderedPrims.size();
            for (int i = start; i < end; ++i) {
                int primNum = primitiveInfo[i].primitiveNumber;
                orderedPrims.push_back(primitives[primNum]);
            }
            node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
            return node;
        } else {
            // Partition primitives based on _splitMethod_
            switch (splitMethod) {
            case SplitMethod::Middle: {
                // Partition primitives through node's midpoint
                Float pmid =
                    (centroidBounds.pMin[dim] + centroidBounds.pMax[dim]) / 2;
                BVHPrimitiveInfo *midPtr = std::partition(
                    &primitiveInfo[start], &primitiveInfo[end - 1] + 1,
                    [dim, pmid](const BVHPrimitiveInfo &pi) {
                        return pi.centroid[dim] < pmid;
                    });
                mid = midPtr - &primitiveInfo[0];
                // For lots of prims with large overlapping bounding boxes, this
                // may fail to partition; in that case don't break and fall
                // through
                // to EqualCounts.
                if (mid != start && mid != end) break;
            }
            case SplitMethod::EqualCounts: {
                // Partition primitives into equally-sized subsets
                mid = (start + end) / 2;
                std::nth_element(&primitiveInfo[start], &primitiveInfo[mid],
                                 &primitiveInfo[end - 1] + 1,
                                 [dim](const BVHPrimitiveInfo &a,
                                       const BVHPrimitiveInfo &b) {
                                     return a.centroid[dim] < b.centroid[dim];
                                 });
                break;
            }
            case SplitMethod::SAH:
            default: {
                // Partition primitives using approximate SAH
                if (nPrimitives <= 2) {
                    // Partition primitives into equally-sized subsets
                    mid = (start + end) / 2;
                    std::nth_element(&primitiveInfo[start], &primitiveInfo[mid],
                                     &primitiveInfo[end - 1] + 1,
                                     [dim](const BVHPrimitiveInfo &a,
                                           const BVHPrimitiveInfo &b) {
                                         return a.centroid[dim] <
                                                b.centroid[dim];
                                     });
                } else {
                    // Allocate _BucketInfo_ for SAH partition buckets
                    PBRT_CONSTEXPR int nBuckets = 12;
                    BucketInfo buckets[nBuckets];

                    // Initialize _BucketInfo_ for SAH partition buckets
                    for (int i = start; i < end; ++i) {
                        int b = nBuckets *
                                centroidBounds.Offset(
                                    primitiveInfo[i].centroid)[dim];
                        if (b == nBuckets) b = nBuckets - 1;
                        CHECK_GE(b, 0);
                        CHECK_LT(b, nBuckets);
                        buckets[b].count++;
                        buckets[b].bounds =
                            Union(buckets[b].bounds, primitiveInfo[i].bounds);
                    }

                    // Compute costs for splitting after each bucket
                    Float cost[nBuckets - 1];
                    for (int i = 0; i < nBuckets - 1; ++i) {
                        Bounds3f b0, b1;
                        int count0 = 0, count1 = 0;
                        for (int j = 0; j <= i; ++j) {
                            b0 = Union(b0, buckets[j].bounds);
                            count0 += buckets[j].count;
                        }
                        for (int j = i + 1; j < nBuckets; ++j) {
                            b1 = Union(b1, buckets[j].bounds);
                            count1 += buckets[j].count;
                        }
                        cost[i] = 1 +
                                  (count0 * b0.SurfaceArea() +
                                   count1 * b1.SurfaceArea()) /
                                      bounds.SurfaceArea();
                    }

                    // Find bucket to split at that minimizes SAH metric
                    Float minCost = cost[0];
                    int minCostSplitBucket = 0;
                    for (int i = 1; i < nBuckets - 1; ++i) {
                        if (cost[i] < minCost) {
                            minCost = cost[i];
                            minCostSplitBucket = i;
                        }
                    }

                    // Either create leaf or split primitives at selected SAH
                    // bucket
                    Float leafCost = nPrimitives;
                    if (nPrimitives > maxPrimsInNode || minCost < leafCost) {
                        BVHPrimitiveInfo *pmid = std::partition(
                            &primitiveInfo[start], &primitiveInfo[end - 1] + 1,
                            [=](const BVHPrimitiveInfo &pi) {
                                int b = nBuckets *
                                        centroidBounds.Offset(pi.centroid)[dim];
                                if (b == nBuckets) b = nBuckets - 1;
                                CHECK_GE(b, 0);
                                CHECK_LT(b, nBuckets);
                                return b <= minCostSplitBucket;
                            });
                        mid = pmid - &primitiveInfo[0];
                    } else {
                        // Create leaf _BVHBuildNode_
                        int firstPrimOffset = orderedPrims.size();
                        for (int i = start; i < end; ++i) {
                            int primNum = primitiveInfo[i].primitiveNumber;
                            orderedPrims.push_back(primitives[primNum]);
                        }
                        node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
                        return node;
                    }
                }
                break;
            }
            }
            node->InitInterior(dim,
                               recursiveBuild(arena, primitiveInfo, start, mid,
                                              totalNodes, orderedPrims),
                               recursiveBuild(arena, primitiveInfo, mid, end,
                                              totalNodes, orderedPrims));
        }
    }
    return node;
}

BVHBuildNode *NBVHAccel::HLBVHBuild(
    MemoryArena &arena, const std::vector<BVHPrimitiveInfo> &primitiveInfo,
    int *totalNodes,
    std::vector<std::shared_ptr<Primitive>> &orderedPrims) const {
    // Compute bounding box of all primitive centroids
    Bounds3f bounds;
    for (const BVHPrimitiveInfo &pi : primitiveInfo)
        bounds = Union(bounds, pi.centroid);

    // Compute Morton indices of primitives
    std::vector<MortonPrimitive> mortonPrims(primitiveInfo.size());
    ParallelFor([&](int i) {
        // Initialize _mortonPrims[i]_ for _i_th primitive
        PBRT_CONSTEXPR int mortonBits = 10;
        PBRT_CONSTEXPR int mortonScale = 1 << mortonBits;
        mortonPrims[i].primitiveIndex = primitiveInfo[i].primitiveNumber;
        Vector3f centroidOffset = bounds.Offset(primitiveInfo[i].centroid);
        mortonPrims[i].mortonCode = EncodeMorton3(centroidOffset * mortonScale);
    }, primitiveInfo.size(), 512);

    // Radix sort primitive Morton indices
    RadixSort(&mortonPrims);

    // Create LBVH treelets at bottom of BVH

    // Find intervals of primitives for each treelet
    std::vector<LBVHTreelet> treeletsToBuild;
    for (int start = 0, end = 1; end <= (int)mortonPrims.size(); ++end) {
#ifdef PBRT_HAVE_BINARY_CONSTANTS
      uint32_t mask = 0b00111111111111000000000000000000;
#else
      uint32_t mask = 0x3ffc0000;
#endif
      if (end == (int)mortonPrims.size() ||
            ((mortonPrims[start].mortonCode & mask) !=
             (mortonPrims[end].mortonCode & mask))) {
            // Add entry to _treeletsToBuild_ for this treelet
            int nPrimitives = end - start;
            int maxBVHNodes = 2 * nPrimitives;
            BVHBuildNode *nodes = arena.Alloc<BVHBuildNode>(maxBVHNodes, false);
            treeletsToBuild.push_back({start, nPrimitives, nodes});
            start = end;
        }
    }

    // Create LBVHs for treelets in parallel
    std::atomic<int> atomicTotal(0), orderedPrimsOffset(0);
    orderedPrims.resize(primitives.size());
    ParallelFor([&](int i) {
        // Generate _i_th LBVH treelet
        int nodesCreated = 0;
        const int firstBitIndex = 29 - 12;
        LBVHTreelet &tr = treeletsToBuild[i];
        tr.buildNodes =
            emitLBVH(tr.buildNodes, primitiveInfo, &mortonPrims[tr.startIndex],
                     tr.nPrimitives, &nodesCreated, orderedPrims,
                     &orderedPrimsOffset, firstBitIndex);
        atomicTotal += nodesCreated;
    }, treeletsToBuild.size());
    *totalNodes = atomicTotal;

    // Create and return SAH BVH from LBVH treelets
    std::vector<BVHBuildNode *> finishedTreelets;
    finishedTreelets.reserve(treeletsToBuild.size());
    for (LBVHTreelet &treelet : treeletsToBuild)
        finishedTreelets.push_back(treelet.buildNodes);
    return buildUpperSAH(arena, finishedTreelets, 0, finishedTreelets.size(),
                         totalNodes);
}

BVHBuildNode *NBVHAccel::emitLBVH(
    BVHBuildNode *&buildNodes,
    const std::vector<BVHPrimitiveInfo> &primitiveInfo,
    MortonPrimitive *mortonPrims, int nPrimitives, int *totalNodes,
    std::vector<std::shared_ptr<Primitive>> &orderedPrims,
    std::atomic<int> *orderedPrimsOffset, int bitIndex) const {
    CHECK_GT(nPrimitives, 0);
    if (bitIndex == -1 || nPrimitives < maxPrimsInNode) {
        // Create and return leaf node of LBVH treelet
        (*totalNodes)++;
        BVHBuildNode *node = buildNodes++;
        Bounds3f bounds;
        int firstPrimOffset = orderedPrimsOffset->fetch_add(nPrimitives);
        for (int i = 0; i < nPrimitives; ++i) {
            int primitiveIndex = mortonPrims[i].primitiveIndex;
            orderedPrims[firstPrimOffset + i] = primitives[primitiveIndex];
            bounds = Union(bounds, primitiveInfo[primitiveIndex].bounds);
        }
        node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
        return node;
    } else {
        int mask = 1 << bitIndex;
        // Advance to next subtree level if there's no LBVH split for this bit
        if ((mortonPrims[0].mortonCode & mask) ==
            (mortonPrims[nPrimitives - 1].mortonCode & mask))
            return emitLBVH(buildNodes, primitiveInfo, mortonPrims, nPrimitives,
                            totalNodes, orderedPrims, orderedPrimsOffset,
                            bitIndex - 1);

        // Find LBVH split point for this dimension
        int searchStart = 0, searchEnd = nPrimitives - 1;
        while (searchStart + 1 != searchEnd) {
            CHECK_NE(searchStart, searchEnd);
            int mid = (searchStart + searchEnd) / 2;
            if ((mortonPrims[searchStart].mortonCode & mask) ==
                (mortonPrims[mid].mortonCode & mask))
                searchStart = mid;
            else {
                CHECK_EQ(mortonPrims[mid].mortonCode & mask,
                         mortonPrims[searchEnd].mortonCode & mask);
                searchEnd = mid;
            }
        }
        int splitOffset = searchEnd;
        CHECK_LE(splitOffset, nPrimitives - 1);
        CHECK_NE(mortonPrims[splitOffset - 1].mortonCode & mask,
                 mortonPrims[splitOffset].mortonCode & mask);

        // Create and return interior LBVH node
        (*totalNodes)++;
        BVHBuildNode *node = buildNodes++;
        BVHBuildNode *lbvh[2] = {
            emitLBVH(buildNodes, primitiveInfo, mortonPrims, splitOffset,
                     totalNodes, orderedPrims, orderedPrimsOffset,
                     bitIndex - 1),
            emitLBVH(buildNodes, primitiveInfo, &mortonPrims[splitOffset],
                     nPrimitives - splitOffset, totalNodes, orderedPrims,
                     orderedPrimsOffset, bitIndex - 1)};
        int axis = bitIndex % 3;
        node->InitInterior(axis, lbvh[0], lbvh[1]);
        return node;
    }
}

BVHBuildNode *NBVHAccel::buildUpperSAH(MemoryArena &arena,
                                      std::vector<BVHBuildNode *> &treeletRoots,
                                      int start, int end,
                                      int *totalNodes) const {
    CHECK_LT(start, end);
    int nNodes = end - start;
    if (nNodes == 1) return treeletRoots[start];
    (*totalNodes)++;
    BVHBuildNode *node = arena.Alloc<BVHBuildNode>();

    // Compute bounds of all nodes under this HLBVH node
    Bounds3f bounds;
    for (int i = start; i < end; ++i)
        bounds = Union(bounds, treeletRoots[i]->bounds);

    // Compute bound of HLBVH node centroids, choose split dimension _dim_
    Bounds3f centroidBounds;
    for (int i = start; i < end; ++i) {
        Point3f centroid =
            (treeletRoots[i]->bounds.pMin + treeletRoots[i]->bounds.pMax) *
            0.5f;
        centroidBounds = Union(centroidBounds, centroid);
    }
    int dim = centroidBounds.MaximumExtent();
    // FIXME: if this hits, what do we need to do?
    // Make sure the SAH split below does something... ?
    CHECK_NE(centroidBounds.pMax[dim], centroidBounds.pMin[dim]);

    // Allocate _BucketInfo_ for SAH partition buckets
    PBRT_CONSTEXPR int nBuckets = 12;
    struct BucketInfo {
        int count = 0;
        Bounds3f bounds;
    };
    BucketInfo buckets[nBuckets];

    // Initialize _BucketInfo_ for HLBVH SAH partition buckets
    for (int i = start; i < end; ++i) {
        Float centroid = (treeletRoots[i]->bounds.pMin[dim] +
                          treeletRoots[i]->bounds.pMax[dim]) *
                         0.5f;
        int b =
            nBuckets * ((centroid - centroidBounds.pMin[dim]) /
                        (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
        if (b == nBuckets) b = nBuckets - 1;
        CHECK_GE(b, 0);
        CHECK_LT(b, nBuckets);
        buckets[b].count++;
        buckets[b].bounds = Union(buckets[b].bounds, treeletRoots[i]->bounds);
    }

    // Compute costs for splitting after each bucket
    Float cost[nBuckets - 1];
    for (int i = 0; i < nBuckets - 1; ++i) {
        Bounds3f b0, b1;
        int count0 = 0, count1 = 0;
        for (int j = 0; j <= i; ++j) {
            b0 = Union(b0, buckets[j].bounds);
            count0 += buckets[j].count;
        }
        for (int j = i + 1; j < nBuckets; ++j) {
            b1 = Union(b1, buckets[j].bounds);
            count1 += buckets[j].count;
        }
        cost[i] = .125f +
                  (count0 * b0.SurfaceArea() + count1 * b1.SurfaceArea()) /
                      bounds.SurfaceArea();
    }

    // Find bucket to split at that minimizes SAH metric
    Float minCost = cost[0];
    int minCostSplitBucket = 0;
    for (int i = 1; i < nBuckets - 1; ++i) {
        if (cost[i] < minCost) {
            minCost = cost[i];
            minCostSplitBucket = i;
        }
    }

    // Split nodes and create interior HLBVH SAH node
    BVHBuildNode **pmid = std::partition(
        &treeletRoots[start], &treeletRoots[end - 1] + 1,
        [=](const BVHBuildNode *node) {
            Float centroid =
                (node->bounds.pMin[dim] + node->bounds.pMax[dim]) * 0.5f;
            int b = nBuckets *
                    ((centroid - centroidBounds.pMin[dim]) /
                     (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
            if (b == nBuckets) b = nBuckets - 1;
            CHECK_GE(b, 0);
            CHECK_LT(b, nBuckets);
            return b <= minCostSplitBucket;
        });
    int mid = pmid - &treeletRoots[0];
    CHECK_GT(mid, start);
    CHECK_LT(mid, end);
    node->InitInterior(
        dim, this->buildUpperSAH(arena, treeletRoots, start, mid, totalNodes),
        this->buildUpperSAH(arena, treeletRoots, mid, end, totalNodes));
    return node;
}

void NBVHAccel::flattenNBVHTree(BVHBuildNode *node, int nodeOffset, int *offset) {
	LinearBVHNode *linearNode = &nodes[nodeOffset];
	linearNode->bounds = node->bounds;
	if (node->nPrimitives > 0) {
		// Create leaf flattened BVH node
		CHECK(!node->children[0] && !node->children[1]);
		CHECK_LT(node->nPrimitives, 65536);
		linearNode->primitivesOffset = node->firstPrimOffset;
		linearNode->nPrimitives = node->nPrimitives;
		linearNode->childOrder = 0;
	}
	else {
		int myOffset = *offset + 1;
		(*offset) += 2;
		// Create interior flattened BVH node
		linearNode->nodeLayout = 0;
		linearNode->axis = node->splitAxis 
			+ (node->splitAxis << 2) 
			+ (node->splitAxis << 4); // three times same axis
		linearNode->nChildren = 2;
		linearNode->childrenOffset = myOffset;
		linearNode->childOrder = 0x44; // as node layout 000 : N1 N0 N1 N0
		flattenNBVHTree(node->children[0], myOffset, offset);
		flattenNBVHTree(node->children[1], myOffset + 1, offset);
	}
}

NBVHAccel::~NBVHAccel() { FreeAligned(nodes); }

bool NBVHAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
	return (this->*usedIntersection)(this, ray, isect);
}

bool NBVHAccel::IntersectRegular(const NBVHAccel* nbvh, const Ray &ray, SurfaceInteraction *isect) const {
	if (!nodes) return false;
	ProfilePhase p(Prof::AccelIntersect);
	bool hit = false;
	Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
	int dirIsNeg[4] = { invDir.x < 0, invDir.y < 0, invDir.z < 0, 0 };
	unsigned char dirIsNegLUT[43];
	for (size_t i = 0; i < 43; i++)
		dirIsNegLUT[i] = dirIsNeg[i % 4] * 4 + dirIsNeg[(i / 4) % 4] * 8 + dirIsNeg[i / 16] * 16;
	// Follow ray through BVH nodes to find primitive intersections
	int toVisitOffset = 0, currentNodeIndex = 0;
	int nodesToVisit[64 * NODE_N];
	while (true) {
		const LinearBVHNode *node = &nodes[currentNodeIndex];
		++traversationSteps;
		++traversationStepsRegular;//TODO will be one .. later
		// Check ray against BVH node
		if (node->bounds.IntersectP(ray, invDir, dirIsNeg)) {
			if (node->childOrder == 0) {
				// Intersect ray with primitives in leaf BVH node
				intersectTests += node->nPrimitives;
				intersectTestsRegular += node->nPrimitives;
				for (int i = 0; i < node->nPrimitives; ++i)
					if (primitives[node->primitivesOffset + i]->Intersect(ray, isect))
						hit = true;
				if (toVisitOffset == 0) break;
				currentNodeIndex = nodesToVisit[--toVisitOffset];
			}
			else {
				unsigned int axisDir = dirIsNegLUT[node->axis];
				int orderLTUindex = 32 * node->nodeLayout + axisDir;

				currentNodeIndex = node->childrenOffset + ((node->childOrder >> (2 * layoutLUT[orderLTUindex])) & 0x3);
				for (unsigned int chilnNo = node->nChildren - 1; chilnNo > 0; --chilnNo)
					nodesToVisit[toVisitOffset++] = node->childrenOffset + ((node->childOrder >> (2 * layoutLUT[orderLTUindex + chilnNo])) & 0x3);
				/*
				// WITHOUT dirIsNeg lookup table :
				unsigned int axisDir = dirIsNeg[node->axis & 0x3]
					+ (dirIsNeg[(node->axis >> 2) & 0x3] << 1)
					+ (dirIsNeg[(node->axis >> 4) & 0x3] << 2);

				unsigned int traverseOrder = layoutLTU[(node->nodeLayout << 3) + axisDir];

				currentNodeIndex = node->childrenOffset + ((node->childOrder >> ((traverseOrder & 0x3) << 1)) & 0x3);
				for (unsigned int chilnNo = node->nChildren - 1; chilnNo > 0; --chilnNo)
					nodesToVisit[toVisitOffset++] = node->childrenOffset + ((node->childOrder >> (((traverseOrder >> (chilnNo << 1)) & 0x3) << 1)) & 0x3);
				*/
				/*
				// SORTED :
				// traversation by sorted distances to child BBs

				// intersection distance to BB (0 for ray inside)
				int indexes[NODE_N];
				Float distances[NODE_N];
				for (size_t i = 0; i < node->nChildren; i++) {
					indexes[i] = i;
					nodes[node->childrenOffset + i].bounds.IntersectP(ray, distances + i);
				}

				// push to stack by the sorted distances (overall minimum is next node to visit)
				// TODO: optimize... slow (just traverse count reference)
				int maximum_i = 0;
				for (size_t i = 0; i < node->nChildren; i++) {
					maximum_i = i;
					for (size_t j = i + 1; j < node->nChildren; j++) {
						if (distances[j] > distances[maximum_i])
							maximum_i = j;
					}
					std::swap(distances[i], distances[maximum_i]);
					std::swap(indexes[i], indexes[maximum_i]);

					if (i == node->nChildren - 1)
						currentNodeIndex = node->childrenOffset + indexes[i];
					else
						nodesToVisit[toVisitOffset++] = node->childrenOffset + indexes[i];
				}
				*/
			}
		}
		else {
			if (toVisitOffset == 0) break;
			currentNodeIndex = nodesToVisit[--toVisitOffset];
		}
	}
	return hit;
}

bool NBVHAccel::IntersectLogging(const NBVHAccel* nbvh, const Ray &ray, SurfaceInteraction *isect) const {
	// logging with binary BVH
	if (!nodes) return false;
	ProfilePhase p(Prof::AccelIntersect);
	bool hit = false;
	Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
	int dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
	// Follow ray through BVH nodes to find primitive intersections
	int toVisitOffset = 0, currentNodeIndex = 0;
	int nodesToVisit[64];
	auto thisBVH = const_cast<NBVHAccel*>(nbvh); // remove const in current bvh for logging
	while (true) {
		const LinearBVHNode *node = &nodes[currentNodeIndex];
		++traversationSteps;
		++traversationStepsRegular;
		// Check ray against BVH node
		if (node->bounds.IntersectP(ray, invDir, dirIsNeg)) {
			//++(thisBVH->hitsNodeR[currentNodeIndex]); // TODO: check, ... delete if optimization only for shadow rays
			if (node->childOrder == 0) {
				// Intersect ray with primitives in leaf BVH node
				intersectTests += node->nPrimitives;
				intersectTestsRegular += node->nPrimitives;
				for (int i = 0; i < node->nPrimitives; ++i)
					if (primitives[node->primitivesOffset + i]->Intersect(ray, isect))
						hit = true;
				if (toVisitOffset == 0) break;
				currentNodeIndex = nodesToVisit[--toVisitOffset];
			}
			else {
				if (dirIsNeg[node->axis & 0x3]) {
                    nodesToVisit[toVisitOffset++] = node->childrenOffset;
					currentNodeIndex = node->childrenOffset + 1;
                } else {
                    nodesToVisit[toVisitOffset++] = node->childrenOffset + 1;
                    currentNodeIndex = node->childrenOffset;
                }
			}
		}
		else {
			if (toVisitOffset == 0) break;
			currentNodeIndex = nodesToVisit[--toVisitOffset];
		}
	}
	return hit;
}

bool NBVHAccel::IntersectP(const Ray &ray) const {
	return (this->*usedIntersectionP)(this, ray);
}

bool NBVHAccel::IntersectPRegular(const NBVHAccel* nbvh, const Ray &ray) const{
	if (!nodes) return false;
	ProfilePhase p(Prof::AccelIntersectP);
	Vector3f invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
	int dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
	std::pair<int, int> nodesToVisit[64]; // node offset start and nodes offset end
	int toVisitOffset = -1, currentNodeIndex = 0;
	while (true) {
		const LinearBVHNode *node = &nodes[currentNodeIndex];
		++traversationSteps;
		++traversationStepsShadow;
		if (node->bounds.IntersectP(ray, invDir, dirIsNeg)) {
			// Process BVH node _node_ for traversal
			if (node->childOrder == 0) {
				for (int i = 0; i < node->nPrimitives; ++i) {
					++intersectTests;
					++intersectTestsShadow;
					if (primitives[node->primitivesOffset + i]->IntersectP(
						ray)) {
						return true;
					}
				}
				if (toVisitOffset == -1) break;
				currentNodeIndex = nodesToVisit[toVisitOffset].first++;
				if (nodesToVisit[toVisitOffset].first == nodesToVisit[toVisitOffset].second) --toVisitOffset;
			} else {
				currentNodeIndex = node->childrenOffset;
				nodesToVisit[++toVisitOffset] = std::pair<int, int>(node->childrenOffset + 1, node->childrenOffset + node->nChildren);
			}
		} else {
			if (toVisitOffset == -1) break;
			currentNodeIndex = nodesToVisit[toVisitOffset].first++;
			if (nodesToVisit[toVisitOffset].first == nodesToVisit[toVisitOffset].second) --toVisitOffset;
		}
	}
	return false;
}

bool NBVHAccel::IntersectPLoggingMulti(const NBVHAccel* nbvh, const Ray &ray) const {
	if (!nodes) return false;
	ProfilePhase p(Prof::AccelIntersectP);
	Vector3f invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
	int dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
	int nodesToVisit[64 * NODE_N];
	int toVisitOffset = 0, currentNodeIndex = 0;
	auto thisBVH = const_cast<NBVHAccel*>(nbvh); // remove const in current bvh for logging
	while (true) {
		const LinearBVHNode *node = &nodes[currentNodeIndex];
		++traversationSteps;
		++traversationStepsShadow;
		if (node->bounds.IntersectP(ray, invDir, dirIsNeg)) {
			// Process BVH node _node_ for traversal
			++(thisBVH->hitsNodeS[currentNodeIndex]);
			if (node->childOrder == 0) {
				for (int i = 0; i < node->nPrimitives; ++i) {
					++intersectTests;
					++intersectTestsShadow;
					if (primitives[node->primitivesOffset + i]->IntersectP(ray)) {
						++(thisBVH->hitsPrimitiveS[node->primitivesOffset + i]);
						return true;
					}
				}
				if (toVisitOffset == 0) break;
				currentNodeIndex = nodesToVisit[--toVisitOffset];
			}
			else {
				// TODO random
				currentNodeIndex = node->childrenOffset;
				for (int i = 1; i < node->nChildren; i++)
					nodesToVisit[toVisitOffset++] = node->childrenOffset + i;
			}
		}
		else {
			if (toVisitOffset == 0) break;
			currentNodeIndex = nodesToVisit[--toVisitOffset];
		}
	}
	return false;
}

bool NBVHAccel::IntersectPLogging(const NBVHAccel* nbvh, const Ray &ray) const{
	// logging with binary BVH
    if (!nodes) return false;
    ProfilePhase p(Prof::AccelIntersectP);
    Vector3f invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
    int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};
	int nodesToVisit[64];
    int toVisitOffset = 0, currentNodeIndex = 0;
	auto thisBVH = const_cast<NBVHAccel*>(nbvh); // remove const in current bvh for logging
    while (true) {
        const LinearBVHNode *node = &nodes[currentNodeIndex];
		++traversationSteps;
		++traversationStepsShadow;
        if (node->bounds.IntersectP(ray, invDir, dirIsNeg)) {
            // Process BVH node _node_ for traversal
			++(thisBVH->hitsNodeS[currentNodeIndex]);
            if (node->childOrder == 0) {
                for (int i = 0; i < node->nPrimitives; ++i) {
					++intersectTests;
					++intersectTestsShadow;
                    if (primitives[node->primitivesOffset + i]->IntersectP(ray)) {
						++(thisBVH->hitsPrimitiveS[node->primitivesOffset + i]);
                        return true;
                    }
                }
				if (toVisitOffset == 0) break;
				currentNodeIndex = nodesToVisit[--toVisitOffset];
			} else {
				//if (rand() < 0.5f * RAND_MAX) { // random shuffle..temporary disabled for accurate results
				if (dirIsNeg[node->axis & 0x3]) {
                    /// second child first
                    nodesToVisit[toVisitOffset++] = node->childrenOffset;
					currentNodeIndex = node->childrenOffset + 1;
					
                } else {
                    nodesToVisit[toVisitOffset++] = node->childrenOffset + 1;
                    currentNodeIndex = node->childrenOffset;
                }
			}
		} else {
			if (toVisitOffset == 0) break;
			currentNodeIndex = nodesToVisit[--toVisitOffset];
		}
	}
    return false;
}


std::shared_ptr<NBVHAccel> CreateNBVHAccelerator(
    const std::vector<std::shared_ptr<Primitive>> &prims, const ParamSet &ps) {
    std::string splitMethodName = ps.FindOneString("splitmethod", "sah");
	bool logging = ps.FindOneBool("logging", true);
    NBVHAccel::SplitMethod splitMethod;
    if (splitMethodName == "sah")
        splitMethod = NBVHAccel::SplitMethod::SAH;
    else if (splitMethodName == "hlbvh")
        splitMethod = NBVHAccel::SplitMethod::HLBVH;
    else if (splitMethodName == "middle")
        splitMethod = NBVHAccel::SplitMethod::Middle;
    else if (splitMethodName == "equal")
        splitMethod = NBVHAccel::SplitMethod::EqualCounts;
    else {
        Warning("BVH split method \"%s\" unknown.  Using \"sah\".",
                splitMethodName.c_str());
        splitMethod = NBVHAccel::SplitMethod::SAH;
    }

    int maxPrimsInNode = ps.FindOneInt("maxnodeprims", 4);
    return std::make_shared<NBVHAccel>(prims, maxPrimsInNode, splitMethod, logging);
}

}  // namespace pbrt
