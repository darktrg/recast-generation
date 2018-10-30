// Recast integration for UE3
// Usage:
//	1. rcSetupGeneration
//	2. rcAddExcludedAreas (optional)
//	3. rcAddWalkableSeeds (optional)
//	4. rcGenerateNavMesh
//
// All vertex coords should be passed in UE3 format

#ifndef RECAST_NAVMESH_H
#define RECAST_NAHMESH_H

#define MAX_VERTS_PER_POLY	6

// Loads navmesh generation config
// Params:
//		cellSize - (in) size of voxel 2D grid
//		cellHeight - (in) height of voxel cell
//		agentMinHeight - (in) min walkable height
//		agentMaxHeight - (in) max walkable height
//		agentMaxSlope - (in) max angle of walkable slope (degrees)
//		agentMaxClimb - (in) max height of walkable step
void __cdecl rcSetupGeneration(
	const float cellSize, const float cellHeight,
	const float agentMinHeight, const float agentMaxHeight,
	const float agentMaxSlope, const float agentMaxClimb,
	const int maxEdgeLenVx, const float maxSimplificationError, const bool simplificationWithZ, const float mergeRegionNormalThreshold);	//HS(TN) 14/04/13

struct rcConvexArea
{
	// 3D coords of 2D convex shape (z won't be used)
	float* verts;
	// number of vertices
	int nverts;
	// area bottom's height
	float minh;
	// area top's height
	float maxh;

	rcConvexArea() : verts(NULL) {}
	rcConvexArea(const rcConvexArea& srcArea) { DeepCopy(srcArea); }
	~rcConvexArea() { delete[] verts; }
	rcConvexArea& operator=(const rcConvexArea& p);
private:
	void DeepCopy(const rcConvexArea& srcArea);
};

// Builds convex area from given vertices
// Area will be always an extruded 2D shape
// Params:
//		verts - (in) vertices
//		nverts - (in) number of vertices
//		area - (out) convex area
// Returns false when ran out of memory
bool __cdecl rcBuildConvexArea(const float* verts, const int nverts, rcConvexArea& area);

// Defines areas excluded from navmesh generation
// Params:
//		areas - (in) convex areas
//		nareas - (in) number of areas
// Returns false when ran out of memory
bool __cdecl rcAddExcludedAreas(const rcConvexArea* areas, const int nareas);

// Defines accessible location on navmesh, used for filtering polys (flood fill from seeds)
// Params:
//		verts - (in) 3D coords of accessible locations
//		nverts - (in) number of vertices
//		polySearchExtent - (in) extent for searching closest polygon on navmesh
// Returns false when ran out of memory
bool __cdecl rcAddWalkableSeeds(const float* verts, const int nverts, const float* polySearchExtent);

// Filter function, allows to reject walkable height spans (2D test)
// Params:
//		x - (in) x coord of span's center
//		y - (in) y coord of span's center
//		minh - (in) z coord of voxelized area's bottom
//		maxh - (in) z coord of voxelized area's top
//		context - (in) pointer to additional data (see rcGenerateNavMesh)
// Returns true when span should be rejected
typedef bool (__cdecl *rcFilterWalkableSpan)(const float x, const float y, const float minh, const float maxh, void* context);

struct rcNavMeshPoly
{
	int verts[MAX_VERTS_PER_POLY*2];	//HS(TN) 08/02/2013 : Use the same encoding of edge like Recast
	float height;
	bool isBorderPoly;
};

// Generates navmesh
// Params:
//		verts - (in) vertices for input geometry
//		nverts - (in) number of vertices
//		faces - (in) CCW faces for input geometry
//		nfaces - (in) number of faces
//		spanFilter - (in) pointer to span filter function (optional)
//		spanFilterContext - (in) pointer to additional data passed to filter func (optional)
//		vertsNavMesh - (out) vertices in navmesh
//		nvertsNavMesh - (out) number of vertices in navmesh
//		polysNavMesh - (out) polys in navmesh
//		npolysNavMesh - (out) number of polys in navmesh
//		vertsDropDown - (out) vertices for possible drop downs
//		nvertsDropDown - (out) number of drop down vertices
// Returns true if navmesh was created
bool __cdecl rcGenerateNavMesh(
	const float* verts, const int nverts,
	const int* faces, const int nfaces, const int nunwalkablefaces,	//HS(TN) 12/04/13 : ntu = number of forced unwalkable faces (listed at the end of the array)
	const rcFilterWalkableSpan spanFilter, void* spanFilterContext,
	float*& vertsNavMesh, int& nvertsNavMesh,
	rcNavMeshPoly*& polysNavMesh, int& npolysNavMesh,
	float*& vertsDropDown, int& nvertsDropDown, int& errorCode,		//HS(TN) 08/10/13 : error code added to debug easily !
	float*& VertsDEBUG, int*& VertsIndicesDEBUG);					//HS(TN) 14/10/13 : to debug intermediate generation state

#endif
