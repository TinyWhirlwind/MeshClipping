#ifndef MESH_CLIPPING_H
#define MESH_CLIPPING_H
#include "mymesh.h"
using namespace vcg;
class MeshClipping
{
public:
	MeshClipping();
	~MeshClipping();

public:
	/**
	 * @brief 设置网格
	 * @param dbMesh 基托
	 * @param olMesh 盖片
	 * @param bMesh 网底
	 * @param gMesh 导板
	 */
	void SetMesh(const CMeshO& dbMesh, const CMeshO& olMesh, const CMeshO& bMesh, const CMeshO& gMesh);
	void SetOffsetDistance(float scale_value);
	CMeshO GetReConstructionMesh();
	/**
	 * @brief 设置移动信息
	 * @param moveDir 移动方向
	 * @param moveDist 移动距离
	 */
	void SetMove(const Point3f& moveDir, float moveDist);
	void SetCutMesh(const Point3f& planeNor, const Point3f& planeCenter);
	/**
	 * @brief clipping mesh
	 */
	void ApplyClipping(const float& offset_distance);
	CMeshO GetMesh();
private:
	class PImpl;
	std::shared_ptr<PImpl> impl_;
};
#endif