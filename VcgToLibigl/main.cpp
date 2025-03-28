#include "TransMesh.h"
#include "mymesh.h"
#include "ConvexHull.h"
#include "MeshReConstruction.h"
#include "mymesh.h"
#include <wrap/io_trimesh/import_stl.h>
#include <wrap/io_trimesh/export_stl.h>
#include <wrap/io_trimesh/export_obj.h>
#include <wrap/io_trimesh/export_ply.h>
#include <igl/writeSTL.h>
#include <igl/copyleft/cgal/mesh_boolean.h>
#include<wrap/io_trimesh/io_ply.h>
#include<vcg/complex/algorithms/refine.h>
#include <vcg/complex/algorithms/isotropic_remeshing.h>
#include <Eigen/core>
int main()
{
#if 0
	CMeshO inputMesh;
	CMeshO inputBottomMesh;
	CMeshO outputMesh;
	int mark = 0;
	std::string input_name = "D:\\data\\ConvexHull\\11.stl";
	std::string bottom_name = "D:\\data\\ConvexHull\\11jig.stl";
	vcg::tri::io::ImporterSTL<CMeshO>::OpenBinary(inputMesh, input_name.c_str(), mark);
	vcg::tri::io::ImporterSTL<CMeshO>::OpenBinary(inputBottomMesh, bottom_name.c_str(), mark);
#endif
	//Eigen::MatrixXd VA;
	//Eigen::MatrixXi FA;
	//Eigen::MatrixXd VC;
	//Eigen::MatrixXi FC;
	//TransMesh::GetLibiglMeshFromVcgData(inputMesh, VA, FA);
	//TransMesh::GetLibiglMeshFromVcgData(inputBottomMesh, VC, FC);
	//int move_iter = 25;
	//int vn = VA.rows();
	//int fn = FA.rows();
	////检测碰撞


	//Eigen::MatrixXd VB = Eigen::MatrixXd(vn * move_iter, 3);
	//Eigen::MatrixXi FB = Eigen::MatrixXi(fn * move_iter, 3);
	//for (int i = 0; i < move_iter; i++)
	//{
	//	for (int j = 0; j < vn; j++)
	//	{
	//		int k = i * vn + j;
	//		VB(k, 0) = VA(j, 0);
	//		VB(k, 1) = VA(j, 1) + 0.03 * i;
	//		VB(k, 2) = VA(j, 2);
	//	}
	//}
	//for (int i = 0; i < move_iter; i++)
	//{
	//	for (int j = 0; j < fn; j++)
	//	{
	//		int k = i * fn + j;
	//		FB(k, 0) = FA(j, 0) + i * vn;
	//		FB(k, 1) = FA(j, 1) + i * vn;
	//		FB(k, 2) = FA(j, 2) + i * vn;
	//	}
	//}
	//Eigen::MatrixXd VD ;
	//Eigen::MatrixXi FD ;
	//igl::writeSTL("D:\\data\\ConvexHull\\move_mesh.stl", VB, FB);
	//igl::copyleft::cgal::mesh_boolean(VC, FC, VB, FB, igl::MESH_BOOLEAN_TYPE_MINUS, VD, FD);
	//igl::writeSTL("D:\\data\\ConvexHull\\result.stl", VD, FD);
	
	
	/*std::vector<CMeshO> mesh;
	mesh.push_back(inputMeshCover);
	QuickHull qh(mesh);
	qh.apply();
	CMeshO outputConVexHull = qh.getMesh();
	TransMesh::GetMeshdifferenceLibigl(inputBottomNet, outputConVexHull, outputMesh);*/
	CMeshO inputMesh;
	CMeshO inputMeshKit;
	CMeshO inputMeshBottom;
	CMeshO outputMesh;
	int mark = 0;
	std::string input_name = "D:\\data\\ConvexHull\\11\\upCover.STL";
	std::string input_name_kit = "D:\\data\\ConvexHull\\11\\upKit.STL";
	std::string input_name_bottom = "D:\\data\\ConvexHull\\11\\bottom.STL";
	vcg::tri::io::ImporterSTL<CMeshO>::OpenBinary(inputMesh, input_name.c_str(), mark);
	vcg::tri::io::ImporterSTL<CMeshO>::OpenBinary(inputMeshKit, input_name_kit.c_str(), mark);
	vcg::tri::io::ImporterSTL<CMeshO>::OpenBinary(inputMeshBottom, input_name_bottom.c_str(), mark);
#if 0 
	inputMesh.face.EnableFFAdjacency();
	inputMesh.vert.EnableVFAdjacency();
	inputMesh.face.EnableVFAdjacency();
	vcg::tri::UpdateTopology<CMeshO>::FaceFace(inputMesh);
	for (int i = 0; i < 3; ++i)
	{
		vcg::tri::Refine<CMeshO, vcg::tri::MidPoint<CMeshO> >(inputMesh, vcg::tri::MidPoint<CMeshO>(&inputMesh), 0.05, false);
	}
	std::string output_name0 = "D:\\data\\ConvexHull\\result0.obj";
	vcg::tri::io::ExporterOBJ<CMeshO>::Save(inputMesh, output_name0.c_str(), mark);
	vcg::tri::UpdateNormal<CMeshO>::PerVertexNormalized(inputMesh);
	vcg::tri::UpdateTopology<CMeshO>::FaceFace(inputMesh);
	vcg::tri::UpdateTopology<CMeshO>::VertexFace(inputMesh);
#endif
#if 0
	std::vector<int> points;
	MeshReConstruction::FilterPoints(inputMesh, points);
	CMeshO::VertexIterator chVi;
	auto viMoving = vcg::tri::Allocator<CMeshO>::AddVertices(outputMesh, points.size());
	for (size_t i = 0; i < points.size(); i++, viMoving++) 
	{
		(*viMoving).P() = inputMesh.vert[points[i]].P();
		(*viMoving).N() = inputMesh.vert[points[i]].N().Normalize();
	}
	std::string output_name = "D:\\data\\ConvexHull\\result2.ply";
	vcg::tri::io::ExporterPLY<CMeshO>::Save(outputMesh, output_name.c_str(), vcg::tri::io::Mask::IOM_VERTNORMAL);
#endif
	//vcg::tri::Clean<CMeshO>::MergeCloseVertex(inputMesh, 0.04);
	std::vector<CMeshO> combine_mesh;
	//combine_mesh.push_back(inputMesh);
	combine_mesh.push_back(inputMeshKit);
	QuickHull qh(combine_mesh);
	qh.apply();
	CMeshO ConvexHull = qh.getMesh();

	std::string output_name = "D:\\data\\ConvexHull\\ConvexHull.obj";
	vcg::tri::io::ExporterOBJ<CMeshO>::Save(ConvexHull, output_name.c_str(), 0);

	MeshReConstruction::EnableAttribute(inputMesh);
	MeshReConstruction::EnableAttribute(inputMeshKit);
	MeshReConstruction::EnableAttribute(inputMeshBottom);
	
	inputMeshKit.face.EnableFFAdjacency();
	inputMeshKit.vert.EnableVFAdjacency();
	vcg::tri::UpdateTopology<CMeshO>::FaceFace(inputMeshKit);
	vcg::tri::Refine<CMeshO, vcg::tri::MidPoint<CMeshO> >(inputMeshKit, vcg::tri::MidPoint<CMeshO>(&inputMeshKit), 0.1, false);
	CMeshO interA = inputMeshKit;
	CMeshO interB = inputMesh;
	CMeshO interC = inputMeshBottom;
	MeshReConstruction::Remeshing(interA);
	MeshReConstruction::Remeshing(interB);
	MeshReConstruction::Remeshing(interC);
	MeshReConstruction::FilterPoints(inputMesh, interA);
	MeshReConstruction::FilterPoints(inputMeshKit, interB);
	MeshReConstruction::FilterPoints(inputMeshKit, interC);
	MeshReConstruction::FilterPoints(inputMeshBottom, interA);
	MeshReConstruction::FilterPoints(inputMeshBottom, interB);
	MeshReConstruction::UpdateMesh(inputMesh);
	MeshReConstruction::UpdateMesh(inputMeshKit);
	MeshReConstruction::UpdateMesh(inputMeshBottom);
	vcg::tri::Append<CMeshO, CMeshO>::Mesh(inputMesh, inputMeshKit);
	vcg::tri::Append<CMeshO, CMeshO>::Mesh(inputMesh, inputMeshBottom);
	output_name = "D:\\data\\ConvexHull\\result0.obj";
	vcg::tri::io::ExporterOBJ<CMeshO>::Save(inputMesh, output_name.c_str(), 0);

	std::vector<int> faces;
	MeshReConstruction::FilterPoints(inputMesh, inputMeshKit);
	for (auto i : faces)
	{
		CMeshO::FaceIterator faceItr = vcg::tri::Allocator<CMeshO>::AddFaces(outputMesh, 1);
		CMeshO::VertexIterator vertexItr = vcg::tri::Allocator<CMeshO>::AddVertices(outputMesh, 3);
		for (int j = 0; j < 3; j++, vertexItr++)
		{
			(*vertexItr).P().Import(inputMesh.face[i].P(j));
			(*faceItr).V(j) = &*vertexItr;
		}
		(*faceItr).N() = inputMesh.face[i].N();
	}
	//std::string output_name = "D:\\data\\ConvexHull\\result2.obj";
	//vcg::tri::io::ExporterOBJ<CMeshO>::Save(outputMesh, output_name.c_str(), mark);

#if 0 
//凸包和网格重建同时缩放
	Eigen::Matrix4f scaleMatrix = Eigen::Matrix4f::Identity();  // 创建单位矩阵
	scaleMatrix(0, 0) = 1.05;
	scaleMatrix(1, 1) = 1.05;
	scaleMatrix(2, 2) = 1.05;

	// 遍历顶点并应用矩阵变换
	for (auto& vertex : ConvexHull.vert)
	{
		if (!vertex.IsD()) { // 确保顶点未被删除
			// 生成 4D 坐标
			Eigen::Vector4f point(vertex.P()[0], vertex.P()[1], vertex.P()[2], 1.0f);
			// 应用缩放矩阵
			Eigen::Vector4f new_point = scaleMatrix * point;
			// 更新顶点位置
			vertex.P()[0] = new_point[0];
			vertex.P()[1] = new_point[1];
			vertex.P()[2] = new_point[2];
		}
	}
#endif
}