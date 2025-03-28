#ifndef SWEEPING_H
#define SWEEPING_H
#include"mymesh.h"
#include "TransMesh.h"
#include <igl/swept_volume.h>
#include <igl/writeOBJ.h>
//#include <wrap/io_trimesh/export_obj.h>

typedef vcg::Point3f Point3f;
class SweepingMesh
{
public:
	static void GenerateSweptVolume(const Point3f& moveDir, float moveDist, CMeshO& mesh, std::vector<MyMesh::VertexHandle> boundaryPts)
	{
#if 0
		CMeshO resultMesh;
		std::map<int, int> point_map;
		for (int i = 0; i < mesh.vn; i++)
		{
			CMeshO::VertexIterator vIter = vcg::tri::Allocator<CMeshO>::AddVertex(resultMesh, mesh.vert[i].P());
			point_map[vcg::tri::Index(mesh, mesh.vert[i])] = vcg::tri::Index(resultMesh, *vIter);
		}
		for (int i = 0; i < resultMesh.vn; i++)
		{
			std::cout << " index: " << vcg::tri::Index(resultMesh, resultMesh.vert[i]) << std::endl;
		}
		for (int i = 0; i < mesh.fn; i++)
		{
			int newV[3];
			for (int j = 0; j < 3; j++)
			{
				newV[j] = point_map[vcg::tri::Index(mesh, mesh.face[i].V(j))];
			}
			vcg::tri::Allocator<CMeshO>::AddFace(resultMesh, newV[0], newV[1], newV[2]);
		}
#endif
		int count = boundaryPts.size();
		/*float avrage_length = 0.0;
		for (int i = 0; i < count; i++)
		{
			avrage_length += (boundaryPts[i].P() - boundaryPts[(i + 1) % count].P()).Norm();
		}
		avrage_length /= count;*/


		float avrage_length = 0.25;
		int num = mesh.vn;
		int add_count = (int)(moveDist / avrage_length);
		avrage_length = moveDist / add_count;
		for (int j = 0; j < add_count; j++)
		{
			for (int i = 0; i < count; i++)
			{
				CVertexO* cv = mesh.vert.data() + boundaryPts[i].idx();
				Point3f addP = cv->P() + moveDir * avrage_length * (j + 1);
				CMeshO::VertexIterator vIter = vcg::tri::Allocator<CMeshO>::AddVertex(mesh, addP);
				//std::cout << "Vertex " << i << " index: " << vcg::tri::Index(mesh, *vIter) << std::endl;
			}
		}
		for (int j = 0; j < add_count - 1; j++)
		{
			for (int i = 0; i < count; i++)
			{
				int v0 = num + j * count + i;
				int v1 = num + (j + 1) * count + i;
				int v2 = num + (j + 1) * count + (i + 1) % count;
				//std::cout << "face vertex if: " << v0 << "," << v1 << "," << v2 << std::endl;
				vcg::tri::Allocator<CMeshO>::AddFace(mesh, num + j * count + i, num + (j + 1) * count + i, num + (j + 1) * count + (i + 1) % count);
				vcg::tri::Allocator<CMeshO>::AddFace(mesh, num + j * count + i, num + (j + 1) * count + (i + 1) % count, num + j * count + (i + 1) % count);
			}
		}
		for (int i = 0; i < count; i++)
		{
			int cur_idx = boundaryPts[i].idx();
			int next_idx = boundaryPts[(i + 1) % count].idx();
			vcg::tri::Allocator<CMeshO>::AddFace(mesh, cur_idx, num + i, num + (i + 1) % count);
			vcg::tri::Allocator<CMeshO>::AddFace(mesh, cur_idx, num + (i + 1) % count, next_idx);
		}
		//std::string output_name = "D:/data/ConvexHull/SweptVolume.obj";
		//vcg::tri::io::ExporterOBJ<CMeshO>::Save(mesh, output_name.c_str(), 0);
	}

	static void MovingSweptVolume(const CMeshO& mesh, CMeshO& swept_volume,const float& offset_distance)
	{
		Eigen::MatrixXi F, SF;
		Eigen::MatrixXd V, SV, VT;
		TransMesh::GetLibiglMeshFromVcgData(mesh, V, F);
		//igl::writeOBJ("D:\\data\\ConvexHull\\output0.obj", V, F);
		const auto& transform = [](const double t)->Eigen::Affine3d
			{
				Eigen::Affine3d T = Eigen::Affine3d::Identity();
				//T.rotate(Eigen::AngleAxisd(t * 2. * igl::PI, Eigen::Vector3d(0, 1, 0)));
				T.translate(Eigen::Vector3d(-0.125 * (3. * igl::PI * t), -0.025 * (3. * igl::PI * t), -0.008 * (3. * igl::PI * t)));
				return T;
			};

		const int grid_size = 60;
		const int time_steps = 100;
		const double isolevel = offset_distance;
		igl::swept_volume(V, F, transform, time_steps, grid_size, isolevel, SV, SF);
		TransMesh::GetVcgMeshFromLibiglData(SV, SF, swept_volume);
		//igl::writeOBJ("D:\\data\\ConvexHull\\output.obj", SV, SF);
		//std::cerr << " finished." << std::endl;
	}
};
#endif
