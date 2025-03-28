#include "MeshClipping.h"
#include "MeshReConstruction.h"
#include "CleanMesh.h"
#include "CutMesh.h"
#include "FillHole.h"
#include "Sweeping.h"
#include "TransMesh.h"
#include "ConvexHull.h"
#include<vcg/complex/algorithms/refine.h>
class MeshClipping::PImpl
{
public:
	PImpl()
	{
		reConstructionMesh_.Clear();
		sweptVolume_.Clear();
		upMesh_.Clear();
		moveDist_ = 5.0f;
		scaleValue_ = 1.02f;
		isSetted = false;
	}
	~PImpl()
	{}

public:
	void deleteInnerPoint();
	void surfaceConstruction();
	void generateSweptVolume(const float& offset_distance);
	bool applyBoolOperation();
public:
	CMeshO dbMesh_;
	CMeshO olMesh_;
	CMeshO bMesh_;
	CMeshO gMesh_;

	CMeshO upMesh_;

	CMeshO reConstructionMesh_;
	CMeshO sweptVolume_;
	Point3f moveDir_;
	Point3f cutPlaneCenter_;
	Point3f cutPlaneDir_;

	float scaleValue_;
	float moveDist_;

	bool isSetted;
};

MeshClipping::MeshClipping()
{
	impl_.reset(new PImpl());
}

MeshClipping::~MeshClipping()
{

}

void MeshClipping::SetMesh(const CMeshO& dbMesh, const CMeshO& olMesh, const CMeshO& bMesh, const CMeshO& gMesh)
{
	impl_->dbMesh_  = dbMesh;
	impl_->olMesh_ = olMesh;
	impl_->bMesh_ = bMesh;
	impl_->gMesh_ = gMesh;
	vcg::tri::Append<CMeshO, CMeshO>::Mesh(impl_->upMesh_, impl_->dbMesh_);
	vcg::tri::Append<CMeshO, CMeshO>::Mesh(impl_->upMesh_, impl_->olMesh_);
	vcg::tri::Append<CMeshO, CMeshO>::Mesh(impl_->upMesh_, impl_->bMesh_);
}

void MeshClipping::SetOffsetDistance(float scale_value)
{
	impl_->scaleValue_ = scale_value;
}

CMeshO MeshClipping::GetReConstructionMesh()
{
	return impl_->reConstructionMesh_;
}

void MeshClipping::SetMove(const Point3f& bDir, float moveDist)
{
	impl_->moveDir_ = bDir;
	impl_->moveDist_ = moveDist;
}

void MeshClipping::SetCutMesh(const Point3f& planeNor, const Point3f& planeCenter)
{
	impl_->cutPlaneDir_ = planeNor;
	impl_->cutPlaneCenter_ = planeCenter;
}

void MeshClipping::ApplyClipping(const float& offset_distance)
{
	if (!impl_->isSetted)
	{
		impl_->isSetted = true;
		impl_->deleteInnerPoint();
		impl_->surfaceConstruction();
	}
	impl_->generateSweptVolume(offset_distance);
	impl_->applyBoolOperation();
	//MeshReConstruction::Remeshing(impl_->gMesh_);
}

CMeshO MeshClipping::GetMesh()
{
	return impl_->gMesh_;
}

void MeshClipping::PImpl::deleteInnerPoint()
{
	std::cout << "Start deleting interior points........." << std::endl;
	olMesh_.EnableAttribute();
	dbMesh_.EnableAttribute();
	vcg::tri::UpdateTopology<CMeshO>::FaceFace(dbMesh_);
	vcg::tri::UpdateNormal<CMeshO>::PerFace(olMesh_);
	vcg::tri::UpdateNormal<CMeshO>::PerFace(dbMesh_);
	for (int i = 0; i < 3; i++)
	{
		vcg::tri::Refine<CMeshO, vcg::tri::MidPoint<CMeshO>>(dbMesh_, vcg::tri::MidPoint<CMeshO>(&dbMesh_), 0.1, false);
		vcg::tri::Refine<CMeshO, vcg::tri::MidPoint<CMeshO>>(olMesh_, vcg::tri::MidPoint<CMeshO>(&olMesh_), 0.1, false);
	}
	CMeshO interA = dbMesh_;
	CMeshO interB = olMesh_;
	interA.EnableAttribute();
	interB.EnableAttribute();
	MeshReConstruction::FilterPoints(olMesh_, interA);
	MeshReConstruction::FilterPoints(dbMesh_, interB);
	MeshReConstruction::UpdateMesh(olMesh_);
	MeshReConstruction::UpdateMesh(dbMesh_);
	vcg::tri::Append<CMeshO, CMeshO>::Mesh(olMesh_, dbMesh_);
	/*std::vector<CMeshO> meshList;
	meshList.push_back(olMesh_);
	meshList.push_back(dbMesh_);
	QuickHull ch(meshList);
	ch.apply();
	CMeshO result = ch.getMesh();*/
	vcg::tri::io::ExporterOBJ<CMeshO>::Save(olMesh_,"D:/data/ConvexHull/olMesh_.obj", 0);
}

void MeshClipping::PImpl::surfaceConstruction()
{
	std::cout << "Start surface reconstruction........." << std::endl;
	olMesh_.EnableAttribute();
	vcg::tri::Clean<CMeshO>::SplitNonManifoldVertex(olMesh_, 0);
	vcg::tri::Clean<CMeshO>::SplitNonManifoldVertex(olMesh_, 1);
	vcg::tri::Clean<CMeshO>::SplitManifoldComponents(olMesh_);
	vcg::tri::UpdateTopology<CMeshO>::FaceFace(olMesh_);
	for (int i = 0; i < 5; i++)
	{
		vcg::tri::Refine<CMeshO, vcg::tri::MidPoint<CMeshO>>(olMesh_, vcg::tri::MidPoint<CMeshO>(&olMesh_), 0.07, false);
	}
	olMesh_.face.DisableVFAdjacency();
	olMesh_.vert.DisableVFAdjacency();
	olMesh_.updateBoxAndNormals();
	CMeshO sampleMesh, reConstructionMesh;
	MeshReConstruction::PoissonDiskSampling(olMesh_, sampleMesh);
	MeshReConstruction::PoissonReConstruction(sampleMesh, reConstructionMesh);
	MeshReConstruction::Remeshing(reConstructionMesh,0.05,10,6);
	vcg::tri::io::ExporterOBJ<CMeshO>::Save(reConstructionMesh, "D:/data/ConvexHull/reConstructionMesh.obj", 0);
	MeshReConstruction::OffsetMesh(reConstructionMesh, scaleValue_);
	CleanMesh::getMaxConnectedComponent(reConstructionMesh, reConstructionMesh_);
	vcg::tri::io::ExporterOBJ<CMeshO>::Save(reConstructionMesh, "D:/data/ConvexHull/offset_reConstructionMesh.obj", 0);
}

void MeshClipping::PImpl::generateSweptVolume(const float& offset_distance)
{
	std::cout << "Start generating swept volume........." << std::endl;
	SweepingMesh::MovingSweptVolume(reConstructionMesh_, sweptVolume_, offset_distance);
}

bool MeshClipping::PImpl::applyBoolOperation()
{
	std::cout << "Start executing Boolean operation........." << std::endl;
	Eigen::MatrixXd VA,VB,VC,VD,Vreusult;
	Eigen::MatrixXi FA,FB,FC,FD,Freusult;
	TransMesh::GetLibiglMeshFromVcgData(gMesh_, VA, FA);
	TransMesh::GetLibiglMeshFromVcgData(sweptVolume_, VB, FB);
	TransMesh::GetLibiglMeshFromVcgData(upMesh_, VC, FC);
	vcg::tri::io::ExporterOBJ<CMeshO>::Save(sweptVolume_, "D:/data/ConvexHull/sweptVolume_.obj",0);
	vcg::tri::io::ExporterOBJ<CMeshO>::Save(upMesh_, "D:/data/ConvexHull/upMesh_.obj", 0);
	vcg::tri::io::ExporterOBJ<CMeshO>::Save(gMesh_, "D:/data/ConvexHull/gMesh_.obj", 0);
	if (TransMesh::GetMeshdifferenceLibigl(VA, FA, VB, FB, VD, FD))
	{
		if (TransMesh::GetMeshdifferenceLibigl(VD, FD, VC, FC, Vreusult, Freusult))
		{
			TransMesh::GetVcgMeshFromLibiglData(Vreusult, Freusult, gMesh_);
			CMeshO result;
			CleanMesh::getMaxConnectedComponent(gMesh_, result);
			gMesh_ = result;
			std::cout << "end" << std::endl;
			return true;
		}
	}
	std::cout << "Boolean operatio error!" << std::endl;
	return false;
}