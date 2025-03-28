#include <iostream>
#include <pcl/io/pcd_io.h>
#include <pcl/io/vtk_lib_io.h>
#include <pcl/point_cloud.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <boost/thread/thread.hpp>
#include <vtkSmartPointer.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkPolyData.h>
#include <pcl/surface/concave_hull.h>
#include <pcl/common/intersections.h>
using namespace std;
using namespace pcl;
int main()
{
    // 读取STL格式模型
    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName("D:/data/ConvexHull/union2.stl");
    reader->Update();
    vtkSmartPointer<vtkPolyData> polydata = reader->GetOutput();

    if (!polydata) {
        std::cerr << "Failed to load STL file." << std::endl;
        return -1;
    }

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>());
    pcl::io::vtkPolyDataToPointCloud(polydata, *cloud);

    // 创建PCLVisualizer显示STL文件
    boost::shared_ptr<pcl::visualization::PCLVisualizer> stl_viewer(new pcl::visualization::PCLVisualizer("STL Viewer"));
    stl_viewer->addModelFromPolyData(polydata, "mesh");

    // 创建新的线程显示STL窗口
    boost::thread stl_thread([&stl_viewer]()
        {
            while (!stl_viewer->wasStopped()) {
                stl_viewer->spinOnce(100);
                boost::this_thread::sleep(boost::posix_time::microseconds(100000));
            }
        });

    PolygonMesh mesh_out;
    boost::shared_ptr<pcl::ConcaveHull<PointXYZ>> concave_hull(new pcl::ConcaveHull<PointXYZ>());
    concave_hull->setInputCloud(cloud);
    concave_hull->setDimension(3);
    concave_hull->setAlpha(0.08);
    concave_hull->reconstruct(mesh_out);
    vtkSmartPointer<vtkPolyData> polydata0 = vtkSmartPointer<vtkPolyData>::New();
    pcl::io::mesh2vtk(mesh_out, polydata0);
    pcl::io::savePolygonFileSTL("D:/data/ConvexHull/output_concave_hull.stl", mesh_out);
    // 创建PCLVisualizer显示concave hull文件
    boost::shared_ptr<pcl::visualization::PCLVisualizer> hull_viewer(new pcl::visualization::PCLVisualizer("Concave Hull Viewer"));
    hull_viewer->addModelFromPolyData(polydata0, "mesh_out");

    std::string outputFilename = "D:/data/ConvexHull/output_concave_hull.stl";
    vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
    writer->SetFileName(outputFilename.c_str());
    writer->SetInputData(polydata0);
    writer->Write();

    // 主要线程处理PCD显示
    while (!hull_viewer->wasStopped()) {
        hull_viewer->spinOnce(100);
        boost::this_thread::sleep(boost::posix_time::microseconds(100000));
    }

    // 保存PCD文件
    pcl::io::savePCDFileASCII("D:/data/ConvexHull/11.pcd", *cloud);

    // 加载点云文件
    pcl::io::loadPCDFile("D:/data/ConvexHull/11.pcd", *cloud);
    std::cout << "Loaded PCD with " << cloud->points.size() << " points." << std::endl;

    // 创建PCLVisualizer显示PCD文件
    boost::shared_ptr<pcl::visualization::PCLVisualizer> pcd_viewer(new pcl::visualization::PCLVisualizer("PCD Viewer"));
    pcd_viewer->addPointCloud(cloud, "cloud");

    // 主要线程处理PCD显示
    while (!pcd_viewer->wasStopped()) {
        pcd_viewer->spinOnce(100);
        boost::this_thread::sleep(boost::posix_time::microseconds(100000));
    }

    // 等待STL显示线程结束
    stl_thread.join();
    return 0;
}

