//
// Created by liguang on 2022/2/19.
//

#ifndef TEST_OPENCV_SINGLE_PICTURE_SOLVER_H
#define TEST_OPENCV_SINGLE_PICTURE_SOLVER_H

#include <opencv2/opencv.hpp>
#include <opencv2/xfeatures2d.hpp>
#include <ceres/ceres.h>
//#include <Eigen/Eigen>
#include <iostream>
#include <chrono>
#include <fstream>

#define PI 3.1415926

//ceres
struct Camera_Dial_Cost {
    Camera_Dial_Cost(double* p, double* d) : pixel_xy(p), device_xyz(d) {}

    template<typename T>
    bool operator()(
            const T* const eulur,
            const T* const translation,
            T* residual
    )const {
        T R[9] =
                {
                        ceres::cos(eulur[1]) * ceres::cos(eulur[2]),
                        ceres::sin(eulur[0]) * ceres::sin(eulur[1]) * ceres::cos(eulur[2]) - ceres::cos(eulur[0]) * ceres::sin(eulur[2]),
                        ceres::cos(eulur[0]) * ceres::sin(eulur[1]) * ceres::cos(eulur[2]) + ceres::sin(eulur[0]) * ceres::sin(eulur[2]),
                        ceres::cos(eulur[1]) * ceres::sin(eulur[2]),
                        ceres::sin(eulur[0]) * ceres::sin(eulur[1]) * ceres::sin(eulur[2]) + ceres::cos(eulur[0]) * ceres::cos(eulur[2]),
                        ceres::cos(eulur[0]) * ceres::sin(eulur[1]) * ceres::sin(eulur[2]) - ceres::sin(eulur[0]) * ceres::cos(eulur[2]),
                        -ceres::sin(eulur[1]),
                        ceres::sin(eulur[0]) * ceres::cos(eulur[1]),
                        ceres::cos(eulur[0]) * ceres::cos(eulur[1])
                };

        //T Fxyz[3] = {
        //	(R[0] * ((T)device_xyz[0] - translation[0]) + R[3] * ((T)device_xyz[1] - translation[1]) + R[6] * ((T)device_xyz[2] - translation[2])),
        //	(R[1] * ((T)device_xyz[0] - translation[0]) + R[4] * ((T)device_xyz[1] - translation[1]) + R[7] * ((T)device_xyz[2] - translation[2])),
        //	(R[2] * ((T)device_xyz[0] - translation[0]) + R[5] * ((T)device_xyz[1] - translation[1]) + R[8] * ((T)device_xyz[2] - translation[2]))
        //};

        T Fxyz[3] = {
                (R[0] * (T)device_xyz[0] + R[1] * (T)device_xyz[1] + R[2] * (T)device_xyz[2] + translation[0]),
                (R[3] * (T)device_xyz[0] + R[4] * (T)device_xyz[1] + R[5] * (T)device_xyz[2] + translation[1]),
                (R[6] * (T)device_xyz[0] + R[7] * (T)device_xyz[1] + R[8] * (T)device_xyz[2] + translation[2])
        };

        residual[0] = (T)(pixel_xy[0]) - (Fxyz[0] / Fxyz[2]);
        residual[1] = (T)(pixel_xy[1]) - (Fxyz[1] / Fxyz[2]);
        return true;
    }

    const double* pixel_xy, * device_xyz;
};
struct Plane_Fit {
    Plane_Fit(double* p) : plane_xyz(p) {}

    template<typename T>
    bool operator()(
            const T* const m,
            T* residual
    )const {

        //T Z_plane = (T)plane_xyz[0] * m[0] + (T)plane_xyz[1] * m[1] + m[2];
        residual[0] = (T)plane_xyz[0] * m[0] + (T)plane_xyz[1] * m[1] + m[2] - (T)plane_xyz[2];
        return true;
    }

    const double* plane_xyz;
};
struct dist
{
    double dx;
    double dy;
};

namespace single_solver
{
    class solver {
    public:
        //Eigen::Matrix3f euler_to_matrix(double eulur[3]);
        dist update_xy(const double* device_xyz, const double* eulur, const double* translation);
        void XwToXc(const double* device_xyz, const double* eulur, const double* translation, double(&Xc)[3]);
        void point2d_to_point3d(double pic_x, double pic_y, double Ceulur[3],double Ctranslation[3],double m[3], double (&p3)[3]);
        void undisrortion(double pixel_x, double pixel_y, double K[9], double(&pixel)[2]);
        void cameraPoseFromHomography(const cv::Mat& H, cv::Mat& pose);
        //solve camera pose
        int solve_camera_pose(std::string file_path,double(&Ceulur)[3],double(&Ctranslation)[3],double(&m)[3],double(&K)[9]);
        //solve 3D point clouds
        int solve_3D_point_clouds(std::string image_path,std::string points3D_path,std::vector<cv::Point>& board_points,
                                  double(&K)[9],double(&Ceulur)[3],double(&Ctranslation)[3],double(&m)[3]);
        //solve homography
        int solve_Homography(std::string image1_path,std::string image2_path,cv::Mat& H);
        //solve Essential Matrix to R T
        int solve_Essential_to_pose(std::string image1_path,std::string image2_path,Eigen::Matrix3d& R,Eigen::Vector3d& T,double(&K)[9]);
        //solve R T from Homography opencv
        //vector<Mat> Rs;
        //vector<Mat> Ts;
        //decomposeHomographyMat(projectionMatrix, cameraMatrix, Rs, Ts, noArray());
        int solve_pose_from_opencv(std::string image1_path,std::string image2_path,cv::Mat& K,
                                   std::vector<cv::Mat>& Rs,std::vector<cv::Mat>& Ts,Eigen::Matrix3f& RR,Eigen::Vector3f& TT);


    };
}


#endif //TEST_OPENCV_SINGLE_PICTURE_SOLVER_H
