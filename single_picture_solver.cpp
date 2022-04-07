//
// Created by liguang on 2022/2/19.
//

#include "single_picture_solver.h"
#include "p3p_kneip.h"
//#include <opencv2/core/eigen.hpp>


namespace single_solver
{
//    Eigen::Matrix3f solver::euler_to_matrix(double *eulur)
//    {
//        Eigen::Matrix3f m3;
//        m3(0, 0) = std::cos(eulur[1]) * std::cos(eulur[2]);
//        m3(0, 1) = std::sin(eulur[0]) * std::sin(eulur[1]) * std::cos(eulur[2]) - std::cos(eulur[0]) * std::sin(eulur[2]);
//        m3(0, 2) = std::cos(eulur[0]) * std::sin(eulur[1]) * std::cos(eulur[2]) + std::sin(eulur[0]) * std::sin(eulur[2]);
//        m3(1, 0) = std::cos(eulur[1]) * std::sin(eulur[2]);
//        m3(1, 1) = std::sin(eulur[0]) * std::sin(eulur[1]) * std::sin(eulur[2]) + std::cos(eulur[0]) * std::cos(eulur[2]);
//        m3(1, 2) = std::cos(eulur[0]) * std::sin(eulur[1]) * std::sin(eulur[2]) - std::sin(eulur[0]) * std::cos(eulur[2]);
//        m3(2, 0) = -std::sin(eulur[1]);
//        m3(2, 1) = std::sin(eulur[0]) * std::cos(eulur[1]);
//        m3(2, 2) = std::cos(eulur[0]) * std::cos(eulur[1]);
//
//        return m3;
//    }

    dist solver::update_xy(const double *device_xyz, const double *eulur, const double *translation)
    {
        double R[9] =
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
        double Fxyz[3] = {
                (R[0] * device_xyz[0] + R[1] * device_xyz[1] + R[2] * device_xyz[2] + translation[0]),
                (R[3] * device_xyz[0] + R[4] * device_xyz[1] + R[5] * device_xyz[2] + translation[1]),
                (R[6] * device_xyz[0] + R[7] * device_xyz[1] + R[8] * device_xyz[2] + translation[2])
        };
        double residual[2];
        residual[0] = (Fxyz[0] / Fxyz[2]);
        residual[1] = (Fxyz[1] / Fxyz[2]);
        dist dist1;
        dist1.dx = residual[0];
        dist1.dy = residual[1];
        return dist1;
    }

    void solver::XwToXc(const double *device_xyz, const double *eulur, const double *translation,double (&Xc)[3])
    {
        double R[9] =
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
        double Fxyz[3] = {
                (R[0] * device_xyz[0] + R[1] * device_xyz[1] + R[2] * device_xyz[2] + translation[0]),
                (R[3] * device_xyz[0] + R[4] * device_xyz[1] + R[5] * device_xyz[2] + translation[1]),
                (R[6] * device_xyz[0] + R[7] * device_xyz[1] + R[8] * device_xyz[2] + translation[2])
        };
        Xc[0] = Fxyz[0];
        Xc[1] = Fxyz[1];
        Xc[2] = Fxyz[2];
    }

    void solver::point2d_to_point3d(double pic_x, double pic_y, double *Ceulur, double *Ctranslation,double *m, double (&p3)[3])
    {
        Eigen::Vector3d pic_xy;
        pic_xy << pic_x, pic_y, 1;

        //?????????? euler
        Eigen::Matrix3d R_xyz ;
       // = euler_to_matrix(Ceulur);

        //??????
        Eigen::Vector3d d_xyz;
        d_xyz << Ctranslation[0], Ctranslation[1], Ctranslation[2];
        //Solve

        Eigen::Matrix3d R1 = R_xyz;
        Eigen::Vector3d T1 = d_xyz;

        double m0 = m[0];
        double m1 = m[1];
        double m2 = m[2];

        double A1 = (R1(0, 0) + R1(0, 2) * m0) - pic_x * (R1(2, 0) + R1(2, 2) * m0);
        double A2 = (R1(0, 1) + R1(0, 2) * m1) - pic_x * (R1(2, 1) + R1(2, 2) * m1);
        double A3 = R1(0, 2) * m2 + T1(0, 0) - pic_x * (R1(2, 2) * m2 + T1(2, 0));

        double B1 = (R1(1, 0) + R1(1, 2) * m0) - pic_y * (R1(2, 0) + R1(2, 2) * m0);
        double B2 = (R1(1, 1) + R1(1, 2) * m1) - pic_y * (R1(2, 1) + R1(2, 2) * m1);
        double B3 = R1(1, 2) * m2 + T1(1, 0) - pic_y * (R1(2, 2) * m2 + T1(2, 0));

        p3[0] = (B3 * A2 - B2 * A3) / (double)(A1 * B2 - A2 * B1);
        p3[1] = (B3 * A1 - B1 * A3) / (double)(A2 * B1 - A1 * B2);
        p3[2] = m0 * p3[0] + m1 * p3[1] + m2;
    }

    void solver::undisrortion(double pixel_x, double pixel_y, double *K, double (&pixel)[2])
    {
        //?????1
        double _u = (pixel_x - K[0]) / K[2];
        double _v = (pixel_y - K[1]) / K[3];
        double _iter_u = _u;
        double _iter_v = _v;
        //cout << "pixelX:  " << _u << " " << "pixelY: " << _v << endl;
        double _r = _u * _u + _v * _v;
        double _du = _u * (K[4] * _r + K[5] * _r * _r + K[8] * _r * _r * _r) + 2 * K[6] * _u * _v + K[7] * (_r + 2 * _u * _u);
        double _dv = _v * (K[4] * _r + K[5] * _r * _r + K[8] * _r * _r * _r) + 2 * K[7] * _u * _v + K[6] * (_r + 2 * _v * _v);
        double eslpse = 1e-8;
        double _iter_du = _du + eslpse;
        double _iter_dv = _dv + eslpse;
        while (abs(_iter_du - _du) + abs(_iter_dv - _dv) > 1e-16) {
            _iter_u = _u - _du;
            _iter_v = _v - _dv;
            _r = _iter_u * _iter_u + _iter_v * _iter_v;
            _iter_du = _du;
            _iter_dv = _dv;
            _du = _iter_u * (K[4] * _r + K[5] * _r * _r + K[8] * _r * _r * _r) + 2 * K[6] * _iter_u * _iter_v + K[7] * (_r + 2 * _iter_u * _iter_u);
            _dv = _iter_v * (K[4] * _r + K[5] * _r * _r + K[8] * _r * _r * _r) + 2 * K[7] * _iter_u * _iter_v + K[6] * (_r + 2 * _iter_v * _iter_v);
            //std::cout <<  "du:" <<_iter_dv << std::endl;
        }
        pixel[0] = _iter_u;
        pixel[1] = _iter_v;
    }

    void solver::cameraPoseFromHomography(const cv::Mat &H, cv::Mat &pose) {
        pose = cv::Mat::eye(3, 4, CV_64FC1); //3x4 matrix
        float norm1 = (float)cv::norm(H.col(0));
        float norm2 = (float)cv::norm(H.col(1));
        float tnorm = (norm1 + norm2) / 2.0f;

        cv::Mat v1 = H.col(0);
        cv::Mat v2 = pose.col(0);

        cv::normalize(v1, v2); // Normalize the rotation

        v1 = H.col(1);
        v2 = pose.col(1);

        cv::normalize(v1, v2);

        v1 = pose.col(0);
        v2 = pose.col(1);

        cv::Mat v3 = v1.cross(v2);  //Computes the cross-product of v1 and v2
        cv::Mat c2 = pose.col(2);
        v3.copyTo(c2);
        pose.col(3) = H.col(2) / tnorm; //vector t [R|t]
    }

    int solver::solve_camera_pose(std::string file_path,double(&Ceulur)[3],double(&Ctranslation)[3],double(&m)[3],double(&K)[9]) {
        //读取数据
        std::ifstream is(file_path, std::ios::in);
        if (!is.is_open()) {
            std::cout << "读取文件错误" << std::endl;
            return false;
        }
        //2D_3D 点对数量
        int N = 0;
        is >> N;
        //相机内参
        K[0] = 0;K[1] = 0;K[2] = 0;K[3] = 0;K[4] = 0;K[5] = 0;K[6] = 0;K[7] = 0;K[8] = 0;
        for (int i = 0; i < 9; i++) {
            is >> K[i];
        }
        //三维控制点坐标
        std::vector<double[3]> Device(N);
        //三维点对应的像素点
        std::vector<double[2]> pixel(N), pinit(N);
        //三维点坐标系转换
        for (int i = 0; i < N; i++) {
            is >> Device[i][0] >> Device[i][1] >> Device[i][2];
            //坐标纠正方向
            //double current_pose[2] = { 0.124976046,-0.06562438 };
            double X[3] = {
                    Device[i][0],
                    Device[i][1],
                    Device[i][2]
            };
            Device[i][2] = X[0];
            Device[i][0] = X[1];
            Device[i][1] = -X[2];
        }
//    cv::Mat camera_k = (cv::Mat_<char>(3, 3) << K[2], 0, K[0],
//            0, K[3], 0,
//            0, 0, 1);
//    cv::Mat camera_dis = (cv::Mat_<char>(5,1) << K[4],K[5],K[6],K[7],K[8] );
        //整张图片进行畸变纠正
        //cv::Mat input_img = cv::imread("D:/Desktop/quan/901/2/Camera90102.jpg");
        //int hei = input_img.rows;
        //int wid = input_img.cols;
        //cv::Mat out_img(hei, wid, CV_16FC3);
        //cv::undistort(input_img, out_img, camera_k, camera_dis);
        //cv::namedWindow("select", cv::WINDOW_NORMAL);
        //cv::resizeWindow("select", 992, 744);
        //cv::imshow("select", input_img);
        //cv::waitKey(0);
        //cv::imshow("select", out_img);
        //cv::waitKey(0);

        // 初始化求解库
        single_solver::solver solver;
        //像素点畸变纠正
        for (int i = 0; i < N; i++) {
            is >> pixel[i][0] >> pixel[i][1];
            pinit[i][0] = pixel[i][0];
            pinit[i][1] = pixel[i][1];
            double pix[2] = { };
            solver.undisrortion(pixel[i][0], pixel[i][1], K, pix);
            pixel[i][0] = pix[0];
            pixel[i][1] = pix[1];
            //std::cout << "my undistortion:" << pix[0] << " " << pix[1] << std::endl;
            //opencv_undisortion
            //std::vector<cv::Point2f> m1(1,cv::Point2f(pixel[i][0], pixel[i][1]));
            //std::vector<cv::Point2f> m2;
            //cv::undistortPoints(m1, m2,camera_k,camera_dis,cv::noArray(),camera_k);
            //std::cout << "opencv undistortion:" << m2 << std::endl;
        }
        is.close();

        //平面平差
        std::cout << "----------------- 平面平差 -----------------" << std::endl;
        //double m[3] = { 0,0,0 };
        //init m -- 平面方程系数
        m[0] = 0;
        m[1] = 0;
        m[2] = 0;
        ceres::Problem Cproblem_plane;
        for (int i = 0; i < Device.size(); i++)
        {
            //if (i != 4 && i != 6 && i != 11 && i != 14)
            //{
            //	continue;
            //}
            Cproblem_plane.AddResidualBlock(
                    new ceres::AutoDiffCostFunction<Plane_Fit, 1, 3>(
                            new Plane_Fit(Device[i])
                    ),
                    nullptr,
                    m
            );
        }
        ceres::Solver::Options options_plane;
        options_plane.linear_solver_type = ceres::DENSE_QR;
        options_plane.minimizer_progress_to_stdout = true;
        options_plane.max_num_iterations = 500;

        ceres::Solver::Summary summary_plane;
        std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
        ceres::Solve(options_plane, &Cproblem_plane, &summary_plane);
        std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();
        std::chrono::duration<double> time_used_plane = std::chrono::duration_cast<std::chrono::duration<double>>(t4 - t3);

        std::cout << "solve time cost = " << time_used_plane.count() << " seconds. " << std::endl;
        std::cout << "m0、m1、m2:            " << m[0] << " " << m[1] << " " << m[2] << std::endl;
        //ransac_p3p_kneip （选取四个点）
        Match2D3Ds corrs;
        for (int i = 0; i < N; i++)
        {
            if (i == 5 || i == 7 || i == 12 || i == 15)
            {
                Match2D3D corr;
                Eigen::Vector3d v3(Device[i][0], Device[i][1], Device[i][2]);
                corr.p3d = v3;
                Eigen::Vector2d v2(pixel[i][0], pixel[i][1]);
                corr.p2d = v2;
                corrs.push_back(corr);
            }

        }
        ResultP3p result;
        p3p_kneip::ransac_pose_p3p_kniep(corrs, &result, 1000, 0.005);

        std::cout << "2D-3D correspondences inliers: " << (100 * result.inliers.size() / corrs.size()) << std::endl;
        std::cout << "Estimated pose: " << std::endl;
        std::cout << result.pose << std::endl;

        Eigen::Matrix3d m33 = result.pose.block<3, 3>(0, 0);
        Eigen::Vector3d v33 = m33.eulerAngles(2, 1, 0);

        std::cout << "p3p_celur:\n" << v33(2) << " " << v33(1) << " " << v33(0) << std::endl;
        std::cout << "p3p_trans:\n" << result.pose(0, 3) << " " << result.pose(1, 3) << " " << result.pose(2, 3) << std::endl;

        //外参平差
        std::cout << "----------------- 外参平差 -----------------" << std::endl;
        Ceulur[0] = v33(2);
        Ceulur[1] = v33(1);
        Ceulur[2] = v33(0);

        Ctranslation[0] = result.pose(0,3);
        Ctranslation[1] = result.pose(1,3);
        Ctranslation[2] = result.pose(2,3);

        //double Ceulur[3] = {}; double Ctranslation[3] = { };
        ceres::Problem Cproblem;
        for (int i = 0; i < N; i++)
        {
            //if (i > 5 && i !=6 && i!=7 && i!= 14)
            //{
            //	continue;
            //}
            Cproblem.AddResidualBlock(
                    new ceres::AutoDiffCostFunction<Camera_Dial_Cost, 2, 3, 3>(
                            new Camera_Dial_Cost(pixel[i], Device[i])
                    ),
                    nullptr,
                    Ceulur, Ctranslation
            );
        }
        ceres::Solver::Options options;
        options.linear_solver_type = ceres::DENSE_QR;
        options.minimizer_progress_to_stdout = true;
        options.max_num_iterations = 500;

        ceres::Solver::Summary summary;
        std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
        ceres::Solve(options, &Cproblem, &summary);
        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        std::chrono::duration<double> time_used = std::chrono::duration_cast< std::chrono::duration<double>>(t2 - t1);

        std::cout << "solve time cost = " << time_used.count() << " seconds. " << std::endl;
        //cout << "Ceulur:            " << Ceulur[0]*180/3.1415926 << " " << Ceulur[1] * 180 / 3.1415926 << " " << Ceulur[2] * 180 / 3.1415926 << endl;
        std::cout << "Ceulur:            " << Ceulur[0] << " " << Ceulur[1] << " " << Ceulur[2] << std::endl;
        std::cout << "Ctranslation:      " << Ctranslation[0] << " " << Ctranslation[1] << " " << Ctranslation[2] << std::endl;

        //Eigen::Matrix3d m3 = euler_to_matrix(Ceulur);
        //std::cout << "ceres matrix :\n" << m3 << std::endl;

        //旋转矩阵转欧拉角
        //Eigen::Vector3d v3 = m3.eulerAngles(2, 1, 0);
        //std::cout << "p3p euler:\n" << v3 << std::endl;

        //Eigen::Matrix3d m33 = result.pose.block<3, 3>(0, 0);
        //Eigen::Vector3d v33 = m33.eulerAngles(2, 1, 0);
        //std::cout << "matrix3d:\n" << m33 << std::endl;
        //std::cout << "p3p euler:\n" << v33 << std::endl;

        //重投影误差
        std::cout << "----------------- 重投影误差 -----------------" << std::endl;
        double repp[2];
        double p1[2];
        double p2[2];
        std::vector<double[3]> p3(N);
        for (int i = 0; i < N; i++) {
            dist dist2 = solver.update_xy(Device[i], Ceulur, Ctranslation);
            repp[0] = dist2.dx;
            repp[1] = dist2.dy;
            double _r = repp[0] * repp[0] + repp[1] * repp[1];
            double xcor = repp[0] + repp[0] * (K[4] * _r + K[5] * _r * _r + K[8] * _r * _r * _r) + 2 * K[6] * repp[0] * repp[1] + K[7] * (_r + 2 * repp[0] * repp[0]);
            double ycor = repp[1] + repp[1] * (K[4] * _r + K[5] * _r * _r + K[8] * _r * _r * _r) + 2 * K[7] * repp[0] * repp[1] + K[6] * (_r + 2 * repp[1] * repp[1]);
            p2[0] = (xcor * K[2] + K[0]);
            p2[1] = (ycor * K[3] + K[1]);

            p1[0] = pinit[i][0];
            p1[1] = pinit[i][1];

            p3[i][0] = p2[0] - p1[0];
            p3[i][1] = p2[1] - p1[1];

            std::cout << "第" << (i + 1) << "个点的残差 dx:" << p3[i][0] << " dy:" << p3[i][1] << std::endl;
        }
        //solvept（二维到三维）
        std::cout << "----------------- 像素点到三维点（平面） -----------------" << std::endl;
        for (int i = 0; i < N; i++)
        {
            // double Ceulur[3] = { v33(2),v33(1),v33(0) };	double Ctranslation[3] = { result.pose(0,3),result.pose(1,3) ,result.pose(2,3) };
            double p3[3] = { };
            solver.point2d_to_point3d(pixel[i][0], pixel[i][1], Ceulur, Ctranslation, m, p3);
            double dx = (Device[i][0] - p3[0])*1000;
            double dy = (Device[i][1] - p3[1])*1000;
            double dz = (Device[i][2] - p3[2])*1000;

            std::cout << "第" << (i + 1) << "个点差值(mm)：dx: " << dx << " dy: " << dy << " dz:" << dz << std::endl;
            //std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(4) << std::sqrt(std::pow(dx, 2) + std::pow(dy, 2) + std::pow(dz, 2)) << std::endl;
        }
    }

    int solver::solve_3D_point_clouds(std::string image_path, std::string points3D_path,
                                      std::vector<cv::Point> &board_points,double(&K)[9],
                                      double(&Ceulur)[3],double(&Ctranslation)[3],double(&m)[3]){
        std::cout << "----------------- 像片生成三维点云 -----------------" << std::endl;
        cv::Mat img_c = imread(image_path, cv::IMREAD_COLOR);
        cv::Mat img_g = imread(image_path, cv::IMREAD_GRAYSCALE);
        if (img_c.empty())
        {
            std::cerr << "打开图片失败" << std::endl;
            return false;
        }
        //2D-3D自定义范围
        cv::Mat img_copy = img_g.clone();
        int	height = img_c.rows;
        int width = img_c.cols;
        int num_pixel = height * width;
        //int num_pixel = 4653600;
        unsigned char* image_data = new unsigned char[height * width * 3];
        memcpy(image_data, img_c.data, height * width * 3 * sizeof(unsigned char));

        int rgb[3] = {};
        int pixelx = 0;
        int pixely = 0;
        int num_3d = 0;

        std::ofstream ofs(points3D_path, std::ios::out);
        // ply
        //ofs << "ply" << std::endl;
        //ofs << "format ascii 1.0" << std::endl;
        //ofs << "element vertex " << num_pixel << std::endl;
        //ofs << "property float x" << std::endl;
        //ofs << "property float y" << std::endl;
        //ofs << "property float z" << std::endl;
        //ofs << "property uchar red" << std::endl;
        //ofs << "property uchar green" << std::endl;
        //ofs << "property uchar blue" << std::endl;
        //ofs << "end_header" << std::endl;

        //init single_solver
        single_solver::solver solver;
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
            {
                if (true)
                    //if (cv::pointPolygonTest(board_points, ptest, false) == 1)
                {
                    double pix[2] = { };
                    solver.undisrortion(j, i, K, pix);
                    double pic_x1 = pix[0];
                    double pic_y1 = pix[1];

                    double p3[3] = { };
                    solver.point2d_to_point3d(pic_x1, pic_y1, Ceulur, Ctranslation, m, p3);

                    rgb[0] = image_data[(i * width + j) * 3];
                    rgb[1] = image_data[(i * width + j) * 3 + 1];
                    rgb[2] = image_data[(i * width + j) * 3 + 2];

                    ofs << p3[0] << " " << p3[1] << " " << p3[2] << " " << (int)rgb[0] << " " << (int)rgb[1] << " " << (int)rgb[2] << " \n";
                    //cout << Xs << " " << Ys << " " << Zs << " " << (int)rgb[0] << " " << (int)rgb[1] << " " << (int)rgb[2] << " \n";
                    //std::cout << "第：" << num_3d + 1 << " 个点" << std::endl;
                    num_3d++;
                }

            }
        }
        std::cout << "写入数据完成!" << std::endl;
        std::cout << "共写入：" << num_3d << "个数据" << std::endl;

        return 0;
    }

    int solver::solve_Homography(std::string image1_path, std::string image2_path, cv::Mat& H) {
        std::cout << "----------------- 计算 Homography -----------------" << std::endl;
        cv::Mat image01 = cv::imread(image2_path, 1);    //右图
        cv::Mat image02 = cv::imread(image1_path, 1);    //左图
        //imshow("p2", image01);
        //imshow("p1", image02);

        if (image01.empty() || image01.empty())
        {
            std::cerr << "打开图片失败" << std::endl;
            return false;
        }
        //灰度图转换
        cv::Mat image1, image2;
        cvtColor(image01, image1, CV_RGB2GRAY);
        cvtColor(image02, image2, CV_RGB2GRAY);


        //提取特征点
        //cv::xfeatures2d::SurfFeatureDetector Detector = cv::xfeatures2d::SurfFeatureDetector::create(2000);
        cv::Ptr<cv::Feature2D> Detector = cv::xfeatures2d::SURF::create();

        std::vector<cv::KeyPoint> keyPoint1, keyPoint2;
        Detector->detect(image1, keyPoint1);
        Detector->detect(image2, keyPoint2);

        //特征点描述，为下边的特征点匹配做准备
        //SurfDescriptorExtractor Descriptor;
        cv::Mat imageDesc1, imageDesc2;
        Detector->compute(image1, keyPoint1, imageDesc1);
        Detector->compute(image2, keyPoint2, imageDesc2);

        cv::FlannBasedMatcher matcher;
        std::vector<std::vector<cv::DMatch> > matchePoints;
        std::vector<cv::DMatch> GoodMatchePoints;

        std::vector<cv::Mat> train_desc(1, imageDesc1);
        matcher.add(train_desc);
        matcher.train();

        matcher.knnMatch(imageDesc2, matchePoints, 2);
        std::cout << "total match points: " << matchePoints.size() << std::endl;

        // Lowe's algorithm,获取优秀匹配点
        for (int i = 0; i < matchePoints.size(); i++)
        {
            if (matchePoints[i][0].distance < 0.4 * matchePoints[i][1].distance)
            {
                GoodMatchePoints.push_back(matchePoints[i][0]);
            }
        }
        std::cout << "good match pointd is:" << GoodMatchePoints.size() << std::endl;
        cv::Mat first_match;
        drawMatches(image02, keyPoint2, image01, keyPoint1, GoodMatchePoints, first_match);
        cv::namedWindow("first_match ", cv::WINDOW_NORMAL);
        cv::resizeWindow("first_match ", 800, 600);
        imshow("first_match ", first_match);

        std::vector<cv::Point2f> imagePoints1, imagePoints2;

        for (int i = 0; i<GoodMatchePoints.size(); i++)
        {
            imagePoints2.push_back(keyPoint2[GoodMatchePoints[i].queryIdx].pt);
            imagePoints1.push_back(keyPoint1[GoodMatchePoints[i].trainIdx].pt);
        }

        //获取图像1到图像2的投影映射矩阵 尺寸为3*3
        cv::Mat homo = findHomography(imagePoints1, imagePoints2, CV_RANSAC);
        //也可以使用getPerspectiveTransform方法获得透视变换矩阵，不过要求只能有4个点，效果稍差
        //Mat   homo=getPerspectiveTransform(imagePoints1,imagePoints2);
        std::cout << "变换矩阵为：\n" << homo << std::endl; //输出映射矩阵
        H = homo;
        return 0;
    }

    int solver::solve_Essential_to_pose(std::string image1_path, std::string image2_path, Eigen::Matrix3d& R,Eigen::Vector3d& T,double(&K)[9]) {
        std::cout << "----------------- 计算 Essential to R T -----------------" << std::endl;
        cv::Mat image01 = cv::imread(image2_path, 1);    //右图
        cv::Mat image02 = cv::imread(image1_path, 1);    //左图

        if (image01.empty() || image01.empty())
        {
            std::cerr << "打开图片失败" << std::endl;
            return false;
        }
        //灰度图转换
        cv::Mat image1, image2;
        cvtColor(image01, image1, CV_RGB2GRAY);
        cvtColor(image02, image2, CV_RGB2GRAY);

        //提取特征点
        //cv::xfeatures2d::SurfFeatureDetector Detector = cv::xfeatures2d::SurfFeatureDetector::create(2000);
        cv::Ptr<cv::Feature2D> Detector = cv::xfeatures2d::SURF::create();

        std::vector<cv::KeyPoint> keyPoint1, keyPoint2;
        Detector->detect(image1, keyPoint1);
        Detector->detect(image2, keyPoint2);

        //特征点描述，为下边的特征点匹配做准备
        //SurfDescriptorExtractor Descriptor;
        cv::Mat imageDesc1, imageDesc2;
        Detector->compute(image1, keyPoint1, imageDesc1);
        Detector->compute(image2, keyPoint2, imageDesc2);

        cv::FlannBasedMatcher matcher;
        std::vector<std::vector<cv::DMatch> > matchePoints;
        std::vector<cv::DMatch> GoodMatchePoints;

        std::vector<cv::Mat> train_desc(1, imageDesc1);
        matcher.add(train_desc);
        matcher.train();

        matcher.knnMatch(imageDesc2, matchePoints, 2);
        std::cout << "total match points: " << matchePoints.size() << std::endl;

        // Lowe's algorithm,获取优秀匹配点
        for (int i = 0; i < matchePoints.size(); i++)
        {
            if (matchePoints[i][0].distance < 0.4 * matchePoints[i][1].distance)
            {
                GoodMatchePoints.push_back(matchePoints[i][0]);
            }
        }
        std::cout << "good match pointd is:" << GoodMatchePoints.size() << std::endl;
        //cv::Mat first_match;
        //drawMatches(image02, keyPoint2, image01, keyPoint1, GoodMatchePoints, first_match);
        //cv::namedWindow("first_match ", cv::WINDOW_NORMAL);
        //cv::resizeWindow("first_match ", 800, 600);
        //imshow("first_match ", first_match);

        std::vector<cv::Point2f> imagePoints1, imagePoints2;

        for (int i = 0; i<GoodMatchePoints.size(); i++)
        {
            imagePoints2.push_back(keyPoint2[GoodMatchePoints[i].queryIdx].pt);
            imagePoints1.push_back(keyPoint1[GoodMatchePoints[i].trainIdx].pt);
        }

        //获取图像1到图像2的投影映射矩阵 尺寸为3*3
        cv::Mat homo = findHomography(imagePoints1, imagePoints2, CV_RANSAC);
        cv::Mat E;
        cv::Mat r;
        cv::Mat t;
        cv::Mat mask;
        double focal = K[0];
        cv::Point2d pp(K[2],K[5]);
        E = findEssentialMat(imagePoints1, imagePoints2, focal, pp, cv::RANSAC, 0.999, 1.0, mask);

        //Find Pright camera matrix from the essential matrix
        //Cheirality check (all points are in front of camera) is performed internally.
        cv::recoverPose(E, imagePoints1, imagePoints2, r, t, focal, pp, mask);

        //TODO: stratify over Pleft
        R << r.at<double>(0,0), r.at<double>(0,1), r.at<double>(0,2),
             r.at<double>(1,0), r.at<double>(1,1), r.at<double>(1,2),
             r.at<double>(2,0), r.at<double>(2,1), r.at<double>(2,2);
        T << t.at<double>(0),
             t.at<double>(1),
             t.at<double>(2);
        return 0;
    }

    int solver::solve_pose_from_opencv(std::string image1_path,std::string image2_path, cv::Mat &K, std::vector<cv::Mat> &Rs,
                                       std::vector<cv::Mat> &Ts,Eigen::Matrix3f& RR,Eigen::Vector3f& TT) {
        std::cout << "----------------- 计算 Homography -----------------" << std::endl;
        cv::Mat image01 = cv::imread(image2_path, 1);    //右图
        cv::Mat image02 = cv::imread(image1_path, 1);    //左图
        //imshow("p2", image01);
        //imshow("p1", image02);

        if (image01.empty() || image01.empty())
        {
            std::cerr << "打开图片失败" << std::endl;
            return false;
        }
        //灰度图转换
        cv::Mat image1, image2;
        cvtColor(image01, image1, CV_RGB2GRAY);
        cvtColor(image02, image2, CV_RGB2GRAY);


        //提取特征点
        //cv::xfeatures2d::SurfFeatureDetector Detector = cv::xfeatures2d::SurfFeatureDetector::create(2000);
        cv::Ptr<cv::Feature2D> Detector = cv::xfeatures2d::SURF::create();

        std::vector<cv::KeyPoint> keyPoint1, keyPoint2;
        Detector->detect(image1, keyPoint1);
        Detector->detect(image2, keyPoint2);

        //特征点描述，为下边的特征点匹配做准备
        //SurfDescriptorExtractor Descriptor;
        cv::Mat imageDesc1, imageDesc2;
        Detector->compute(image1, keyPoint1, imageDesc1);
        Detector->compute(image2, keyPoint2, imageDesc2);

        cv::FlannBasedMatcher matcher;
        std::vector<std::vector<cv::DMatch> > matchePoints;
        std::vector<cv::DMatch> GoodMatchePoints;

        std::vector<cv::Mat> train_desc(1, imageDesc1);
        matcher.add(train_desc);
        matcher.train();

        matcher.knnMatch(imageDesc2, matchePoints, 2);
        std::cout << "total match points: " << matchePoints.size() << std::endl;

        // Lowe's algorithm,获取优秀匹配点
        for (int i = 0; i < matchePoints.size(); i++)
        {
            if (matchePoints[i][0].distance < 0.4 * matchePoints[i][1].distance)
            {
                GoodMatchePoints.push_back(matchePoints[i][0]);
            }
        }
        std::cout << "good match pointd is:" << GoodMatchePoints.size() << std::endl;
        //cv::Mat first_match;
        //drawMatches(image02, keyPoint2, image01, keyPoint1, GoodMatchePoints, first_match);
        //cv::namedWindow("first_match ", cv::WINDOW_NORMAL);
        //cv::resizeWindow("first_match ", 800, 600);
        //imshow("first_match ", first_match);

        std::vector<cv::Point2f> imagePoints1, imagePoints2;

        for (int i = 0; i<GoodMatchePoints.size(); i++)
        {
            imagePoints2.push_back(keyPoint2[GoodMatchePoints[i].queryIdx].pt);
            imagePoints1.push_back(keyPoint1[GoodMatchePoints[i].trainIdx].pt);
        }

        //获取图像1到图像2的投影映射矩阵 尺寸为3*3
        cv::Mat homo = findHomography(imagePoints1, imagePoints2, CV_RANSAC);
        std::cout << "H: \n" << homo << std::endl;
        //homography 分解
        int results = cv::decomposeHomographyMat(homo,K,Rs,Ts,cv::noArray());
        for (int i = 0; i < results; ++i) {
            cv::Mat R = Rs.at(i).t();
            cv::Mat t = Ts.at(i).t();

            cv::Matx34f pose1 = cv::Matx34f::eye();
            cv::Matx34f pose2 = cv::Matx34f(R.at<double>(0,0), R.at<double>(0,1), R.at<double>(0,2), t.at<double>(0),
                             R.at<double>(1,0), R.at<double>(1,1), R.at<double>(1,2), t.at<double>(1),
                             R.at<double>(2,0), R.at<double>(2,1), R.at<double>(2,2), t.at<double>(2));
            std::vector<cv::Point2f> p1,p2;
            p1.push_back(imagePoints1[0]);
            p2.push_back(imagePoints2[0]);
            cv::Mat point3D;
            cv::triangulatePoints(pose1,pose2,p2,p1,point3D);
            cv::Mat rvecLeft;
            cv:Rodrigues(pose2.get_minor<3, 3>(0, 0), rvecLeft);
            cv::Mat tvecLeft(pose2.get_minor<3, 1>(0, 3).t());
            Eigen::Matrix3f eR;
            eR << R.at<double>(0,0), R.at<double>(0,1), R.at<double>(0,2),
                  R.at<double>(1,0), R.at<double>(1,1), R.at<double>(1,2),
                  R.at<double>(2,0), R.at<double>(2,1), R.at<double>(2,2);
                    //cv::cv2eigen(rvecLeft,eR);
            Eigen::Vector3f eT;
            eT << t.at<double>(0),
                    t.at<double>(1),
                    t.at<double>(2);
            Eigen::Vector3f p3(point3D.at<float>(0,0),point3D.at<float>(0,1),point3D.at<float>(0,2));
            Eigen::Vector3f px1 = eR * p3 + eT;
            //std::cout << "px1: \n" << px1 << std::endl;

            if (px1[2] > 0 && p3[2] > 0)
            {
                std::cout << "find right R T" <<std::endl;
                std::cout << "R:" << "\n" << R <<std::endl;
                std::cout << "T:" << "\n" << t <<std::endl;
                RR = eR;
                TT = eT;
                break;
            }
            //std::vector<cv::Point2f> projectedOnLeft;
            //cv::projectPoints(point3D, rvecLeft, tvecLeft, K, cv::Mat(), projectedOnLeft);
            //std::cout << (i+1) << " point2d: " << projectedOnLeft[0].x << " " << projectedOnLeft[0].y << std::endl;
        }
        return 0;
    }

}