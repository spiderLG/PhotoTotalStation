#pragma once
#include <iostream>
#include <set>
#include <vector>
#include <omp.h>

#include <Eigen/Eigen>

struct Match2D3D
{
	Eigen::Vector3d p3d;
	Eigen::Vector2d p2d;
};
typedef std::vector<Match2D3D> Match2D3Ds;

typedef Eigen::Matrix<double, 3, 4> pose;
typedef std::vector<pose> PutativePoses;
typedef Eigen::Matrix<double, 5, 1> Vec5d;
typedef Eigen::Matrix<double, 4, 1> Vec4d;

struct ResultP3p
{
	Eigen::Matrix<double, 3, 4> pose;
	std::vector<int> inliers;
};

namespace p3p_kneip
{
	//�����cos()���Ĵη���
	void solve_quartic_roots(Eigen::Matrix<double, 5, 1> const& factors, Eigen::Matrix<double, 4, 1>* real_roots);
	//solve_p3p_kneip
	void pose_p3p_kneip(Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3, 
		Eigen::Vector3d f1, Eigen::Vector3d f2, Eigen::Vector3d f3,std::vector<Eigen::Matrix<double,3,4>>* solutions);
	void compute_p3p(Match2D3Ds const& corresp, PutativePoses* poses);
	//ransac_p3p_kniep
	void ransac_pose_p3p_kniep(Match2D3Ds const& corresp, ResultP3p* result,int max_iterations,double threshold);
	void find_inliers(Match2D3Ds const& corresp, pose const& pose, std::vector<int>* inliers,double threshold);
}