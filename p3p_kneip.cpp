#include "p3p_kneip.h"
#define Eigen_EPSILON_EQ(x,v,eps) (((v - eps) <= x) && (x <= (v + eps)))
namespace p3p_kneip {
	void solve_quartic_roots(Eigen::Matrix<double, 5, 1> const& factors, Eigen::Matrix<double, 4, 1>* real_roots)
	{
		double const A = factors[0];
		double const B = factors[1];
		double const C = factors[2];
		double const D = factors[3];
		double const E = factors[4];

		double const A2 = A * A;
		double const B2 = B * B;
		double const A3 = A2 * A;
		double const B3 = B2 * B;
		double const A4 = A3 * A;
		double const B4 = B3 * B;

		double const alpha = -3.0 * B2 / (8.0 * A2) + C / A;
		double const beta = B3 / (8.0 * A3) - B * C / (2.0 * A2) + D / A;
		double const gamma = -3.0 * B4 / (256.0 * A4) + B2 * C / (16.0 * A3) - B * D / (4.0 * A2) + E / A;

		double const alpha2 = alpha * alpha;
		double const alpha3 = alpha2 * alpha;
		double const beta2 = beta * beta;

		std::complex<double> P(-alpha2 / 12.0 - gamma, 0.0);
		std::complex<double> Q(-alpha3 / 108.0 + alpha * gamma / 3.0 - beta2 / 8.0, 0.0);
		std::complex<double> R = -Q / 2.0 + std::sqrt(Q * Q / 4.0 + P * P * P / 27.0);

		std::complex<double> U = std::pow(R, 1.0 / 3.0);
		std::complex<double> y = (U.real() == 0.0)
			? -5.0 * alpha / 6.0 - std::pow(Q, 1.0 / 3.0)
			: -5.0 * alpha / 6.0 - P / (3.0 * U) + U;

		std::complex<double> w = std::sqrt(alpha + 2.0 * y);
		std::complex<double> part1 = -B / (4.0 * A);
		std::complex<double> part2 = 3.0 * alpha + 2.0 * y;
		std::complex<double> part3 = 2.0 * beta / w;

		std::complex<double> complex_roots[4];
		complex_roots[0] = part1 + 0.5 * (w + std::sqrt(-(part2 + part3)));
		complex_roots[1] = part1 + 0.5 * (w - std::sqrt(-(part2 + part3)));
		complex_roots[2] = part1 + 0.5 * (-w + std::sqrt(-(part2 - part3)));
		complex_roots[3] = part1 + 0.5 * (-w - std::sqrt(-(part2 - part3)));

		for (int i = 0; i < 4; ++i)
			(*real_roots)[i] = complex_roots[i].real();
	}
	void pose_p3p_kneip(Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3, Eigen::Vector3d f1, Eigen::Vector3d f2, Eigen::Vector3d f3, std::vector<Eigen::Matrix<double, 3, 4>>* solutions)
	{
		/* Check if points are co-linear. In this case return no solution. */
		double const colinear_threshold = 1e-10;
		if ((p2 - p1).cross(p3 - p1).norm() < colinear_threshold) {
			solutions->clear();
			return;
		}

		/* Normalize directions if necessary. */
		double const normalize_epsilon = 1e-10;
		if (!Eigen_EPSILON_EQ(f1.norm(), 1.0, normalize_epsilon))
			f1.normalize();
		if (!Eigen_EPSILON_EQ(f2.norm(), 1.0, normalize_epsilon))
			f2.normalize();
		if (!Eigen_EPSILON_EQ(f3.norm(), 1.0, normalize_epsilon))
			f3.normalize();
		//std::cout << "f1:\n" << f1 << std::endl;
		//std::cout << "f2:\n" << f2 << std::endl;
		//std::cout << "f3:\n" << f3 << std::endl;
		/* Create camera frame. */
		Eigen::Matrix3d T;
		{
			Eigen::Vector3d e1 = f1;
			Eigen::Vector3d e3 = (f1.cross(f2)).normalized();
			Eigen::Vector3d e2 = e3.cross(e1);//?
			T.row(0) = e1;
			T.row(1) = e2;
			T.row(2) = e3;
			//std::copy(e1.begin(), e1.end(), T.begin() + 0);
			//std::copy(e2.begin(), e2.end(), T.begin() + 3);
			//std::copy(e3.begin(), e3.end(), T.begin() + 6);
			f3 = T * f3;
			//std::cout << "e1:\n" << e1 << std::endl;
			//std::cout << "e2:\n" << e2 << std::endl;
			//std::cout << "e3:\n" << e3 << std::endl;

			//std::cout << "f3:\n" << f3 << std::endl;
		}

		/* Change camera frame and point order if f3[2] > 0. */
		if (f3[2] > 0.0)
		{
			std::swap(p1, p2);
			std::swap(f1, f2);

			Eigen::Vector3d e1 = f1;
			Eigen::Vector3d e3 = f1.cross(f2).normalized();
			Eigen::Vector3d e2 = e3.cross(e1);
			T.row(0) = e1;
			T.row(1) = e2;
			T.row(2) = e3;
			//std::copy(e1.begin(), e1.end(), T.begin() + 0);
			//std::copy(e2.begin(), e2.end(), T.begin() + 3);
			//std::copy(e3.begin(), e3.end(), T.begin() + 6);
			f3 = T * f3;
		}

		/* Create world frame. */
		Eigen::Matrix3d N;
		{
			Eigen::Vector3d n1 = (p2 - p1).normalized();
			Eigen::Vector3d n3 = n1.cross(p3 - p1).normalized();
			Eigen::Vector3d n2 = n3.cross(n1);
			N.row(0) = n1;
			N.row(1) = n2;
			N.row(2) = n3;
	/*		std::copy(n1.begin(), n1.end(), N.begin() + 0);
			std::copy(n2.begin(), n2.end(), N.begin() + 3);
			std::copy(n3.begin(), n3.end(), N.begin() + 6);*/
		}
		p3 = N * (p3 - p1);
		//std::cout << "f1:" << f1 << std::endl;
		//std::cout << "f2:" << f2 << std::endl;
		//std::cout << "f3:" << f3 << std::endl;
		/* Extraction of known parameters. */
		double d_12 = (p2 - p1).norm();
		double f_1 = f3[0] / f3[2];
		//std::cout << "f_1:" << f_1 << std::endl;
		double f_2 = f3[1] / f3[2];
		double p_1 = p3[0];
		double p_2 = p3[1];

		double cos_beta = f1.dot(f2);
		//double b = 1.0 / (1.0 - std::pow(cos_beta)) - 1;
		double b = 1.0 / (1.0 - std::pow(cos_beta,2)) - 1;
		////std::cout << "cos_beta:" << b << std::endl;


		if (cos_beta < 0.0)
			b = -std::sqrt(b);
		else
			b = std::sqrt(b);

		/* Temporary precomputed variables. */
		double f_1_pw2 = std::pow(f_1,2);
		//std::cout << "f_1_pw2:" << f_1_pw2 << std::endl;
		double f_2_pw2 = std::pow(f_2,2);
		double p_1_pw2 = std::pow(p_1,2);
		double p_1_pw3 = p_1_pw2 * p_1;
		double p_1_pw4 = p_1_pw3 * p_1;
		double p_2_pw2 = std::pow(p_2,2);
		double p_2_pw3 = p_2_pw2 * p_2;
		double p_2_pw4 = p_2_pw3 * p_2;
		double d_12_pw2 = std::pow(d_12,2);
		double b_pw2 = std::pow(b,2);

		/* Factors of the 4th degree polynomial. */
		Vec5d factors;
		factors[0] = -f_2_pw2 * p_2_pw4 - p_2_pw4 * f_1_pw2 - p_2_pw4;

		factors[1] = 2.0 * p_2_pw3 * d_12 * b
			+ 2.0 * f_2_pw2 * p_2_pw3 * d_12 * b
			- 2.0 * f_2 * p_2_pw3 * f_1 * d_12;

		factors[2] = -f_2_pw2 * p_2_pw2 * p_1_pw2
			- f_2_pw2 * p_2_pw2 * d_12_pw2 * b_pw2
			- f_2_pw2 * p_2_pw2 * d_12_pw2
			+ f_2_pw2 * p_2_pw4
			+ p_2_pw4 * f_1_pw2
			+ 2.0 * p_1 * p_2_pw2 * d_12
			+ 2.0 * f_1 * f_2 * p_1 * p_2_pw2 * d_12 * b
			- p_2_pw2 * p_1_pw2 * f_1_pw2
			+ 2.0 * p_1 * p_2_pw2 * f_2_pw2 * d_12
			- p_2_pw2 * d_12_pw2 * b_pw2
			- 2.0 * p_1_pw2 * p_2_pw2;

		factors[3] = 2.0 * p_1_pw2 * p_2 * d_12 * b
			+ 2.0 * f_2 * p_2_pw3 * f_1 * d_12
			- 2.0 * f_2_pw2 * p_2_pw3 * d_12 * b
			- 2.0 * p_1 * p_2 * d_12_pw2 * b;

		factors[4] = -2.0 * f_2 * p_2_pw2 * f_1 * p_1 * d_12 * b
			+ f_2_pw2 * p_2_pw2 * d_12_pw2
			+ 2.0 * p_1_pw3 * d_12
			- p_1_pw2 * d_12_pw2
			+ f_2_pw2 * p_2_pw2 * p_1_pw2
			- p_1_pw4
			- 2.0 * f_2_pw2 * p_2_pw2 * p_1 * d_12
			+ p_2_pw2 * f_1_pw2 * p_1_pw2
			+ f_2_pw2 * p_2_pw2 * d_12_pw2 * b_pw2;

		//std::cout << "factors[0]:" << factors[0] << " factors[1]:" << factors[1] << " factors[2]:" << factors[2]
			//<< " factors[3]:" << factors[3] << " factors[4]:" << factors[4] << std::endl;
		//factors[0]:-203.927 factors[1] : 383.494 factors[2] : -157.803 factors[3] : -54.3482 factors[4] : 24.1584

		/* Solve for the roots of the polynomial. */
		Vec4d real_roots;
		solve_quartic_roots(factors, &real_roots);

		//std::cout << "real_roots:" << real_roots << std::endl;

		/* Back-substitution of each solution. */
		solutions->clear();
		solutions->resize(4);
		for (int i = 0; i < 4; ++i)
		{
			double cot_alpha = (-f_1 * p_1 / f_2 - real_roots[i] * p_2 + d_12 * b)
				/ (-f_1 * real_roots[i] * p_2 / f_2 + p_1 - d_12);

			double cos_theta = real_roots[i];
			//std::cout << "cos_theta:" <<  cos_theta << std::endl;
			double sin_theta = std::sqrt(1.0 - std::pow(real_roots[i],2));
			double sin_alpha = std::sqrt(1.0 / (std::pow(cot_alpha,2) + 1));
			double cos_alpha = std::sqrt(1.0 - std::pow(sin_alpha,2));

			if (cot_alpha < 0.0)
				cos_alpha = -cos_alpha;

			Eigen::Vector3d C(
				d_12 * cos_alpha * (sin_alpha * b + cos_alpha),
				cos_theta * d_12 * sin_alpha * (sin_alpha * b + cos_alpha),
				sin_theta * d_12 * sin_alpha * (sin_alpha * b + cos_alpha));

			C = p1 + N.transpose() * C;

			Eigen::Matrix3d R;
			R(0,0) = -cos_alpha; R(0,1) = -sin_alpha * cos_theta; R(0,2) = -sin_alpha * sin_theta;
			R(1,0) = sin_alpha;  R(1,1) = -cos_alpha * cos_theta; R(1,2) = -cos_alpha * sin_theta;
			R(2,0) = 0.0;        R(2,1) = -sin_theta;             R(2,2) = cos_theta;
			R = N.transpose()*(R.transpose())*(T);
			//std::cout << "R:" << R << std::endl;

			/* Convert camera position and cam-to-world rotation to pose. */

			//R = R.transpose();
			R.transposeInPlace();

			//std::cout << "RT:" << R << std::endl;

			C = -R * C;

			solutions->at(i) << R(0, 0), R(0, 1), R(0, 2), C(0),
				R(1, 0), R(1, 1), R(1, 2), C(1),
				R(2, 0), R(2, 1), R(2, 2), C(2);


			//solutions->at(i) = R.hstack(C);
		}
	}
	void compute_p3p(Match2D3Ds const& corresp, PutativePoses* poses)
	{
		if (corresp.size() < 3)
			throw std::invalid_argument("At least 3 correspondences required");

		/* Draw 3 unique random numbers. */
		std::set<int> result;
		while (result.size() < 3)
			result.insert(std::rand() % corresp.size());

		std::set<int>::const_iterator iter = result.begin();
		Match2D3D const& c1(corresp[*iter++]);
		Match2D3D const& c2(corresp[*iter++]);
		Match2D3D const& c3(corresp[*iter]);

		//std::cout << "c1:" << c1.p3d << " c2:" << c2.p3d << std::endl;



		//
		pose_p3p_kneip(
			Eigen::Vector3d(c1.p3d), Eigen::Vector3d(c2.p3d), Eigen::Vector3d(c3.p3d),
			(Eigen::Vector3d(c1.p2d[0], c1.p2d[1], 1.0)),
			(Eigen::Vector3d(c2.p2d[0], c2.p2d[1], 1.0)),
			(Eigen::Vector3d(c3.p2d[0], c3.p2d[1], 1.0)),
			poses);
	}

	void ransac_pose_p3p_kniep(Match2D3Ds const& corresp, ResultP3p* result, int max_iterations, double threshold)
	{
		//std::atomic<int> num_iterations;
#pragma omp parallel
		{
			std::vector<int> inliers;
			inliers.reserve(corresp.size());
#pragma omp for
			for (int i = 0; i < max_iterations; i++)
			{
				PutativePoses poses;
				compute_p3p(corresp, &poses);
				/* Check all putative solutions and count inliers. */
				for (int j = 0; j < poses.size(); j++)
				{
					find_inliers(corresp, poses[j], &inliers,threshold);
#pragma omp critical
					if (inliers.size() > result->inliers.size())
					{
						result->pose = poses[j];
						std::swap(result->inliers, inliers);
						inliers.reserve(corresp.size());
					}
				}
			}
		}
	}

	void find_inliers(Match2D3Ds const& corresp, pose const& pose, std::vector<int>* inliers, double threshold)
	{
		inliers->resize(0);
		double const square_threshold = std::pow(threshold, 2);
		for (int i = 0; i < corresp.size(); i++)
		{
			Match2D3D const c = corresp[i];
			Vec4d p3d(c.p3d(0), c.p3d(1), c.p3d(2),1);
			Eigen::Vector3d p2d = pose * p3d;
			double square_error = std::pow((p2d(0) / p2d(2) - c.p2d(0)), 2) +
				std::pow((p2d(1) / p2d(2) - c.p2d(1)), 2);
			if (square_error < square_threshold)
			{
				inliers->push_back(i);
			}
		}
	}

}