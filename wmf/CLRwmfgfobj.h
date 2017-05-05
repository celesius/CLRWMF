/*
 * CLRwmfgfobj.h
 *
 *  Created on: 2017年5月4日
 *      Author: clover
 */

#ifndef WMF_CLRWMFGFOBJ_H_
#define WMF_CLRWMFGFOBJ_H_

#include <opencv2/opencv.hpp>

/*
    I: [288x384x3 double]
        r: 10
      eps: 1.0000e-04
        N: [288x384 double]
 mean_I_r: [288x384 double]
 mean_I_g: [288x384 double]
 mean_I_b: [288x384 double]
 var_I_rr: [288x384 double]
 var_I_rg: [288x384 double]
 var_I_rb: [288x384 double]
 var_I_gg: [288x384 double]
 var_I_gb: [288x384 double]
 var_I_bb: [288x384 double]
 invSigma: {288x384 cell}
*/

class CLRwmf_gfobj {
public:
	CLRwmf_gfobj();
	virtual ~CLRwmf_gfobj();

	cv::Mat m_image;
	int m_r;
	double m_eps;
	cv::Mat m_N;
	cv::Mat m_mean_I_r;
	cv::Mat m_mean_I_g;
	cv::Mat m_mean_I_b;
	cv::Mat m_var_I_rr;
	cv::Mat m_var_I_rg;
	cv::Mat m_var_I_rb;
	cv::Mat m_var_I_gg;
	cv::Mat m_var_I_gb;
	cv::Mat m_var_I_bb;
	//cv::Mat m_invSigma;
	std::vector<cv::Mat> m_invSigma;
};

#endif /* WMF_CLRWMFGFOBJ_H_ */
