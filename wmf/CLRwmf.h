/*
 * CLRwmf.h
 *
 *  Created on: 2017年5月3日
 *      Author: clover
 */

#ifndef WMF_CLRWMF_H_
#define WMF_CLRWMF_H_
#include <opencv2/opencv.hpp>
#include "CLRwmfgfobj.h"


class CLRwmf {
	CLRwmf_gfobj* m_gfobj;

public:
	CLRwmf();
	virtual ~CLRwmf();
	void run_wmf(const cv::Mat& src_mat,const cv::Mat& disp);
};

#endif /* WMF_CLRWMF_H_ */
