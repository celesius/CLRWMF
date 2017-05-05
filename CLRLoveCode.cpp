//============================================================================
// Name        : CLRLoveCode.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <memory>
#include "CLRClass.h"
#include <opencv2/opencv.hpp>
#include "wmf/CLRwmf.h"
using namespace std;

void printff(unique_ptr<CLRClass> &cls){
	cls->print_para();
}

void local(){
	unique_ptr<CLRClass> clr_ptr(new CLRClass(20));
	//printff(clr_ptr);
	CLRClass* clr_ptr2 = new CLRClass(30);
}


int main() {
	//cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!
	//unique_ptr<CLRClass> clr_ptr(new CLRClass(10));

	uchar data[] = {1,2,3,
				  4,5,6,
				  7,8,9};

	cv::Mat _mat = cv::Mat(9,9,CV_8UC1,cv::Scalar::all(1));
	cv::Mat img = cv::imread("/home/clover/workspace/Matlab/wmf/matlab/matlab_wmf_release_v1/img_stereo/tsukuba_left.png");
	cv::Mat disp = cv::imread("/home/clover/workspace/Matlab/wmf/matlab/matlab_wmf_release_v1/img_stereo/tsukuba_boxagg.png");
	cv::cvtColor(disp, disp, CV_BGR2GRAY);
	//cv::boxFilter(_mat,_mat,-1,cv::Size(3,3));
	//std::cout<< _mat << std::endl;


	std::cout<< " start " <<std::endl;
	unique_ptr<CLRwmf> clr_ptr(new CLRwmf());
	clr_ptr->run_wmf(img, disp);
	//cv::imshow("disp", disp);
	//cv::waitKey(0);

	cv::Mat data_mat = cv::Mat(3,3,CV_8UC1, data);
	cv::Mat data2_mat = cv::Mat(3,3,CV_8UC1, data);
	cv::multiply(data_mat, data2_mat, data_mat);
	cv::multiply(data_mat, data2_mat, data_mat);

	std::cout << data_mat << std::endl;

//	double	data_db[] = {0.1,0.2,0.3,
//				  0.4,0.5,0.6,
//				  0.7,0.8,0.9};
	/*
	double	data_db[] = {2,0,0,
				  0,1,3,
				  10,8,1};

	cv::Mat db = cv::Mat(3,3,CV_64FC1, data_db);
	cv::Mat db_inv ;//= db.inv();
	cv::invert(db, db_inv);
	cv::Mat mult_db = db.mul(db_inv);
	//cv::multiply(db, db_inv, mult_db);
	std::cout << db << std::endl;
	std::cout << db_inv << std::endl;
	std::cout << db*db_inv << std::endl;
*/
	double d_3[3] = {10,2,3};
	double d_33[9] = {1,2,3,4,5,6,7,8,9};
	cv::Mat mat_3 = cv::Mat(1,3, CV_64FC1, d_3);
	cv::Mat mat_33 = cv::Mat(3,3, CV_64FC1, d_33);

	std::cout << mat_3 * 3<< std::endl;
	std::cout << mat_3.ptr<double>(0)[0] <<std::endl;
	std::cout << mat_3.ptr<double>(0)[1] <<std::endl;
	std::cout << mat_3.ptr<double>(0)[2] <<std::endl;

	//cv::Mat eye_mat = cv::Mat::eye(cv::Size(3,3), CV_64FC1);
	//std::cout << eye_mat << std::endl;

	//clr_ptr->run_wmf(_mat);

	/*
	cv::Mat data_mat = cv::Mat(3,3,CV_8UC1, data);
	std::cout<< data_mat <<std::endl;
	cv::Mat double_mat;
	data_mat.convertTo(double_mat, CV_64F, 1.0/255.0 );
	std::cout<< double_mat <<std::endl;
*/

	std::cout<< " end " <<std::endl;
	return 0;
}
