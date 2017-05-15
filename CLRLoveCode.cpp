//============================================================================
// Name        : CLRLoveCode.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <memory>
#include <unistd.h>
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

struct clr_struct{
	int a;
	int b;
	int c;
	int d;
	int e = 15;
	struct clr_struct* cc;
};

typedef clr_struct CLRStruct ;

int main() {
	//cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!
	//unique_ptr<CLRClass> clr_ptr(new CLRClass(10));
	cv::Mat blob = cv::imread("/home/clover/ifvrftp/usr/jiangbo/BinoStereoImageDataset/2017.3.28THU/StereoDataset/blob/image182.jpg");
	cv::resize(blob, blob, cv::Size(720, 404));
	cv::imwrite("blob.jpg",blob);

	double d_3[3] = {10,2,3};
	double d_33[9] = {1,2,3,4,5,6,7,8,9};
	cv::Mat mat_3 = cv::Mat(1,3, CV_64FC1, d_3);
	cv::Mat mat_33 = cv::Mat(3,3, CV_64FC1, d_33);

	std::cout << mat_3 * 3<< std::endl;
	std::cout << mat_3.ptr<double>(0)[0] <<std::endl;
	std::cout << mat_3.ptr<double>(0)[1] <<std::endl;
	std::cout << mat_3.ptr<double>(0)[2] <<std::endl;


	/*
	printf("float %d\n", sizeof(float*));
	printf("%d\n", sizeof(int*));
	printf("%d\n", sizeof(void*));
	printf("%d\n", sizeof(double*));
	printf("--\n");
	printf("%d\n", sizeof(int));
	printf("%d\n", sizeof(CLRStruct));
	printf("%d\n", sizeof(CLRStruct*));

	CLRStruct* aaa = (CLRStruct *)malloc(sizeof(CLRStruct));
	aaa->cc = aaa;
	aaa->a = 1;
	aaa->b = 2;
	aaa->c = 3;
	aaa->d = 4;
	aaa->e = 15;
	printf("%d\n", sizeof(aaa));
	printf("d = %d\n", aaa->cc->d);
	aaa->d = 3;
	printf("d = %d\n", aaa->cc->d);
	//free(aaa->cc);

	CLRStruct* bbb = (CLRStruct *)malloc(sizeof(CLRStruct));
	bbb->cc = bbb;
	bbb->a = 4;
	bbb->b = 3;
	bbb->c = 2;
	bbb->d = 1;
	bbb->e = 15;
*/

	//cv::Mat eye_mat = cv::Mat::eye(cv::Size(3,3), CV_64FC1);
	//std::cout << eye_mat << std::endl;

	cv::Mat image = cv::imread("/home/clover/workspace/Matlab/wmf/matlab/matlab_wmf_release_v1/img_stereo/tsukuba_left.png");
	cv::Mat disp = cv::imread("/home/clover/workspace/Matlab/wmf/matlab/matlab_wmf_release_v1/img_stereo/tsukuba_boxagg.png", -1);

	unique_ptr<CLRwmf> clr_ptr(new CLRwmf());
	double t = cv::getTickCount();
	cv::Mat out = clr_ptr->run_wmf(image, disp);
	t = cv::getTickCount() - t;
	t = t/(cv::getTickFrequency());
	printf("t = %f\n",t);

	t = cv::getTickCount();
	usleep(1*1000*1000);
	//cv::waitKey(1*1000);
	t = cv::getTickCount() - t;
	t = t/(cv::getTickFrequency());
	printf("t = %f\n",t);


	cv::imshow("result", out);
	cv::waitKey(0);
	/*
	 */

	std::cout<< " end " <<std::endl;
	return 0;
}
