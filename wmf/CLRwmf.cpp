/*
 * CLRwmf.cpp
 *
 *  Created on: 2017年5月3日
 *      Author: clover
 */

#include "CLRwmf.h"
#include <fstream>

template <class T> void show_bit(const cv::Mat& mat){
	cv::Mat show = mat.clone();
	for(int y = 0; y< mat.rows; y++){
		for(int x = 0; x< mat.cols; x++){
			if(mat.ptr<T>(y)[x] > 0.5)
				show.ptr<uchar>(y)[x] = 255;
		}
	}
	cv::imshow("bit", show);
}

template <class T> void save_mat_to_file(const cv::Mat& src, std::string file_name){
//void save_mat_to_file(){
	std::ofstream outf;
	outf.open(file_name);

	for(int y = 0; y < src.rows; y++){
		const T* src_rows = src.ptr<T>(y);
		for(int x = 0; x < src.cols; x++){
			outf<< src_rows[x] << " ";
		}
		outf << std::endl;
	}
	outf.close();
}

void save_vector_to_file(std::vector<cv::Mat>& v, std::string file_name){

	std::ofstream outf;
	outf.open(file_name);

	for(int i = 0;i<v.size();i++){
		outf<< v[i];
		if(i%752 == 0)
			outf << std::endl;
	}

	outf.close();
}


template <class T> void repmat_row(const cv::Mat& imSrc, cv::Mat& imDst, int scalar_y, int rep_rows_num)
{
	imDst = cv::Mat(scalar_y, imSrc.cols, imSrc.type(), cv::Scalar::all(0));
	for(int y = 0; y < scalar_y; y++){
		const T* src_rows = imSrc.ptr<T>(rep_rows_num);
		T* dst_rows = imDst.ptr<T>(y);
		for(int x = 0; x < imDst.cols; x++){
			dst_rows[x] = src_rows[x];
		}
	}
}

template <class T> void repmat_col(const cv::Mat& imSrc, cv::Mat& imDst, int scalar_x, int rep_col_num)
{
	imDst = cv::Mat(imSrc.rows, scalar_x, imSrc.type(), cv::Scalar::all(0));
	for(int y = 0; y < imDst.rows; y++){
		const T* src_rows_data = imSrc.ptr<T>(y);
		T* dst_rows_data = imDst.ptr<T>(y);
		for(int x = 0; x < imDst.cols; x++){
			dst_rows_data[x] = src_rows_data[rep_col_num];
		}
	}
}

template <class T> void mat_row_cpy(const cv::Mat& imSrc, cv::Mat& imDsi,
		int from_start_row, int from_end_row, int to_start_row, int to_end_row){
	try {
		if(imDsi.empty())
			throw "imDsi.empty()";
		if(imSrc.type() != imDsi.type())
			throw "imSrc.type() != imDsi.type()";
		if(imSrc.cols != imDsi.cols)
			throw "imSrc.cols != imDsi.cols";
		//if(imSrc.rows != imDsi.rows)
		//	throw "imSrc.rows != imDsi.rows";
		if((from_end_row - from_start_row) != (to_end_row - to_start_row))
			throw "(from_end_row - from_start_row) != (to_end_row - to_start_row)";

		cv::Mat temp_mat = cv::Mat(from_end_row - from_start_row, imSrc.cols, imSrc.type(), cv::Scalar::all(0));

		for(int y = from_start_row; y < from_end_row; y++){
			const T* imSrc_row_data = imSrc.ptr<T>(y);
			T* temp_mat_row_data = temp_mat.ptr<T>(y - from_start_row);
			for(int x = 0; x < imSrc.cols; x ++){
				temp_mat_row_data[x] = imSrc_row_data[x];
			}
		}

		for(int y = to_start_row; y < to_end_row; y++){
			T* dsi_row_data = imDsi.ptr<T>(y);
			T* temp_mat_row_data = temp_mat.ptr<T>(y - to_start_row);
			for(int x = 0; x < imDsi.cols; x ++){
				dsi_row_data[x] = temp_mat_row_data[x];
			}
		}
	} catch (char* msg) {
		std::cout << msg <<std::endl;
	}
}

template <class T> void mat_col_cpy(const cv::Mat& imSrc, cv::Mat& imDsi,
		int from_start_col, int from_end_col, int to_start_col, int to_end_col){
	try {
		if(imSrc.rows != imDsi.rows)
			throw "imSrc.rows != imDsi.rows";
		if(imSrc.type() != imDsi.type())
			throw "imSrc.type() != imDsi.type()";
		if((from_end_col - from_start_col) != (to_end_col - to_start_col))
			throw "(from_end_col - from_start_col) != (to_end_col - to_start_col)";
		cv::Mat temp_mat = cv::Mat(imSrc.rows, from_end_col - from_start_col, imSrc.type(), cv::Scalar::all(0));
		for(int y = 0; y < imSrc.rows; y++){
			const T* src_row_data = imSrc.ptr<T>(y);
			T* temp_row_data = temp_mat.ptr<T>(y);
			for(int x = from_start_col; x < from_end_col; x++){
				temp_row_data[x - from_start_col] = src_row_data[x];
			}
		}

		for(int y = 0; y < imDsi.rows; y++){
			T* dsi_row_data = imDsi.ptr<T>(y);
			T* temp_row_data = temp_mat.ptr<T>(y);
			for(int x = to_start_col; x < to_end_col; x++){
				dsi_row_data[x] = temp_row_data[x - to_start_col];
			}
		}

	} catch (char* msg) {
		std::cout<< msg << std::endl;
	}
}

template <class T> cv::Mat cumsum(const cv::Mat& imSrc, int r)
{
	cv::Mat imDsi = imSrc.clone(); //cv::Mat(imSrc.size(), imSrc.type(), cv::Scalar::all(0));
	if(r == 1)
		for(int y = 1; y < imSrc.rows; y++){
			//const T* last_src_row_data = imSrc.ptr<T>(y - 1);
			//const T* current_src_row_data = imSrc.ptr<T>(y);
			T* last_dsi_row_data = imDsi.ptr<T>(y - 1);
			T* current_dsi_row_data = imDsi.ptr<T>(y);

			for(int x = 0; x < imSrc.cols; x++){
				current_dsi_row_data[x] = last_dsi_row_data[x] + current_dsi_row_data[x];
			}
		}
	else
		for(int y = 0; y < imSrc.rows; y++){
			T* current_dsi_row_data = imDsi.ptr<T>(y);
			for(int x = 1; x < imSrc.cols; x++){
				current_dsi_row_data[x] = current_dsi_row_data[x - 1] + current_dsi_row_data[x];
			}
		}
	return imDsi;
}

template <class T> cv::Mat box_filter(const cv::Mat& imSrc, int r){

	int hei = imSrc.rows;
	int wid = imSrc.cols;
	cv::Mat imDsi = cv::Mat(hei, wid, imSrc.type(), cv::Scalar::all(0));
	cv::Mat temp1 = cv::Mat(hei, wid, imSrc.type(), cv::Scalar::all(0));
	cv::Mat temp2 = cv::Mat(hei, wid, imSrc.type(), cv::Scalar::all(0));
	cv::Mat imCum;
	imCum =  cumsum<T>(imSrc, 1);

	mat_row_cpy<T>(imCum, imDsi, r, 2*r+1, 0, r+1); // imDst(1:r+1, :) = imCum(1+r:2*r+1, :) //c语言都是 0~47

	//imDst(r+2:hei-r, :) = imCum(2*r+2:hei, :) - imCum(1:hei-2*r-1, :);
	mat_row_cpy<T>(imCum, temp1, 2*r+1, hei, 0, hei - (2*r+1));
	mat_row_cpy<T>(imCum, temp2, 0, hei-2*r-1, 0, hei - (2*r+1));
	temp1 = temp1 - temp2;
	mat_row_cpy<T>(temp1, imDsi, 0, hei - 2*r - 1, r+1, hei - r);
	//end
	temp1 = cv::Scalar::all(0);
	temp2 = cv::Scalar::all(0);
	cv::Mat repMat;
	//imDst(hei-r+1:hei, :) = repmat(imCum(hei, :), [r, 1]) - imCum(hei-2*r:hei-r-1, :);
	repmat_row<T>(imCum,repMat,r,hei - 1);
	temp1 = cv::Mat(r, imCum.cols, imCum.type(), cv::Scalar::all(0) );
	mat_row_cpy<T>(imCum, temp1, hei - 2*r - 1, hei - r - 1, 0, r);
	temp1 = repMat - temp1;
	mat_row_cpy<T>(temp1, imDsi, 0, r, hei-r, hei);
	//end
	//------------------------
	imCum = cumsum<T>(imDsi, 2);
	//imDst(:, 1:r+1) = imCum(:, 1+r:2*r+1);
	mat_col_cpy<T>(imCum, imDsi, r, 2*r + 1, 0, r+1);
	//end
	//imDst(:, r+2:wid-r) = imCum(:, 2*r+2:wid) - imCum(:, 1:wid-2*r-1);
	temp1 = cv::Mat(imCum.rows, wid - (2*r + 2 - 1) , imCum.type(), cv::Scalar::all(0));
	mat_col_cpy<T>(imCum, temp1, 2*r + 1, wid, 0, wid - 2*r - 1);
	temp2 = cv::Mat(imCum.rows, wid - 2*r - 1, imCum.type(), cv::Scalar::all(0));
	mat_col_cpy<T>(imCum, temp2, 0, wid - 2*r - 1, 0, wid - 2*r - 1);
	temp1 = temp1 - temp2;
	mat_col_cpy<T>(temp1, imDsi, 0, wid - 2*r - 1, r + 1, wid - r);
	//end
	//imDst(:, wid-r+1:wid) = repmat(imCum(:, wid), [1, r]) - imCum(:, wid-2*r:wid-r-1);
	repmat_col<T>(imCum, temp1, r, wid - 1);
	temp2 = cv::Mat(imCum.rows, r, imCum.type(), cv::Scalar::all(0));
	mat_col_cpy<T>(imCum, temp2, wid - 2*r - 1, wid - r - 1  , 0, r);
	temp1 = temp1 - temp2;
	mat_col_cpy<T>(temp1, imDsi, 0, r, wid-r, wid);
	//end
#if 0
#endif
	return imDsi;
	/*
	[hei, wid] = size(imSrc);
	imDst = zeros(size(imSrc));

	%cumulative sum over Y axis
	imCum = cumsum(imSrc, 1);
	%difference over Y axis
	imDst(1:r+1, :) = imCum(1+r:2*r+1, :);
	imDst(r+2:hei-r, :) = imCum(2*r+2:hei, :) - imCum(1:hei-2*r-1, :);
	imDst(hei-r+1:hei, :) = repmat(imCum(hei, :), [r, 1]) - imCum(hei-2*r:hei-r-1, :);

	%cumulative sum over X axis
	imCum = cumsum(imDst, 2);
	%difference over Y axis
	imDst(:, 1:r+1) = imCum(:, 1+r:2*r+1);
	imDst(:, r+2:wid-r) = imCum(:, 2*r+2:wid) - imCum(:, 1:wid-2*r-1);
	imDst(:, wid-r+1:wid) = repmat(imCum(:, wid), [1, r]) - imCum(:, wid-2*r:wid-r-1);

*/
}

void split_bgr(const cv::Mat& bgr_im, cv::Mat& b_im , cv::Mat& g_im , cv::Mat& r_im ){
	std::vector<cv::Mat> mult_mat;
	cv::split(bgr_im, mult_mat );

	if(mult_mat.size() > 1){
		b_im = mult_mat[0].clone();
		g_im = mult_mat[1].clone();
		r_im = mult_mat[2].clone();
	}else{
		b_im = bgr_im.clone();
		g_im = bgr_im.clone();
		r_im = bgr_im.clone();
	}
}

template <class T> cv::Mat get_var_I(const cv::Mat& double_img_a, const cv::Mat& double_img_b,
		const cv::Mat& double_mean_img_a, const cv::Mat& double_mean_img_b,
		int r, const cv::Mat& N ) {
	cv::Mat mult_I;
	cv::Mat mult_mean_I;

	cv::multiply(double_img_a, double_img_b, mult_I);
	cv::multiply(double_mean_img_a, double_mean_img_b, mult_mean_I);
	return box_filter<double>(mult_I,r)/N - mult_mean_I;
}

template <class T> void fill_invSigma(CLRwmf_gfobj* a_gobj)
{
	a_gobj->m_invSigma.clear();
	cv::Mat sigma_mat;
	cv::Mat eye_mat = cv::Mat::eye(cv::Size(3,3), CV_64FC1);
	double sigma_array[9];
	int hei = a_gobj->m_image.rows;
	int wid = a_gobj->m_image.cols;
	for(int y = 0; y < hei; y++){
		T* var_rr_row = a_gobj->m_var_I_rr.ptr<T>(y);
		T* var_rg_row = a_gobj->m_var_I_rg.ptr<T>(y);
		T* var_rb_row = a_gobj->m_var_I_rb.ptr<T>(y);
		T* var_gg_row = a_gobj->m_var_I_gg.ptr<T>(y);
		T* var_gb_row = a_gobj->m_var_I_gb.ptr<T>(y);
		T* var_bb_row = a_gobj->m_var_I_bb.ptr<T>(y);
		for(int x = 0; x < wid; x++){
			sigma_array[0] = var_rr_row[x];
			sigma_array[1] = var_rg_row[x];
			sigma_array[2] = var_rb_row[x];
			sigma_array[3] = var_rg_row[x];
			sigma_array[4] = var_gg_row[x];
			sigma_array[5] = var_gb_row[x];
			sigma_array[6] = var_rb_row[x];
			sigma_array[7] = var_gb_row[x];
			sigma_array[8] = var_bb_row[x];
			sigma_mat = cv::Mat(3,3, CV_64FC1, sigma_array);
			cv::Mat temp;
			cv::invert(sigma_mat + a_gobj->m_eps*eye_mat, temp );
			a_gobj->m_invSigma.push_back(temp);
		}
	}
save_vector_to_file(a_gobj->m_invSigma, "vector.txt");
}


void guidedfilter_color_precompute(const cv::Mat& image, int r, double eps , CLRwmf_gfobj* a_gobj){
	a_gobj->m_image = image.clone();
	a_gobj->m_r = r;
	a_gobj->m_eps = eps;
	int hei = image.rows;
	int wid = image.cols;
	cv::Mat one_box = cv::Mat(hei, wid, CV_64FC1, cv::Scalar::all(1.0));
	a_gobj->m_N = box_filter<double>(one_box,r);
	cv::Mat img_r;
	cv::Mat img_g;
	cv::Mat img_b;
	split_bgr(image, img_b, img_g, img_r);

	cv::Mat box_r = box_filter<double>(img_r, r);
	cv::Mat box_g = box_filter<double>(img_g, r);
	cv::Mat box_b = box_filter<double>(img_b, r);
	cv::divide(box_r, a_gobj->m_N, a_gobj->m_mean_I_r);
	cv::divide(box_g, a_gobj->m_N, a_gobj->m_mean_I_g);
	cv::divide(box_b, a_gobj->m_N, a_gobj->m_mean_I_b);

	//	cv::multiply(img_r, img_r, mult_rr);
	//	cv::multiply(a_gobj->m_mean_I_r, a_gobj->m_mean_I_r, mult_mean_rr);
	//	box_filter<double>(mult_rr,r)/a_gobj->m_N - mult_mean_rr;
	a_gobj->m_var_I_rr =  get_var_I<double>(img_r, img_r,
			a_gobj->m_mean_I_r, a_gobj->m_mean_I_r,
			r,a_gobj->m_N);

	a_gobj->m_var_I_rg =  get_var_I<double>(img_r, img_g,
			a_gobj->m_mean_I_r, a_gobj->m_mean_I_g,
			r,a_gobj->m_N);

	a_gobj->m_var_I_rb =  get_var_I<double>(img_r, img_b,
			a_gobj->m_mean_I_r, a_gobj->m_mean_I_b,
			r,a_gobj->m_N);

	a_gobj->m_var_I_gg =  get_var_I<double>(img_g, img_g,
			a_gobj->m_mean_I_g, a_gobj->m_mean_I_g,
			r,a_gobj->m_N);

	a_gobj->m_var_I_gb =  get_var_I<double>(img_g, img_b,
			a_gobj->m_mean_I_g, a_gobj->m_mean_I_b,
			r,a_gobj->m_N);

	a_gobj->m_var_I_bb =  get_var_I<double>(img_b, img_b,
			a_gobj->m_mean_I_b, a_gobj->m_mean_I_b,
			r,a_gobj->m_N);
	//save_mat_to_file<double>(a_gobj->m_var_I_bb,"var_bb.txt");
	fill_invSigma<double>(a_gobj);

//	save_mat_to_file<double>(a_gobj->m_var_I_rr,"var_rr.txt");


/*
	cv::multiply(img_r, img_g, mult_rg);
	cv::multiply(a_gobj->m_mean_I_r, a_gobj->m_mean_I_g, mult_mean_rg);
	a_gobj->m_var_I_rr = box_filter<double>(mult_rg,r)/a_gobj->m_N - mult_mean_rg;
*/
}

template <class T> cv::Mat get_mean_Ip(const cv::Mat& color_I, const cv::Mat& p, int r, const cv::Mat& N){
	//cv::Mat mean_Ip;
	cv::Mat mult_mat;
	cv::Mat div_mat;
	cv::multiply(color_I, p, mult_mat);
	div_mat = box_filter<T>(mult_mat, r);
	cv::divide(div_mat, N, div_mat);
	return div_mat;
}

//输入数据为disp分层的二指图，每次只有一个深度值输入
template <class T> cv::Mat guidedfilter_color_runfilter(const cv::Mat& disp_bit, CLRwmf_gfobj* a_gfobj){
	int hei = disp_bit.rows;
	int wid = disp_bit.cols;
	int r = a_gfobj->m_r;
	cv::Mat mean_p = box_filter<T>(disp_bit, r);
	cv::divide(mean_p, a_gfobj->m_N, mean_p);
	cv::Mat img_r;
	cv::Mat img_g;
	cv::Mat img_b;
	split_bgr(a_gfobj->m_image, img_b, img_g, img_r);
	cv::Mat mean_Ip_r = get_mean_Ip<T>(img_r, disp_bit, r, a_gfobj->m_N);
	cv::Mat mean_Ip_g = get_mean_Ip<T>(img_g, disp_bit, r, a_gfobj->m_N);
	cv::Mat mean_Ip_b = get_mean_Ip<T>(img_b, disp_bit, r, a_gfobj->m_N);
	cv::Mat mult_mean_Ir_p;
	cv::Mat mult_mean_Ig_p;
	cv::Mat mult_mean_Ib_p;
	cv::multiply(a_gfobj->m_mean_I_r, mean_p, mult_mean_Ir_p);
	cv::multiply(a_gfobj->m_mean_I_g, mean_p, mult_mean_Ig_p);
	cv::multiply(a_gfobj->m_mean_I_b, mean_p, mult_mean_Ib_p);
	cv::Mat cov_Ip_r = mean_Ip_r - mult_mean_Ir_p;
	cv::Mat cov_Ip_g = mean_Ip_g - mult_mean_Ig_p;
	cv::Mat cov_Ip_b = mean_Ip_b - mult_mean_Ib_p;

	cv::Mat sigma;
	int num = 0;

	cv::Mat a_r = cv::Mat(hei, wid, CV_64FC1, cv::Scalar::all(0));
	cv::Mat a_g = cv::Mat(hei, wid, CV_64FC1, cv::Scalar::all(0));
	cv::Mat a_b = cv::Mat(hei, wid, CV_64FC1, cv::Scalar::all(0));

	double temp_array[3] = {0};
	cv::Mat cov_Ip;

	for(int y = 0; y < hei; y++){
		T* cov_Ip_r_row = cov_Ip_r.ptr<T>(y);
		T* cov_Ip_g_row = cov_Ip_g.ptr<T>(y);
		T* cov_Ip_b_row = cov_Ip_b.ptr<T>(y);
		T* a_r_row = a_r.ptr<T>(y);
		T* a_g_row = a_g.ptr<T>(y);
		T* a_b_row = a_b.ptr<T>(y);
		for(int x = 0; x < wid; x++){
			num = (y*wid) + x;
			sigma = a_gfobj->m_invSigma[num];
			temp_array[0] = cov_Ip_r_row[x];
			temp_array[1] = cov_Ip_g_row[x];
			temp_array[2] = cov_Ip_b_row[x];
			cov_Ip = cv::Mat(1,3,CV_64FC1,temp_array);
			cov_Ip = cov_Ip * sigma;
			a_r_row[x] = cov_Ip.ptr<T>(0)[0];
			a_g_row[x] = cov_Ip.ptr<T>(0)[1];
			a_b_row[x] = cov_Ip.ptr<T>(0)[2];
			//temp_mat.data<double>(0);
		}
	}

	cv::Mat b_mat = mean_p - a_r.mul(a_gfobj->m_mean_I_r) - a_g.mul(a_gfobj->m_mean_I_g) - a_b.mul(a_gfobj->m_mean_I_b);
	cv::Mat q = box_filter<double>(a_r, r).mul(img_r) +
				box_filter<double>(a_g, r).mul(img_g) +
				box_filter<double>(a_b, r).mul(img_b) +
				box_filter<double>(b_mat, r).mul(a_gfobj->m_N);
	return q;
}

void slice_disp(const cv::Mat& disp, uchar find_disp, cv::Mat& slice_bit_mat, uchar scalar)
{
	int hei = disp.rows;
	int wid = disp.cols;
	slice_bit_mat = cv::Mat(hei, wid, CV_64FC1, cv::Scalar::all(0));
	for(int y = 0; y < hei; y++){
		const uchar* disp_row = disp.ptr<uchar>(y);
		double* bit_row = slice_bit_mat.ptr<double>(y);
		for(int x = 0; x < wid; x++){
			if(disp_row[x]/scalar == find_disp){
				bit_row[x] = 1.0;
			}
		}
	}
}

void fill_filter_out(cv::Mat& disp_out, cv::Mat& ref_mat, uchar disp_will_fill, uchar scalar){
	int hei = disp_out.rows;
	int wid = disp_out.cols;
	for(int y = 0; y < hei; y++){
		uchar* disp_rows = disp_out.ptr<uchar>(y);
		double* ref_mat_rows = ref_mat.ptr<double>(y);
		for(int x = 0; x < wid; x++){
			uchar disp = disp_rows[x];
			double ref = ref_mat_rows[x];
			if(ref > 0.5 && disp == 0){
				disp_rows[x] = disp_will_fill*scalar;
			}
		}
	}

}

CLRwmf::CLRwmf() {
	// TODO Auto-generated constructor stub
	m_gfobj = new CLRwmf_gfobj();
}

CLRwmf::~CLRwmf() {
	// TODO Auto-generated destructor stub
	delete m_gfobj;
}

void CLRwmf::run_wmf(const cv::Mat& src_mat,const cv::Mat& disp)
{
	//cv::Mat box = box_filter<uchar>(src_mat, 4);
	cv::Mat src_double;
	src_mat.convertTo(src_double, CV_64FC3, 1.0/255.0);

	int disp_num = 16;
	uchar scalar = 16;

	guidedfilter_color_precompute(src_double, 10, 1.0e-4,m_gfobj);
	cv::Mat imgOL;
	cv::Mat imgAccum = cv::Mat(src_mat.rows, src_mat.cols, CV_64FC1, cv::Scalar::all(0));

	cv::Mat idxSelected;
	cv::Mat disp_out = cv::Mat(src_mat.size(), CV_8UC1, cv::Scalar::all(0));


	for(int i = 1; i < disp_num + 1; i++ )
	{
		cv::Mat slice_mat;
		slice_disp(disp, i, slice_mat,scalar);
		//show_bit(slice_mat);
		imgOL = guidedfilter_color_runfilter<double>(slice_mat, m_gfobj);
		imgAccum = imgAccum + imgOL;
		fill_filter_out(disp_out, imgAccum, i, scalar);
		//show_bit<uchar>(disp_out );
		//cv::waitKey(0);
	}
	disp_out = disp_out*16;
	cv::imshow("dd", disp_out);
	//std::cout << " src: " <<std::endl;
	//std::cout << src_mat <<std::endl;
	//std::cout << " box: "  << std::endl;
	//std::cout << box << std::endl;

	//cv::imshow("div", c);
	cv::waitKey(0);
}
