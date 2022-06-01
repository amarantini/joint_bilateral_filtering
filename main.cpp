#include "pfm.h"
#include "Vector.hpp"

#include <iostream>
#include <math.h>
#include <cstring>
#include <vector>
#include <algorithm>
#include <numeric>

using namespace std;

// indexing 1d image array
int imIdx(int i, int j, int width, int channels = 3) {
	return channels * width * i + channels * j;
}

void MeanSTD(std::vector<float> samples, float& mean,float& std)
{
	int size = samples.size();
	mean = std::accumulate(samples.begin(), samples.end(), 0.0)/size;

    float sqDiff = 0;
    for (int i = 0; i < size; i++)
        sqDiff += (samples[i] - mean) *
                  (samples[i] - mean);
	std = sqrt(sqDiff / size);
    return;
}


float* OutlierFilter(float* image, unsigned width, unsigned height) {
	// C_output
	float* img_out = new float[width*height*3];
	memcpy(img_out, image, sizeof(float)*width*height*3);
	float sigma_p = 5;
	float k = 5;
	for(int i=0; i<width; i++){
		for(int j=0; j<height; j++){
			vector<float> samples_r;
			vector<float> samples_b;
			vector<float> samples_g;
			for(int m=i-2*sigma_p; m<=i+2*sigma_p; m++){
				for(int n=j-2*sigma_p; n<=j+2*sigma_p; n++){
					if(m>=0 && m<width && n>=0 && n<height){
						samples_r.push_back(image[imIdx(n,m,width)]);
						samples_b.push_back(image[imIdx(n,m,width)+1]);
						samples_g.push_back(image[imIdx(n,m,width)+2]);
					} 
				}
			}
			float mean_r,std_r,mean_g,std_g,mean_b,std_b;
			MeanSTD(samples_r, mean_r,std_r);
			MeanSTD(samples_g, mean_g,std_g);
			MeanSTD(samples_b, mean_b,std_b);
			img_out[imIdx(j,i,width)] = clamp(image[imIdx(j,i,width)],mean_r-k*std_r, mean_r+k*std_r);
			img_out[imIdx(j,i,width)+1] = clamp(image[imIdx(j,i,width)+1],mean_g-k*std_g, mean_g+k*std_g);
			img_out[imIdx(j,i,width)+2] = clamp(image[imIdx(j,i,width)+2],mean_b-k*std_b, mean_b+k*std_b);
		}
	}
	return img_out;
}

float* GaussianFilter(float* image, unsigned width, unsigned height) {
	// C_output
	float* gauss = new float[width*height*3];
	memcpy(gauss, image, sizeof(float)*width*height*3);
	// Parameters
	float sigma_p = 5;
	// int window_size = 4 * sigma_p + 1;

	// TODO 1
	// Implement Gaussian Filtering	
	for(int i=0; i<width; i++){
		for(int j=0; j<height; j++){
			for(int k=0; k<3; k++) {
				float sum_gauss_weight = 0;
				float sum_pixel_c = 0;
				for(int m=i-2*sigma_p; m<=i+2*sigma_p; m++){
					for(int n=j-2*sigma_p; n<=j+2*sigma_p; n++){
						if(m>=0 && m<width && n>=0 && n<height){
							float G_ij = exp(-1*((i-m)*(i-m)+(j-n)*(j-n))/2/pow(sigma_p,2));
							sum_gauss_weight += G_ij;
							sum_pixel_c += G_ij*image[imIdx(n,m,width)+k];
						} 
					}
				}
				gauss[imIdx(j,i,width)+k] = sum_pixel_c / sum_gauss_weight;
			}
		}
	}

	cout << "Gaussian Filtering -- DONE\n";
	return gauss;
}

float* BilateralFilter(float* image, unsigned width, unsigned height) {
	// C_output
	float* bilateral = new float[width*height*3];
	memcpy(bilateral, image, sizeof(float)*width*height*3);

	// Parameters
	float sigma_p = 5;
	float sigma_c = 0.1;

	// TODO 2
	// Implement Bilateral Filtering 
	for(int i=0; i<width; i++){
		for(int j=0; j<height; j++){
			float sum_bilateral_weight = 0;
			float sum_pixel_c1 = 0;
			float sum_pixel_c2 = 0;
			float sum_pixel_c3 = 0;
			Vector3f C1 = Vector3f(image[imIdx(j,i,width)],image[imIdx(j,i,width)+1],image[imIdx(j,i,width)+1]);
			for(int m=i-2*sigma_p; m<=i+2*sigma_p; m++){
				for(int n=j-2*sigma_p; n<=j+2*sigma_p; n++){
					if(m>=0 && m<width && n>=0 && n<height){
						Vector3f C2 = Vector3f(image[imIdx(n,m,width)],image[imIdx(n,m,width)+1],image[imIdx(n,m,width)+2]);
						float B_ij = exp(-1*((i-m)*(i-m)+(j-n)*(j-n))/2/pow(sigma_p,2) - pow((C1-C2).norm(),2)/2/pow(sigma_c,2));
						sum_bilateral_weight += B_ij;
						sum_pixel_c1 += B_ij*image[imIdx(n,m,width)];
						sum_pixel_c2 += B_ij*image[imIdx(n,m,width)+1];
						sum_pixel_c3 += B_ij*image[imIdx(n,m,width)+2];
					}
				}
			}
			if(sum_bilateral_weight!=0){
				bilateral[imIdx(j,i,width)] = sum_pixel_c1 /  (sum_bilateral_weight);
				bilateral[imIdx(j,i,width)+1] = sum_pixel_c2 /  (sum_bilateral_weight);
				bilateral[imIdx(j,i,width)+2] = sum_pixel_c3 /  (sum_bilateral_weight);
			}
		}
	}
	
	cout << "Bilateral Filtering -- DONE\n";
	return bilateral;
}

float* JointBilateralFilter(float* image, float* normal, float* position, unsigned width, unsigned height) {
	// C_output
	float* joint = new float[width*height*3];
	memcpy(joint, image, sizeof(float)*width*height*3);

	// Parameters
	float sigma_p = 5;
	float sigma_c = 0.3;
	float sigma_n = 0.1;
	float sigma_d = 0.1;

	// TODO 3
	// Implement Joint Bilateral Filtering
	for(int i=0; i<width; i++){
		for(int j=0; j<height; j++){
			float sum_joint_weight = 0;
			float sum_pixel_c1 = 0;
			float sum_pixel_c2 = 0;
			float sum_pixel_c3 = 0;
			Vector3f C1 = Vector3f(image[imIdx(j,i,width)],image[imIdx(j,i,width)+1],image[imIdx(j,i,width)+2]);
			Vector3f N1 = Vector3f(normal[imIdx(j,i,width)],normal[imIdx(j,i,width)+1],normal[imIdx(j,i,width)+2]);
			Vector3f P1 = Vector3f(position[imIdx(j,i,width)],position[imIdx(j,i,width)+1],position[imIdx(j,i,width)+2]);
			for(int m=i-2*sigma_p; m<=i+2*sigma_p; m++){
				for(int n=j-2*sigma_p; n<=j+2*sigma_p; n++){
					if(m>=0 && m<width && n>=0 && n<height){
						Vector3f C2 = Vector3f(image[imIdx(n,m,width)],image[imIdx(n,m,width)+1],image[imIdx(n,m,width)+2]);
						Vector3f N2 = Vector3f(normal[imIdx(n,m,width)],normal[imIdx(n,m,width)+1],normal[imIdx(n,m,width)+2]);
						Vector3f P2 = Vector3f(position[imIdx(n,m,width)],position[imIdx(n,m,width)+1],position[imIdx(n,m,width)+2]);
						float D_normal = acos(std::clamp(dotProduct(N1,N2),(float)0.0,(float)1.0));
						Vector3f P21 = P2 - P1;
						float D_plane = dotProduct(N1,(P2-P1).normalized());
						if(P21.x==0 && P21.y==0 && P21.z==0){
							D_plane = 0;
						}
						float J_ij = exp(-1*((i-m)*(i-m)+(j-n)*(j-n))/2/pow(sigma_p,2) 
										- pow((C1-C2).norm(),2)/2/pow(sigma_c,2)
										- pow(D_normal,2)/2/pow(sigma_n,2)
										- pow(D_plane,2)/2/pow(sigma_d,2));
						sum_joint_weight += J_ij;
						sum_pixel_c1 += J_ij*image[imIdx(n,m,width)];
						sum_pixel_c2 += J_ij*image[imIdx(n,m,width)+1];
						sum_pixel_c3 += J_ij*image[imIdx(n,m,width)+2];
					}
				}
			}
			// if(sum_joint_weight!=0) {
				joint[imIdx(j,i,width)] = sum_pixel_c1 / (sum_joint_weight+0.000001);
				joint[imIdx(j,i,width)+1] = sum_pixel_c2 / (sum_joint_weight+0.000001);
				joint[imIdx(j,i,width)+2] = sum_pixel_c3 / (sum_joint_weight+0.000001);
				
			// }	
		}
	}

	cout << "Joint Bilateral Filtering -- DONE\n";
	return joint;
}



int main() {
    unsigned *w = new unsigned;
    unsigned *h = new unsigned;

	// Input buffers
	float* imageBuffer = read_pfm_file3("im/s2.pfm", w, h);      		 // load image buffer - 3 channels
    float* normalBuffer = read_pfm_file3("im/s2_normal.pfm", w, h);      // load normal buffer - 3 channels
	float* positionBuffer = read_pfm_file3("im/s2_position.pfm", w, h);  // load position buffer - 3 channels
	float* imageBuffer_filtered = OutlierFilter(imageBuffer, *w, *h);  //filter outliers
	
	float* gausssianFiltered = GaussianFilter(imageBuffer, *w, *h);
	write_pfm_file3("im/s2_gaussian.pfm", gausssianFiltered, *w, *h);
	delete[] gausssianFiltered;
	
	float* bilateralFiltered = BilateralFilter(imageBuffer_filtered, *w, *h);
	write_pfm_file3("im/s2_bilateral_f.pfm", bilateralFiltered, *w, *h);
	delete[] bilateralFiltered;
	
	float* jointBilateralFiltered = JointBilateralFilter(imageBuffer_filtered, normalBuffer, positionBuffer, *w, *h);
	write_pfm_file3("im/s2_jointBilateral_f1.pfm", jointBilateralFiltered, *w, *h);
	delete[] jointBilateralFiltered;


	// 
	delete[] imageBuffer_filtered;
	delete[] imageBuffer;
	delete[] normalBuffer;
	delete[] positionBuffer;
	delete w;
	delete h;
	
	return 0;
}