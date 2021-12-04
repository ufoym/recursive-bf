#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image.h"
#include "stb/stb_image_write.h"
#include "../include/rbf.hpp"
#include <stdio.h>
#include <time.h>


class Timer {
private:
	unsigned long begTime;
public:
	void start() { begTime = clock(); }
	float elapsedTime() { return float((unsigned long)clock() - begTime) / CLOCKS_PER_SEC; }
};

int main(int argc, char*argv[])
{
	if (argc != 5)
	{
		printf("Usage:\n");
		printf("--------------------------------------------------------------------\n\n");
		printf("rbf filename_out filename_in (only support ppm images) \n");
		printf("    sigma_spatial(e.g., 0.03) sigma_range(e.g., 0.1)\n\n");
		printf("--------------------------------------------------------------------\n");
		return(-1);
	}
	else
	{
		const int n = 100;
		const char * filename_out = argv[1];
		const char * filename_in = argv[2];
		float sigma_spatial = static_cast<float>(atof(argv[3]));
		float sigma_range = static_cast<float>(atof(argv[4]));

		int width, height, channel;
		unsigned char * img = stbi_load(filename_in, &width, &height, &channel, 0);
		unsigned char * img_out = 0;
		Timer timer;

		timer.start();
		for (int i = 0; i < n; ++i)
			recursive_bf(img, img_out, sigma_spatial, sigma_range, width, height, channel);
		printf("Internal Buffer: %2.5fsecs\n", timer.elapsedTime() / n);


		float * buffer = new float[(width * height* channel + width * height
									+ width * channel + width) * 2];
		timer.start();
		for (int i = 0; i < n; ++i)
			recursive_bf(img, img_out, sigma_spatial, sigma_range, width, height, channel, buffer);
		printf("External Buffer: %2.5fsecs\n", timer.elapsedTime() / n);
		delete[] buffer;

		stbi_write_bmp(filename_out, width, height, channel, img_out);
		delete[] img;
	}
}
