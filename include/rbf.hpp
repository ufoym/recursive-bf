#ifndef INCLUDE_RBF
#define INCLUDE_RBF
#include <math.h>
#include <string.h>
#define QX_DEF_CHAR_MAX 255

/* ======================================================================

RecursiveBF: A lightweight library for recursive bilateral filtering.

-------------------------------------------------------------------------

Intro:      Recursive bilateral filtering (developed by Qingxiong Yang) 
            is pretty fast compared with most edge-preserving filtering 
            methods.

            -   computational complexity is linear in both input size and 
                dimensionality
            -   takes about 43 ms to process a one mega-pixel color image
                (i7 1.8GHz & 4GB memory)
            -   about 18x faster than Fast high-dimensional filtering 
                using the permutohedral lattice
            -   about 86x faster than Gaussian kd-trees for fast high-
                dimensional filtering


Usage:      // ----------------------------------------------------------
            // Basic Usage
            // ----------------------------------------------------------

            unsigned char * img = ...;                    // input image
            unsigned char * img_out = 0;            // output image
            int width = ..., height = ..., channel = ...; // image size
            recursive_bf(img, img_out, 
                         sigma_spatial, sigma_range, 
                         width, height, channel);

            // ----------------------------------------------------------
            // Advanced: using external buffer for better performance
            // ----------------------------------------------------------

            unsigned char * img = ...;                    // input image
            unsigned char * img_out = 0;            // output image
            int width = ..., height = ..., channel = ...; // image size
            float * buffer = new float[                   // external buf
                                 ( width * height* channel 
                                 + width * height
                                 + width * channel 
                                 + width) * 2];
            recursive_bf(img, img_out, 
                         sigma_spatial, sigma_range, 
                         width, height, channel, 
                         buffer);
            delete[] buffer;


Notice:     Large sigma_spatial/sigma_range parameter may results in 
            visible artifact which can be removed by an additional 
            filter with small sigma_spatial/sigma_range parameter.

-------------------------------------------------------------------------

Reference:  Qingxiong Yang, Recursive Bilateral Filtering,
            European Conference on Computer Vision (ECCV) 2012, 399-413.

====================================================================== */

inline void recursive_bf(
    unsigned char * img_in, 
    unsigned char *& img_out, 
    float sigma_spatial, float sigma_range, 
    int width, int height, int channel, 
    float * buffer /*= 0*/);

// ----------------------------------------------------------------------

inline void _recursive_bf(
    unsigned char * img,
    float sigma_spatial, float sigma_range, 
    int width, int height, int channel,
    float * buffer = 0)
{
    const int width_height = width * height;
    const int width_channel = width * channel;
    const int width_height_channel = width * height * channel;

    bool is_buffer_internal = (buffer == 0);
    if (is_buffer_internal)
        buffer = new float[(width_height_channel + width_height 
                            + width_channel + width) * 2];

    float * img_out_f = buffer;
    float * img_temp = &img_out_f[width_height_channel];
    float * map_factor_a = &img_temp[width_height_channel];
    float * map_factor_b = &map_factor_a[width_height]; 
    float * slice_factor_a = &map_factor_b[width_height];
    float * slice_factor_b = &slice_factor_a[width_channel];
    float * line_factor_a = &slice_factor_b[width_channel];
    float * line_factor_b = &line_factor_a[width];
    
    //compute a lookup table
    float range_table[QX_DEF_CHAR_MAX + 1];
    float inv_sigma_range = 1.0f / (sigma_range * QX_DEF_CHAR_MAX);
    for (int i = 0; i <= QX_DEF_CHAR_MAX; i++) 
        range_table[i] = static_cast<float>(exp(-i * inv_sigma_range));

    float alpha = static_cast<float>(exp(-sqrt(2.0) / (sigma_spatial * width)));
    float ypr, ypg, ypb, ycr, ycg, ycb;
    float fp, fc;
    float inv_alpha_ = 1 - alpha;
    for (int y = 0; y < height; y++)
    {
        float * temp_x = &img_temp[y * width_channel];
        unsigned char * in_x = &img[y * width_channel];
        unsigned char * texture_x = &img[y * width_channel];
        *temp_x++ = ypr = *in_x++; 
        *temp_x++ = ypg = *in_x++; 
        *temp_x++ = ypb = *in_x++;
        unsigned char tpr = *texture_x++; 
        unsigned char tpg = *texture_x++;
        unsigned char tpb = *texture_x++;

        float * temp_factor_x = &map_factor_a[y * width];
        *temp_factor_x++ = fp = 1;

        // from left to right
        for (int x = 1; x < width; x++) 
        {
            unsigned char tcr = *texture_x++; 
            unsigned char tcg = *texture_x++; 
            unsigned char tcb = *texture_x++;
            unsigned char dr = abs(tcr - tpr);
            unsigned char dg = abs(tcg - tpg);
            unsigned char db = abs(tcb - tpb);
            int range_dist = (((dr << 1) + dg + db) >> 2);
            float weight = range_table[range_dist];
            float alpha_ = weight*alpha;
            *temp_x++ = ycr = inv_alpha_*(*in_x++) + alpha_*ypr; 
            *temp_x++ = ycg = inv_alpha_*(*in_x++) + alpha_*ypg; 
            *temp_x++ = ycb = inv_alpha_*(*in_x++) + alpha_*ypb;
            tpr = tcr; tpg = tcg; tpb = tcb;
            ypr = ycr; ypg = ycg; ypb = ycb;
            *temp_factor_x++ = fc = inv_alpha_ + alpha_*fp;
            fp = fc;
        }
        *--temp_x; *temp_x = 0.5f*((*temp_x) + (*--in_x));
        *--temp_x; *temp_x = 0.5f*((*temp_x) + (*--in_x));
        *--temp_x; *temp_x = 0.5f*((*temp_x) + (*--in_x));
        tpr = *--texture_x; 
        tpg = *--texture_x; 
        tpb = *--texture_x;
        ypr = *in_x; ypg = *in_x; ypb = *in_x;

        *--temp_factor_x; *temp_factor_x = 0.5f*((*temp_factor_x) + 1);
        fp = 1;

        // from right to left
        for (int x = width - 2; x >= 0; x--) 
        {
            unsigned char tcr = *--texture_x; 
            unsigned char tcg = *--texture_x; 
            unsigned char tcb = *--texture_x;
            unsigned char dr = abs(tcr - tpr);
            unsigned char dg = abs(tcg - tpg);
            unsigned char db = abs(tcb - tpb);
            int range_dist = (((dr << 1) + dg + db) >> 2);
            float weight = range_table[range_dist];
            float alpha_ = weight * alpha;

            ycr = inv_alpha_ * (*--in_x) + alpha_ * ypr; 
            ycg = inv_alpha_ * (*--in_x) + alpha_ * ypg; 
            ycb = inv_alpha_ * (*--in_x) + alpha_ * ypb;
            *--temp_x; *temp_x = 0.5f*((*temp_x) + ycr);
            *--temp_x; *temp_x = 0.5f*((*temp_x) + ycg);
            *--temp_x; *temp_x = 0.5f*((*temp_x) + ycb);
            tpr = tcr; tpg = tcg; tpb = tcb;
            ypr = ycr; ypg = ycg; ypb = ycb;

            fc = inv_alpha_ + alpha_*fp;
            *--temp_factor_x; 
            *temp_factor_x = 0.5f*((*temp_factor_x) + fc);
            fp = fc;
        }
    }
    alpha = static_cast<float>(exp(-sqrt(2.0) / (sigma_spatial * height)));
    inv_alpha_ = 1 - alpha;
    float * ycy, * ypy, * xcy;
    unsigned char * tcy, * tpy;
    memcpy(img_out_f, img_temp, sizeof(float)* width_channel);

    float * in_factor = map_factor_a;
    float*ycf, *ypf, *xcf;
    memcpy(map_factor_b, in_factor, sizeof(float) * width);
    for (int y = 1; y < height; y++)
    {
        tpy = &img[(y - 1) * width_channel];
        tcy = &img[y * width_channel];
        xcy = &img_temp[y * width_channel];
        ypy = &img_out_f[(y - 1) * width_channel];
        ycy = &img_out_f[y * width_channel];

        xcf = &in_factor[y * width];
        ypf = &map_factor_b[(y - 1) * width];
        ycf = &map_factor_b[y * width];
        for (int x = 0; x < width; x++)
        {
            unsigned char dr = abs((*tcy++) - (*tpy++));
            unsigned char dg = abs((*tcy++) - (*tpy++));
            unsigned char db = abs((*tcy++) - (*tpy++));
            int range_dist = (((dr << 1) + dg + db) >> 2);
            float weight = range_table[range_dist];
            float alpha_ = weight*alpha;
            for (int c = 0; c < channel; c++) 
                *ycy++ = inv_alpha_*(*xcy++) + alpha_*(*ypy++);
            *ycf++ = inv_alpha_*(*xcf++) + alpha_*(*ypf++);
        }
    }
    int h1 = height - 1;
    ycf = line_factor_a;
    ypf = line_factor_b;
    memcpy(ypf, &in_factor[h1 * width], sizeof(float) * width);
    for (int x = 0; x < width; x++) 
        map_factor_b[h1 * width + x] = 0.5f*(map_factor_b[h1 * width + x] + ypf[x]);

    ycy = slice_factor_a;
    ypy = slice_factor_b;
    memcpy(ypy, &img_temp[h1 * width_channel], sizeof(float)* width_channel);
    int k = 0; 
    for (int x = 0; x < width; x++) {
        for (int c = 0; c < channel; c++) {
            int idx = (h1 * width + x) * channel + c;
            img_out_f[idx] = 0.5f*(img_out_f[idx] + ypy[k++]) / map_factor_b[h1 * width + x];
        }
    }

    for (int y = h1 - 1; y >= 0; y--)
    {
        tpy = &img[(y + 1) * width_channel];
        tcy = &img[y * width_channel];
        xcy = &img_temp[y * width_channel];
        float*ycy_ = ycy;
        float*ypy_ = ypy;
        float*out_ = &img_out_f[y * width_channel];

        xcf = &in_factor[y * width];
        float*ycf_ = ycf;
        float*ypf_ = ypf;
        float*factor_ = &map_factor_b[y * width];
        for (int x = 0; x < width; x++)
        {
            unsigned char dr = abs((*tcy++) - (*tpy++));
            unsigned char dg = abs((*tcy++) - (*tpy++));
            unsigned char db = abs((*tcy++) - (*tpy++));
            int range_dist = (((dr << 1) + dg + db) >> 2);
            float weight = range_table[range_dist];
            float alpha_ = weight*alpha;

            float fcc = inv_alpha_*(*xcf++) + alpha_*(*ypf_++);
            *ycf_++ = fcc;
            *factor_ = 0.5f * (*factor_ + fcc);

            for (int c = 0; c < channel; c++)
            {
                float ycc = inv_alpha_*(*xcy++) + alpha_*(*ypy_++);
                *ycy_++ = ycc;
                *out_ = 0.5f * (*out_ + ycc) / (*factor_);
                *out_++;
            }
            *factor_++;
        }
        memcpy(ypy, ycy, sizeof(float) * width_channel);
        memcpy(ypf, ycf, sizeof(float) * width);
    }

    for (int i = 0; i < width_height_channel; ++i)
        img[i] = static_cast<unsigned char>(img_out_f[i]);

    if (is_buffer_internal)
        delete[] buffer;
}


inline void recursive_bf(
    unsigned char * img_in,
    unsigned char *& img_out,
    float sigma_spatial, float sigma_range,
    int width, int height, int channel,
    float * buffer = 0)
{
    if (img_out == 0)
        img_out = new unsigned char[width * height * channel];
    for (int i = 0; i < width * height * channel; ++i)
        img_out[i] = img_in[i];
    _recursive_bf(img_out, sigma_spatial, sigma_range, width, height, channel, buffer);
}

#endif // INCLUDE_RBF
