Recursive bilateral filtering (developed by Qingxiong Yang) is pretty fast compared with most edge-preserving filtering methods
- computational complexity is linear in both input size and dimensionality:
- takes about 43 ms to process a one megapixel color image (i7 1.8GHz & 4GB mem)
- about 18x faster than *Fast high-dimensional filtering using the permutohedral lattice*
- about 86x faster than *Gaussian kd-trees for fast high-dimensional filtering*

## Results
<table>
<tr>
<td><img src="https://cloud.githubusercontent.com/assets/2270240/26041579/7d7c034e-3960-11e7-9549-912685043e39.jpg" width="300px"><br/><p align="center">Original Image</p></td>
<td><img src="https://cloud.githubusercontent.com/assets/2270240/26041586/8b4afb42-3960-11e7-9bd8-62bbb924f1e9.jpg" width="300px"><br/><p align="center">OpenCV's BF (896ms)</p></td>
<td><img src="https://cloud.githubusercontent.com/assets/2270240/26041590/8d08c16c-3960-11e7-8a0c-95a77d6d9085.jpg" width="300px"><br/><p align="center">RecursiveBF (18ms)</p></td>
</tr>
<tr>
<td></td>
<td><img src="https://cloud.githubusercontent.com/assets/2270240/26041583/86ea7b22-3960-11e7-8ded-5109b76966ca.jpg" width="300px"><br/><p align="center">Gaussian Blur</p></td>
<td><img src="https://cloud.githubusercontent.com/assets/2270240/26041584/88dfc9b4-3960-11e7-8c9d-2634eac098d0.jpg" width="300px"><br/><p align="center">Median Blur</p></td>
</tr></table>

For more details of the algorithm, please refer to the original paper

    @inproceedings{yang2012recursive,
        title={Recursive bilateral filtering},
        author={Yang, Qingxiong},
        booktitle={European Conference on Computer Vision},
        pages={399--413},
        year={2012},
        organization={Springer}
    }

Optionally, you can cite this repo

    @misc{Ming2017,
      author = {Ming Yang},
      title = {A lightweight C++ library for recursive bilateral filtering},
      year = {2017},
      publisher = {GitHub},
      journal = {GitHub repository},
      howpublished = {\url{https://github.com/ufoym/RecursiveBF}}
    }
