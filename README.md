# ScaleSCAN: Scalable Density-based Graph Clustering
The codes of ScaleSCAN developed by Tomokatsu Takahashi for DEXA'18.
If you use our software, **please cite the following paper**.


```
Hiroaki Shiokawa, Tomokatsu Takahashi, Hiroyuki Kitagawa,
"ScaleSCAN: Scalable Density-based Graph Clustering,"
In Proceedings of the 29th International Conference on Database and Expert Systems Applications (DEXA 2018), pp.18-34, Regensburg, Germany, September 2018.
```
```
@inproceedings{ShiokawaTK18,
  author    = {Hiroaki Shiokawa and Tomokatsu Takahashi and Hiroyuki Kitagawa},
  title     = {{ScaleSCAN: Scalable Density-Based Graph Clustering}},
  booktitle = {Proceedings of the 29th International Conference on Database and Expert Systems Applications (DEXA 2018)},
  pages     = {18--34},
  year      = {2018},
  month	    = {9},
  url       = {https://doi.org/10.1007/978-3-319-98809-2\_2},
  doi       = {10.1007/978-3-319-98809-2\_2},
}
```

### LICENSE
This software is released under the MIT license. Please read `LICENSE.txt` for details.

## How to Use
### Requirements
ScaleSCAN requires the following softwares.
* gcc Version 4.8.5 (or later)
* Intel AVX2
* OpemMP library

We have confirmed that our software works on the following environments.
* Intel(R) Xeon(R) E5-2690
* Intel Core i5-8210Y

### Build
1. Install *gcc*.
2. Run `make`.
3. If you can find `scalescan` and `csr_gen`, the build has been successfully completed.

### Usage
#### Input file
Input file must be formatted as a list of edges included in a given graph; each line represents a pair of node IDs that are connected via an edge. The nodes in the same line must be spanned by `TAB` or `space` like as follows:
``` sample_graph.txt
1	90
1	109
1	207
1	282
1	321
1	699
...
```

#### File conversion
`scalescan` reads the given graph by the CSR format, and this requires a file conversion process. 
To covert the input file (`sample_graph.txt`) into the CSR format (`sample_graph.bin`), you should run `csr_gen` like as follows:
``` csr_gen
$ ./csr_gen sample_graph.txt sample_graph.bin
```

#### Clustering
Finally, we can run the clustering by `scalescan` like as follows:
```
$ ./scalescan sample_graph.bin <T> <EPS> <MU>
```

where `<T>` is # of threads, `<EPS>` is the epsilon value (0<= `<EPS>` <= 1), and `<MU>` is the mu value (`<MU>` >= 2).

## Reference
* Hiroaki Shiokawa, Tomokatsu Takahashi, Hiroyuki Kitagawa,"ScaleSCAN: Scalable Density-based Graph Clustering," In Proceedings of the 29th International Conference on Database and Expert Systems Applications (DEXA 2018), pp.18-34, Regensburg, Germany, September 2018.
