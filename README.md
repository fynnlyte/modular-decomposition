```
$ cd PATH/TO/MD/build
$ cmake ..
$ make  
```

or (for four kernels) `$ make -j 4`  

execute MD:

```
$ ./mod_dec -/graphname.file	> out_tmp.dot
$ xdot out_tmp.dot
```

some exmpl-graphs are in ./MD/build/graphs
out_tmp is in graphViz format. 
to modify it into png with:

```
$ dot -Tpng out_tmp -o out_tmp.png
$ eog out_tmp.png
```





