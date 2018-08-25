# modular-decomposition

This Code is part of Adrian Fritz's Master Thesis: "A Heuristic for Cograph-Editing", Saarland University, 2015. The algorithm itselft is based on the following paper:

Christian Capelle, Michel Habib, Fabien Montgolfier. Graph Decompositions and Factorizing
Permutations. Discrete Mathematics and Theoretical Computer Science, DMTCS, 2002, 5,
pp.55-70.


Copyright (C) 2015 Adrian Fritz, 2018 Fynn Leitow (only done some minor bugfixes)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.


### Compilation

Ensure that boost::graph v.1.55 or higher is installed:

$ cd PATH/TO/MD/build
$ cmake ..
$ make  
or (for four cores) $ make -j 4  


### Usage

$ ./mod_dec -your_graph_as_string > out.dot

convert to PNG via:

$ dot -Tpng out.dot -o out.png

Because I used this inside a Java-Project, the format of the graph is the same as a .toString() - call on a Graph in the JGraphT - Library. An example can be found at the bottom of the modDecopm.cpp - file. Feel free to use it directly on a boost::graph or implement your own parser.
