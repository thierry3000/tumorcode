/**
This file is part of tumorcode project.
(http://www.uni-saarland.de/fak7/rieger/homepage/research/tumor/tumor.html)

Copyright (C) 2016  Michael Welter and Thierry Fredrich

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
//=======================================================================
// Copyright 2001 Jeremy G. Siek, Andrew Lumsdaine, Lie-Quan Lee, 
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
#include <boost/config.hpp>
#include <vector>
#include <list>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <iterator>
#include <iostream>

namespace boost
{
  struct edge_component_t
  {
    enum
    { num = 555 };
    typedef edge_property_tag kind;
  }
  edge_component;
}

int
main()
{
  using namespace boost;
  typedef adjacency_list < vecS, vecS, directedS,
    no_property, property < edge_component_t, std::size_t > >graph_t;
  typedef graph_traits < graph_t >::vertex_descriptor vertex_t;
  graph_t g(11);
  add_edge(0, 2, g);
  add_edge(1, 0, g);
  add_edge(1, 3, g);
  add_edge(1, 5, g);
  add_edge(1, 7, g);
  add_edge(2, 10, g);
  add_edge(3, 4, g);
  add_edge(3, 9, g);
  add_edge(4, 10, g);
  add_edge(5, 6, g);
  add_edge(6, 10, g);
  add_edge(7, 8, g);
  add_edge(8, 10, g);

  property_map < graph_t, edge_component_t >::type
    component = get(edge_component, g);

  std::size_t num_comps = biconnected_components(g, component);
  std::cerr << "Found " << num_comps << " biconnected components.\n";

  std::vector<vertex_t> art_points;
  articulation_points(g, std::back_inserter(art_points));
  std::cerr << "Found " << art_points.size() << " articulation points.\n";

  std::cout << "graph A {\n" << "  node[shape=\"circle\"]\n";

  for (std::size_t i = 0; i < art_points.size(); ++i) {
    std::cout << (char)(art_points[i] + 'A') 
              << " [ style=\"filled\", fillcolor=\"red\" ];" 
              << std::endl;
  }

  graph_traits < graph_t >::edge_iterator ei, ei_end;
  for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
    std::cout << (char)(source(*ei, g) + 'A') << " -- " 
              << (char)(target(*ei, g) + 'A')
              << "[label=\"" << component[*ei] << "\"]\n";
  std::cout << "}\n";

  return 0;
}
