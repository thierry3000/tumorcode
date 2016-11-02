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
#ifndef _KD_TREE_H
#define _KD_TREE_H

#include "helpers-vec.h"
#include "mempool.h"
#include "timer.h"

#include <algorithm>
#include <vector>

#include "tbb/parallel_sort.h"
#include "tbb/task_group.h"
#include "tbb/parallel_reduce.h"
#include "omp.h"

#define KDDEBUG(x) // x



template<int d, class T, class BT, class GetBBox>
class KdTree
{
  KdTree(const KdTree&);
  KdTree& operator=(const KdTree&);

public:
  typedef BBoxd<BT,d> BBox;

  KDDEBUG(static inline void debugOut(int level, const string &s)
  {
    std::cout << string(level*2,' ') << s << std::endl;
  })

  struct Node : public FixedSizeAllocated<Node>
  {
    BBox bbox;
    bool isLeaf;
    union {
      struct {
        Node *left, *right;
      } inner;
      struct {
        int a, b;
      } leaf;
    };
  public:
    enum { LEAF, INNER };
    Node(const BBox &bb, int type) : isLeaf(type == LEAF), bbox(bb)
    {
      inner.left = inner.right = NULL;
    }
    ~Node()
    {
      if (!isLeaf)
      {
        delete inner.left;
        delete inner.right;
      }
    };
  };


  struct DataCompare
  {
    const GetBBox &getBBox;
    int axis;
    DataCompare(int axis, const GetBBox &getBBox) : axis(axis), getBBox(getBBox) {}
    bool operator()(const T &a, const T &b) const
    {
      const BBox bba = getBBox(a);
      const BBox bbb = getBBox(b);
      int ca = (bba.max[axis]+bba.min[axis])/2;
      int cb = (bbb.max[axis]+bba.min[axis])/2;
      return ca < cb;
    }
  };

  struct ParallelBBoxFunc
  {
    BBox b;
    const KdTree* tree;
    ParallelBBoxFunc(const KdTree* tree) : tree(tree) {}
    ParallelBBoxFunc(ParallelBBoxFunc &other, tbb::split) : tree(other.tree) {}
    void operator()(const tbb::blocked_range<int> &r)
    {
      BBox tmp = b;
      for(int i=r.begin(); i!=r.end(); ++i) {
        tmp.Add(tree->getBBox(tree->objects[i]));
      }
      b = tmp;
    }
    void join( ParallelBBoxFunc &rhs ) { b.Add(rhs.b); }
  };

  struct RecBuildFunctor
  {
    KdTree* tree;
    int start, end, level;
    BBox box;
    Node* &node;
    RecBuildFunctor(KdTree* tree, int start, int end, int level, const BBox &box, Node* &node) :
      tree(tree), start(start), end(end), level(level), box(box), node(node)
      {}
    bool use_multithreading() const { return omp_get_max_threads()>1; }
    bool use_mt_pernode() const { return use_multithreading() && ((1 << level) <= omp_get_num_procs()); }
    
    std::pair<Node*,Node*> split(int splitIndex, const BBox &subboxa, const BBox &subboxb)
    {
      Node *a = NULL, *b = NULL;
      RecBuildFunctor builda(tree, start, splitIndex, level+1, subboxa, a);
      RecBuildFunctor buildb(tree, splitIndex, end, level+1, subboxb, b);
      if (use_multithreading())
      {
        tbb::task_group g;
        g.run(builda);
        g.run(buildb);
        g.wait();
      }
      else
      {
        builda();
        buildb();
      }
      return std::make_pair(a, b);
    }

    Node* operator()()
    {
      node = NULL;
      const int cnt = end-start;
      myAssert(cnt >= 0);
      if( cnt == 0 ) return NULL;

      const int MIN_CNT = 4;
      const int MAX_LEVEL = 100;
      const float MIN_EFFICIENCY = 0.25;

      if( cnt>MIN_CNT && level<MAX_LEVEL)
      {
        int axis = MajorAxis<float,d>(box.max - box.min);
        double box_len = (box.max[axis] - box.min[axis]);
        if (box_len > tree->min_box_len)
        {
          partitionData(axis);
          //int trialIndicies[3] = { start + cnt * 0.5, start + cnt * 0.33, start + cnt *0.66 };
          //FOR_EACH(int splitIndex, trialIndicies, 3)
          int splitIndex = start+cnt*0.5;
          {
            BBox subboxa = computeBBox(start, splitIndex);
            BBox subboxb = computeBBox(splitIndex, end);
            BBox intersection = BBox(subboxa).Clip(subboxb);
            double q = double(intersection.max[axis] - intersection.min[axis])/box_len;
            float efficiency = 1. - q;

            //if(efficiency > MIN_EFFICIENCY)
            {
              KDDEBUG(debugOut(level, boost::str(boost::format("level %i, split %i of %i") % level % splitIndex % cnt));)
              Node *nodea, *nodeb; 
              tie(nodea, nodeb) = split(splitIndex, subboxa, subboxb);
              if(nodea && nodeb)
              {
                node = new Node(box, Node::INNER);
                node->inner.left = nodea;
                node->inner.right = nodeb;
                return node;
              }
              else if(nodea)
              {
                KDDEBUG(debugOut(level, "(single left)");)
                node = nodea;
                return node;
              }
              else if(nodeb)
              {
                KDDEBUG(debugOut(level, "(single right)");)
                node = nodeb;
                return node;
              }
            }
          }
        }
      }

      KDDEBUG(debugOut(level, boost::str(boost::format("level %i, leaf %i") % level % cnt));)
      // if splitting failed
      node = new Node(box, Node::LEAF);
      node->leaf.a = start;
      node->leaf.b = end;
      return node;
    }

    int partitionData(int axis)
    {
      if (use_mt_pernode())
        std::sort(tree->objects.begin()+start, tree->objects.begin()+end, DataCompare(axis, tree->getBBox));
      else
       tbb::parallel_sort(tree->objects.begin()+start, tree->objects.begin()+end, DataCompare(axis, tree->getBBox));
      return (end-start)/2;
    }

    BBox computeBBox(int start, int end) const
    {
      ParallelBBoxFunc f(tree);
      if (use_mt_pernode())
        tbb::parallel_reduce(tbb::blocked_range<int>(start, end), f);
      else
        f(tbb::blocked_range<int>(start, end));
      return f.b;
    }

    RecBuildFunctor(KdTree* tree, Node* &node) : start(0), end(tree->objects.size()), node(node), level(0), tree(tree)
    {
      box = computeBBox(0, tree->objects.size());
    }
  };

  void initialize(std::vector<T> &objects_, double min_box_len_, const GetBBox &getBBox_)
  {
    KDDEBUG(my::Time t_;)
    KDDEBUG(printf("kdtree init with %i objects\n", objects_.size());)
    getBBox = getBBox_;
    min_box_len = min_box_len_;
    objects_.swap(objects);
    RecBuildFunctor(this, rootNode)();
    KDDEBUG(printf("tree construct time: %f\n",(my::Time()-t_).to_ms());)
  }

  int debug_search;

  void regionSearchRecursive(const BBox &bbox, std::vector<T> &res, int level, Node* node, bool check_bbox) const
  {
    if (!node) return;
    if(!bbox.Overlaps(node->bbox)) return;
    if(node->isLeaf)
    {
      //std::cout << std::string(level,' ') << node->bbox << std::endl;
      for(int i=node->leaf.a; i<node->leaf.b; ++i)
      {
        if (check_bbox)
        {
          const BBox b = getBBox(objects[i]);
          if(b.Overlaps(bbox))
            res.push_back(objects[i]);
        }
        else res.push_back(objects[i]);
      }
    }
    else
    {
      regionSearchRecursive(bbox, res, level+1, node->inner.left, check_bbox);
      regionSearchRecursive(bbox, res, level+1, node->inner.right, check_bbox);
    }
  }

  void regionSearch(const BBox &bbox, std::vector<T> &res, bool check_bbox = true) const
  {
    //printf("new search\n");
    res.clear();
    regionSearchRecursive(bbox, res, 0, rootNode, check_bbox);
    //printf("search set: %i\n", res.size());
  }

  void pointSearch(const Int3 &p, std::vector<T> &res, bool check_bbox = true) const
  {
    regionSearch(BBox().Add(p), res, check_bbox);
  }

  void clear()
  {
    delete rootNode;
    rootNode = NULL;
    clear_and_free_memory<>(objects);
  }

  void write_gnuplot(std::ostream &os, Node* node, int level) const
  {
    if (!node) return;
    const BBox &bb = node->bbox;
    double w = double(((std::size_t)node) * 167845367 % 50)/100.0 + level;
    ::write_gnuplot(os, bb, w);

    if(!node->isLeaf)
    {
      write_gnuplot(os, node->inner.left, level+1);
      write_gnuplot(os, node->inner.right, level+1);
    }
    else
    {
      for (int i=node->leaf.a; i<node->leaf.b; ++i)
      {
        ::write_gnuplot(os, getBBox(objects[i]), w+0.01);
      }
    }
  }

  void write_gnuplot(std::ostream &os) const
  {
    write_gnuplot(os, rootNode, 0);
  }

  std::size_t memory_consumption(const Node* n) const
  {
    std::size_t m = sizeof(Node);
    if (!n->isLeaf) {
      m += memory_consumption(n->inner.left) + memory_consumption(n->inner.right);
    }
    return m;
  }

  std::size_t memory_consumption() const
  {
    std::size_t m = sizeof(*this);
    m += objects.capacity() * sizeof(T);
    m += memory_consumption(rootNode);
    return m;
  }

  KdTree() : rootNode(NULL) 
  {
  }

  ~KdTree() 
  {
    clear();
  }

private:
  Node* rootNode;
  std::vector<T> objects;
  double min_box_len;
  GetBBox getBBox;
};

#endif