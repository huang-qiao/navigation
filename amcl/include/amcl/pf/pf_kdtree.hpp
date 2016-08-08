#pragma once

#include "pf_vector.hpp"

#include <memory>
#include <vector>

namespace amcl {
namespace pf {
namespace kdtree {

// forward declaration
struct Node;
class Tree;
using NodePtr = std::shared_ptr<Node>;
using TreePtr = std::shared_ptr<Tree>;
using Key = int[3];

// Info for a node in the tree
struct Node
{

  Node():leaf(false), depth(0)
  {
    children[0] = nullptr;
    children[1] = nullptr;
  }

  // Depth in the tree
  bool leaf;
  size_t depth;

  // Pivot dimension and value
  size_t pivot_dim;
  double pivot_value;

  // The key for this node
  Key key;

  // The value for this node
  double value;

  // The cluster label (leaf nodes)
  int cluster;

  // Child nodes
  NodePtr children[2];
};

// A kd tree
class Tree
{
public:
  // Create a tree
  static TreePtr CreateTree(size_t max_size);

  // Clear all entries from the tree
  void clear();

  // Insert a pose into the tree
  void insert(const Pose &pose, const double &value);

  // Cluster the leaves in the tree
  void cluster();

  // Determine the probability estimate for the given pose
  double getProb(const Pose &pose);

  // Determine the cluster label for the given pose
  int getCluster(const Pose &pose);

  inline std::vector<NodePtr> getNodes() { return nodes_; }
  inline int leafCount() { return leaf_count_; }
  inline size_t nodeCount() { return nodes_.size(); }
  bool isEqual(const Key& key_a, const Key& key_b);

  // Insert a node into the tree
  NodePtr insertNode(NodePtr &parent, NodePtr &node, const Key &key, const double &value);

  // Recursive node search, return nullptr if not found
  NodePtr findNode(const NodePtr entry_node, const Key &key);

  void clusterNode(const NodePtr node, const int& depth);

private:
  // The number of leaf nodes in the tree
  int leaf_count_;
  std::vector<NodePtr> nodes_;

  // Cell size
  double size_[3];

  // The root node of the tree
  NodePtr root_;

  // The number of nodes in the tree
  size_t node_max_count;

}; // class Tree

using TreePtr = std::shared_ptr<Tree>;

} // namespace kdtree
} // namespace pf
} // namespace amcl
