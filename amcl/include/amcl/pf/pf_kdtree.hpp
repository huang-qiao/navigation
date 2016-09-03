#pragma once

#include "pf_vector.hpp"

// Info for a node in the tree
struct Node {
  // Depth in the tree
  bool leaf_;
  int depth;

  // Pivot dimension and value
  int pivot_dim;
  double pivot_value;

  // The key for this node
  int key[3];

  // The value for this node
  double value;

  // The cluster label (leaf nodes)
  int cluster;

  // Child nodes
  struct Node *children[2];
};

// A kd tree
struct KdTree {
  // Cell size
  double size[3];

  // The root node of the tree
  Node *root;

  // The number of nodes in the tree
  int node_count, node_max_count;
  Node *nodes;

  // The number of leaf nodes in the tree
  int leaf_count;

  KdTree(const int &max_size);
  virtual ~KdTree();

  // Clear all entries from the tree
  void clear();

  // Cluster the leaves in the tree
  void cluster();

  // Determine the probability estimate for the given pose
  double getProb(Pose pose);

  // Determine the cluster label for the given pose
  int getCluster(Pose pose);

  // Insert a pose into the tree
  void insert(Pose pose, double value);

private:
  // Recursive node search
  Node *findNode(Node *node, int key[]);

  // Recursively label nodes in this cluster
  void clusterNode(Node *node, int depth);

  // Compare keys to see if they are equal
  bool isEqual(int key_a[], int key_b[]);

  // Insert a node into the tree
  Node *insertNode(Node *parent, Node *node, int key[], double value);

  // Recursive node printing
  void printNode(Node *node);
};

// Destroy a tree
// extern void pf_kdtree_delete(KdTree *self);

#ifdef INCLUDE_RTKGUI

// Draw the tree
extern void pf_kdtree_draw(pf_kdtree_t *self, rtk_fig_t *fig);

#endif
