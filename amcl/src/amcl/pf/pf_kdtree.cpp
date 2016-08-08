#include <cassert>
#include <cmath>
#include <cstdlib>

#include <memory>
#include <string>
#include <vector>

#include "pf_vector.hpp"
#include "pf_kdtree.hpp"

#include <cstring>
#include <iostream>
#define __FILENAME__ \
  (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#define TRACE_FUNC                                                     \
  do {                                                                 \
    std::cout << __FILENAME__ << ":" << __LINE__ << " in " << __func__ \
              << std::endl;                                            \
  } while (0);

#define TRACE_FUNC_ENTER                                               \
  do {                                                                 \
    std::cout << __FILENAME__ << ":" << __LINE__ << " in " << __func__ \
              << ": [ENTER]" << std::endl;                             \
  } while (0);

#define TRACE_FUNC_EXIT                                                \
  do {                                                                 \
    std::cout << __FILENAME__ << ":" << __LINE__ << " in " << __func__ \
              << ": [EXIT]" << std::endl;                              \
  } while (0);


using namespace amcl::pf;
using namespace amcl::pf::kdtree;

TreePtr Tree::CreateTree(size_t max_size) {
  TRACE_FUNC_ENTER
  TreePtr tree = std::make_shared<Tree>();

  // [FIXME] should cell size be 0.05? should it correspond to map resolution?
  tree->size_[0] = 0.50;
  tree->size_[1] = 0.50;
  tree->size_[2] = 10 * M_PI / 180;

  tree->root_ = nullptr;

  tree->node_max_count = max_size;

  /* [FIXME] is memory reservation necessary?
  // fill the node queue
  tree->nodes_.resize(tree->node_max_count);
  for (auto node_ptr : tree->nodes_) {
    node_ptr = std::make_shared<Node>();
  }
  */

  tree->leaf_count_ = 0;
  TRACE_FUNC_EXIT
  return tree;
}

void Tree::clear() {
  TRACE_FUNC_ENTER
  root_.reset();
  leaf_count_ = 0;
  TRACE_FUNC_EXIT
}

// Insert a pose into the tree.
// void pf_kdtree_insert(pf_kdtree_t *self, pf_vector_t pose, double value)
void Tree::insert(const Pose &pose, const double &value) {
  //TRACE_FUNC_ENTER
  Key key;

  key[0] = floor(pose.v[0] / size_[0]);
  key[1] = floor(pose.v[1] / size_[1]);
  key[2] = floor(pose.v[2] / size_[2]);

  // self->root = pf_kdtree_insert_node(self, NULL, self->root, key, value);
  NodePtr null_node = nullptr;
  insertNode(null_node, root_, key, value);
  //TRACE_FUNC_EXIT
}

////////////////////////////////////////////////////////////////////////////////
// Determine the probability estimate for the given pose. TODO: this
// should do a kernel density estimate rather than a simple histogram.
// double pf_kdtree_get_prob(pf_kdtree_t *self, pf_vector_t pose)
double Tree::getProb(const Pose &pose) {
  TRACE_FUNC_ENTER
  Key key;

  key[0] = floor(pose.v[0] / size_[0]);
  key[1] = floor(pose.v[1] / size_[1]);
  key[2] = floor(pose.v[2] / size_[2]);

  NodePtr node = findNode(root_, key);

  if (node == nullptr) {
    return 0.0;
  }

  return node->value;
  TRACE_FUNC_EXIT
}

////////////////////////////////////////////////////////////////////////////////
// Determine the cluster label for the given pose
// int pf_kdtree_get_cluster(pf_kdtree_t *self, pf_vector_t pose)
int Tree::getCluster(const Pose &pose) {
  //TRACE_FUNC_ENTER
  // int key[3];
  Key key;
  // pf_kdtree_node_t *node;

  key[0] = floor(pose.v[0] / size_[0]);
  key[1] = floor(pose.v[1] / size_[1]);
  key[2] = floor(pose.v[2] / size_[2]);

  NodePtr node = findNode(root_, key);

  if (node == nullptr) {
    //TRACE_FUNC_EXIT
    return -1;
  }

  //TRACE_FUNC_EXIT
  return node->cluster;
}

////////////////////////////////////////////////////////////////////////////////
// Compare keys to see if they are equal
// int pf_kdtree_equal(pf_kdtree_t *self, int key_a[], int key_b[])
bool Tree::isEqual(const Key &key_a, const Key &key_b) {
  // double a, b;

  if (key_a[0] != key_b[0]) return false;
  if (key_a[1] != key_b[1]) return false;
  if (key_a[2] != key_b[2]) return false;

  /* TODO: make this work (pivot selection needs fixing, too)
  // Normalize angles
  a = key_a[2] * self->size[2];
  a = atan2(sin(a), cos(a)) / self->size[2];
  b = key_b[2] * self->size[2];
  b = atan2(sin(b), cos(b)) / self->size[2];

  if ((int) a != (int) b)
    return 0;
  */
  return true;
}

// Insert a node into the tree
// pf_kdtree_node_t *pf_kdtree_insert_node(
//  pf_kdtree_t *self, pf_kdtree_node_t *parent,
//  pf_kdtree_node_t *node, int key[], double value)
NodePtr Tree::insertNode(NodePtr &parent, NodePtr &node,
                      const Key &key, const double &value) {
  //TRACE_FUNC_ENTER
  //std::cout << std::endl;
  //std::cout << "insertNode: key(" << key[0] << "," << key[1] << "," << key[2] << "), value(" << value << ")" << std::endl;
  // If the node doesnt exist yet...
  if (!node) {
    //std::cout << "!node" << std::endl;
    // [FIXME] should it be an assert statement? or just return false?
    assert(nodes_.size() < node_max_count);

    // node = self->nodes + self->node_count++;
    node = std::make_shared<Node>();
    // memset(node, 0, sizeof(pf_kdtree_node_t));

    node->leaf = true;

    if (!parent) {
      node->depth = 0;
    } else {
      node->depth = parent->depth + 1;
    }

    for (size_t i = 0; i < 3; i++) { node->key[i] = key[i]; }

    node->value = value;
    leaf_count_ += 1;
    nodes_.push_back(node);

  } else if (node->leaf) {
    // std::cout << "node is LEAF" << std::endl;
    // If the node exists, and it is a leaf node...

    // If the keys are equal, increment the value
    if (isEqual(key, node->key)) {
      // std::cout << "key(" << key[0] << "," << key[1] << "," << key[2] << ") == "
      //           << "node->key(" << node->key[0] << "," << node->key[1] << "," << node->key[2] << ")" << std::endl;
      // std::cout << "Increase Value" << std::endl;
      node->value += value;
    } else {
      //std::cout << "key(" << key[0] << "," << key[1] << "," << key[2] << ") != "
      //          << "node->key(" << node->key[0] << "," << node->key[1] << "," << node->key[2] << ")" << std::endl;
      // The keys are not equal, so split this node

      // Find the dimension with the largest variance and do a mean

      // split
      int max_split = 0;
      node->pivot_dim = -1;
      for (size_t i = 0; i < 3; i++) {
        int split = abs(key[i] - node->key[i]);
        if (split > max_split) {
          max_split = split;
          node->pivot_dim = i;
        }
      }
      assert(node->pivot_dim >= 0);

      node->pivot_value =
          (key[node->pivot_dim] + node->key[node->pivot_dim]) / 2.0;

      std::shared_ptr<Node> new_node_0 = nullptr;
      std::shared_ptr<Node> new_node_1 = nullptr;
      if (key[node->pivot_dim] < node->pivot_value) {
        node->children[0] = insertNode(node, new_node_0, key, value);
        node->children[1] =
            insertNode(node, new_node_1, node->key, node->value);
      } else {
        node->children[0] =
            insertNode(node, new_node_0, node->key, node->value);
        node->children[1] = insertNode(node, new_node_1, key, value);
      }

      node->leaf = false;
      leaf_count_ -= 1;
    }
  } else {
    // std::cout << "node is PARENT" << std::endl;
    // If the node exists, and it has children...
    assert(node->children[0] != nullptr);
    assert(node->children[1] != nullptr);

    if (key[node->pivot_dim] < node->pivot_value) {
      insertNode(node, node->children[0], key, value);
    } else {
      insertNode(node, node->children[1], key, value);
    }
  }
  //TRACE_FUNC_EXIT
  return node;
}

////////////////////////////////////////////////////////////////////////////////
// Recursive node search
// pf_kdtree_node_t *pf_kdtree_find_node(pf_kdtree_t *self, pf_kdtree_node_t
// *node, int key[])
NodePtr Tree::findNode(const NodePtr entry_node, const Key &key) {
  //TRACE_FUNC_ENTER
  if (entry_node->leaf) {
    // printf("find  : leaf %p %d %d %d\n", node, node->key[0], node->key[1],
    // node->key[2]);
    //TRACE_FUNC
    // If the keys are the same...
    if (isEqual(key, entry_node->key))
      return entry_node;
    else
      return nullptr;
  } else {
    // printf("find  : brch %p %d %f\n", node, node->pivot_dim,
    // node->pivot_value);
    assert(entry_node->children[0] != nullptr);
    assert(entry_node->children[1] != nullptr);

    // If the keys are different...
    if (key[entry_node->pivot_dim] < entry_node->pivot_value)
      return findNode(entry_node->children[0], key);
    else
      return findNode(entry_node->children[1], key);
  }
  // TRACE_FUNC_EXIT
  return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
// Recursive node printing
/*
void pf_kdtree_print_node(pf_kdtree_t *self, pf_kdtree_node_t *node)
{
  if (node->leaf)
  {
    printf("(%+02d %+02d %+02d)\n", node->key[0], node->key[1], node->key[2]);
    printf("%*s", node->depth * 11, "");
  }
  else
#ifdef INCLUDE_RTKGUI

////////////////////////////////////////////////////////////////////////////////
// Draw the tree
void pf_kdtree_draw(pf_kdtree_t *self, rtk_fig_t *fig) {
  if (self->root != NULL) pf_kdtree_draw_node(self, self->root, fig);
  return;
}

////////////////////////////////////////////////////////////////////////////////
// Recursively draw nodes
void pf_kdtree_draw_node(pf_kdtree_t *self, pf_kdtree_node_t *node,
                         rtk_fig_t *fig) {
  double ox, oy;
  char text[64];

  if (node->leaf) {
    ox = (node->key[0] + 0.5) * self->size[0];
    oy = (node->key[1] + 0.5) * self->size[1];

    rtk_fig_rectangle(fig, ox, oy, 0.0, self->size[0], self->size[1], 0);

    // snprintf(text, sizeof(text), "%0.3f", node->value);
    // rtk_fig_text(fig, ox, oy, 0.0, text);

    snprintf(text, sizeof(text), "%d", node->cluster);
    rtk_fig_text(fig, ox, oy, 0.0, text);
  } else {
    assert(node->children[0] != NULL);
    assert(node->children[1] != NULL);
    pf_kdtree_draw_node(self, node->children[0], fig);
    pf_kdtree_draw_node(self, node->children[1], fig);
  }

  return;
}

#endif

  {
    printf("(%+02d %+02d %+02d) ", node->key[0], node->key[1], node->key[2]);
    pf_kdtree_print_node(self, node->children[0]);
    pf_kdtree_print_node(self, node->children[1]);
  }
  return;
}
*/

////////////////////////////////////////////////////////////////////////////////
// Cluster the leaves in the tree
//void pf_kdtree_cluster(pf_kdtree_t *self)
void Tree::cluster()
{
  TRACE_FUNC_ENTER
  std::vector<NodePtr> queue;
TRACE_FUNC
  // Put all the leaves in a queue
  // debug only
  std::cout << "nodes_.size = " << nodes_.size() << std::endl;
  for (auto node : nodes_) {
    if (node->leaf) {
      node->cluster = -1;
      assert(queue.size() < nodes_.size());
      queue.push_back(node);

      // TESTING; remove
      assert(node == findNode(root_, node->key));
    }
  }
TRACE_FUNC
  int cluster_count = 0;

  // Do connected components for each node
  // [TODO] should the queue be reversed first?
  for (auto node : queue) {
    // If this node has already been labelled, skip it
    if (node->cluster >= 0) continue;
    // Assign a label to this cluster
    node->cluster = cluster_count++;

    // Recursively label nodes in this cluster
    clusterNode(node, 0);
  }
TRACE_FUNC_EXIT
  return;
}

////////////////////////////////////////////////////////////////////////////////
// Recursively label nodes in this cluster
//void pf_kdtree_cluster_node(pf_kdtree_t *self, pf_kdtree_node_t *node, int depth)
void Tree::clusterNode(const NodePtr node, const int &depth)
{
  int i;
  Key nkey;

  for (i = 0; i < 3 * 3 * 3; i++) {
    // [FIXME] what's following code piece doing?
    nkey[0] = node->key[0] + (i / 9) - 1;
    nkey[1] = node->key[1] + ((i % 9) / 3) - 1;
    nkey[2] = node->key[2] + ((i % 9) % 3) - 1;

    NodePtr nnode = findNode(root_, nkey);
    if (nnode == nullptr) continue;

    assert(nnode->leaf);

    // This node already has a label; skip it.  The label should be
    // consistent, however.
    if (nnode->cluster >= 0) {
      assert(nnode->cluster == node->cluster);
      continue;
    }

    // Label this node and recurse
    nnode->cluster = node->cluster;

    clusterNode(nnode, depth + 1);
  }
  return;
}
