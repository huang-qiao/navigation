#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "pf_kdtree.hpp"
#include "pf_vector.hpp"
#include "util.hpp"

#ifdef INCLUDE_RTKGUI

// Recursively draw nodes
static void pf_kdtree_draw_node(pf_kdtree_t *self, pf_kdtree_node_t *node,
                                rtk_fig_t *fig);

#endif

KdTree::KdTree(const int &max_size) {
  size[0] = 0.50;
  size[1] = 0.50;
  size[2] = (10 * M_PI / 180);

  root = NULL;

  node_count = 0;
  node_max_count = max_size;
  nodes = (Node *)malloc(node_max_count * sizeof(Node));

  leaf_count = 0;
}

KdTree::~KdTree() { free(nodes); }

////////////////////////////////////////////////////////////////////////////////
/* Destroy a tree
void pf_kdtree_delete(KdTree *self) {
  free(self->nodes);
  free(self);
  return;
}
*/
////////////////////////////////////////////////////////////////////////////////
// Clear all entries from the tree
void KdTree::clear() {
  root = NULL;
  leaf_count = 0;
  node_count = 0;
}

////////////////////////////////////////////////////////////////////////////////
// Insert a pose into the tree.
void KdTree::insert(Pose pose, double value) {
  int key[3];

  key[0] = floor(pose.v[0] / size[0]);
  key[1] = floor(pose.v[1] / size[1]);
  key[2] = floor(pose.v[2] / size[2]);

  root = insertNode(NULL, root, key, value);

  /* Test code
  //
  // printf("find %d %d %d\n", key[0], key[1], key[2]);
  std::cout << "find " << key[0] << " " << key[1] << " " << key[2] << std::endl;
  assert(pf_kdtree_find_node(self, root, key) != NULL);

  // pf_kdtree_print_node(self, root);

  printf("\n");

  for (int i = 0; i < node_count; i++) {
    pf_kdtree_node_t *node = nodes + i;
    if (node->leaf) {
      // printf("find %d %d %d\n", node->key[0], node->key[1], node->key[2]);
      std::cout << "find " << key[0] << " " << key[1] << " " << key[2]
                << std::endl;
      assert(pf_kdtree_find_node(self, root, node->key) == node);
    }
  }
  std::cout << std::endl << std::endl; // printf("\n\n");
  */

  return;
}

////////////////////////////////////////////////////////////////////////////////
// Determine the probability estimate for the given pose. TODO: this
// should do a kernel density estimate rather than a simple histogram.
double KdTree::getProb(Pose pose) {
  int key[3];
  Node *node;

  key[0] = floor(pose.v[0] / size[0]);
  key[1] = floor(pose.v[1] / size[1]);
  key[2] = floor(pose.v[2] / size[2]);

  node = findNode(root, key);
  if (node == NULL)
    return 0.0;
  return node->value;
}

////////////////////////////////////////////////////////////////////////////////
// Determine the cluster label for the given pose
int KdTree::getCluster(Pose pose) {
  int key[3];
  Node *node;

  key[0] = floor(pose.v[0] / size[0]);
  key[1] = floor(pose.v[1] / size[1]);
  key[2] = floor(pose.v[2] / size[2]);

  node = findNode(root, key);
  if (node == NULL)
    return -1;
  return node->cluster;
}

////////////////////////////////////////////////////////////////////////////////
// Compare keys to see if they are equal
bool KdTree::isEqual(int key_a[], int key_b[]) {
  // double a, b;

  if (key_a[0] != key_b[0])
    return false;
  if (key_a[1] != key_b[1])
    return false;

  if (key_a[2] != key_b[2])
    return false;

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

////////////////////////////////////////////////////////////////////////////////
// Insert a node into the tree
Node *KdTree::insertNode(Node *parent, Node *node, int key[], double value) {
  int i;
  int split, max_split;

  // If the node doesnt exist yet...
  if (node == NULL) {
    assert(node_count < node_max_count);
    node = nodes + node_count++;
    memset(node, 0, sizeof(Node));

    node->leaf_ = true;

    if (parent == NULL)
      node->depth = 0;
    else
      node->depth = parent->depth + 1;

    for (i = 0; i < 3; i++)
      node->key[i] = key[i];

    node->value = value;
    leaf_count += 1;
  }

  // If the node exists, and it is a leaf node...
  else if (node->leaf_) {
    // If the keys are equal, increment the value
    if (isEqual(key, node->key)) {
      node->value += value;
    }

    // The keys are not equal, so split this node
    else {
      // Find the dimension with the largest variance and do a mean
      // split
      max_split = 0;
      node->pivot_dim = -1;
      for (i = 0; i < 3; i++) {
        split = abs(key[i] - node->key[i]);
        if (split > max_split) {
          max_split = split;
          node->pivot_dim = i;
        }
      }
      assert(node->pivot_dim >= 0);

      node->pivot_value =
          (key[node->pivot_dim] + node->key[node->pivot_dim]) / 2.0;

      if (key[node->pivot_dim] < node->pivot_value) {
        node->children[0] = insertNode(node, NULL, key, value);
        node->children[1] = insertNode(node, NULL, node->key, node->value);
      } else {
        node->children[0] = insertNode(node, NULL, node->key, node->value);
        node->children[1] = insertNode(node, NULL, key, value);
      }

      node->leaf_ = false;
      leaf_count -= 1;
    }
  }

  // If the node exists, and it has children...
  else {
    assert(node->children[0] != NULL);
    assert(node->children[1] != NULL);

    if (key[node->pivot_dim] < node->pivot_value)
      insertNode(node, node->children[0], key, value);
    else
      insertNode(node, node->children[1], key, value);
  }

  return node;
}

////////////////////////////////////////////////////////////////////////////////
// Recursive node search
Node *KdTree::findNode(Node *node, int key[]) {
  if (node->leaf_) {
    // printf("find  : leaf %p %d %d %d\n", node, node->key[0], node->key[1],
    // node->key[2]);

    // If the keys are the same...
    if (isEqual(key, node->key))
      return node;
    else
      return NULL;
  } else {
    // printf("find  : brch %p %d %f\n", node, node->pivot_dim,
    // node->pivot_value);

    assert(node->children[0] != NULL);
    assert(node->children[1] != NULL);

    // If the keys are different...
    if (key[node->pivot_dim] < node->pivot_value)
      return findNode(node->children[0], key);
    else
      return findNode(node->children[1], key);
  }

  return NULL;
}

////////////////////////////////////////////////////////////////////////////////
// Recursive node printing
//
void KdTree::printNode(Node *node) {
#if 0
  if (node->leaf) {
    // printf("(%+02d %+02d %+02d)\n", node->key[0], node->key[1],
    // node->key[2]);
    std::cout << "(" << node->key[0] << " " << node->key[1] << " "
              << node->key[2] << ")" << std::endl;
    // printf("%*s", node->depth * 11, "");
    std::cout << "depth = " << node->depth << std::endl;
  } else {
    // printf("(%+02d %+02d %+02d) ", node->key[0], node->key[1], node->key[2]);
    std::cout << "(" << node->key[0] << " " << node->key[1] << " "
              << node->key[2] << ")" << std::endl;
    pf_kdtree_print_node(self, node->children[0]);
    pf_kdtree_print_node(self, node->children[1]);
  }
  return;
#endif
}
//

////////////////////////////////////////////////////////////////////////////////
// Cluster the leaves in the tree
void KdTree::cluster() {
  int i;
  int queue_count, cluster_count;
  Node **queue, *node;

  queue_count = 0;
  queue = (Node **)malloc(node_count * sizeof(queue[0]));

  // Put all the leaves in a queue
  for (i = 0; i < node_count; i++) {
    node = nodes + i;
    if (node->leaf_) {
      node->cluster = -1;
      assert(queue_count < node_count);
      queue[queue_count++] = node;

      // TESTING; remove
      assert(node == findNode(root, node->key));
    }
  }

  cluster_count = 0;

  // Do connected components for each node
  while (queue_count > 0) {
    node = queue[--queue_count];

    // If this node has already been labelled, skip it
    if (node->cluster >= 0)
      continue;

    // Assign a label to this cluster
    node->cluster = cluster_count++;

    // Recursively label nodes in this cluster
    clusterNode(node, 0);
  }

  free(queue);
  return;
}

////////////////////////////////////////////////////////////////////////////////
// Recursively label nodes in this cluster
void KdTree::clusterNode(Node *node, int depth) {
  int i;
  int nkey[3];
  Node *nnode;

  for (i = 0; i < 3 * 3 * 3; i++) {
    nkey[0] = node->key[0] + (i / 9) - 1;
    nkey[1] = node->key[1] + ((i % 9) / 3) - 1;
    nkey[2] = node->key[2] + ((i % 9) % 3) - 1;

    nnode = findNode(root, nkey);
    if (nnode == NULL)
      continue;

    assert(nnode->leaf_);

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

#ifdef INCLUDE_RTKGUI

////////////////////////////////////////////////////////////////////////////////
// Draw the tree
void pf_kdtree_draw(pf_kdtree_t *self, rtk_fig_t *fig) {
  if (self->root != NULL)
    pf_kdtree_draw_node(self, self->root, fig);
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
