#ifndef _TREEMAP_H
#define _TREEMAP_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

/* Both key and value data types should be a base numerical type. */
typedef int32_t tree_node_key_t;
typedef int32_t tree_node_value_t;

/* Tree map node structure */
struct tree_node {
    tree_node_key_t key;      /* Key */
    tree_node_value_t value;  /* Value */
    struct tree_node *parent; /* Parent node; root node's parent == NULL */
    struct tree_node *left;   /* Left subtree */
    struct tree_node *right;  /* Right subtree */
};
typedef struct tree_node tree_node_t;

tree_node_t *tree_node_init(tree_node_key_t key, tree_node_value_t value);
void tree_node_destroy(tree_node_t *node);

/**
 * Tree map structure
 * Designed so that the root can be moved in place.
 * Empty tree: *root == NULL
 * Non-empty tree: *root == root_node, root_node->parent == NULL
 */
struct tree_map {
    tree_node_t *root; /* Root node */
    size_t num_nodes;  /* Total number of nodes in tree */
};
typedef struct tree_map tree_map_t;

tree_map_t *tree_map_init();
void tree_map_destroy(tree_map_t *tree);

/* Tree map functions */
void tree_map_set(tree_map_t *tree, tree_node_key_t key,
                  tree_node_value_t value);
tree_node_t *tree_map_get(tree_map_t *tree, tree_node_key_t key);
tree_node_key_t *tree_map_keyset(tree_map_t *tree);
bool tree_map_set_node(tree_node_t **node, tree_map_t *tree, tree_node_key_t key);
bool tree_map_set_keyset(tree_node_key_t **keyset, tree_map_t *tree);

#endif /* _TREEMAP_H */
