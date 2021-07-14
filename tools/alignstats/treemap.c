#include "treemap.h"
#include "err.h"
#include "logging.h"
#include <stdio.h>

/**
 * Basic binary search tree-based map.
 */

/**
 * Find where to insert a new node with a given key.
 */
tree_node_t *_insert_pos(tree_node_t *root, tree_node_key_t key)
{
    /* Tree empty */
    tree_node_t *curr_node = NULL;

    if (root != NULL) {
        /* Search for insertion position */
        curr_node = root;

        while (key != curr_node->key) {
            if (key < curr_node->key) {        /* Go left */
                if (curr_node->left == NULL) { /* At leaf node */
                    break;
                }
                curr_node = curr_node->left;
            } else {                            /* Go right */
                if (curr_node->right == NULL) { /* At leaf node */
                    break;
                }
                curr_node = curr_node->right;
            }
        }
    }

    return curr_node;
}

/**
 * Create new tree node.
 */
tree_node_t *tree_node_init(tree_node_key_t key, tree_node_value_t value)
{
    tree_node_t *node = calloc(1, sizeof(tree_node_t));
    die_on_alloc_fail(node);

    node->key = key;
    node->value = value;
    node->parent = NULL;
    node->left = NULL;
    node->right = NULL;

    return node;
}

/**
 * Free *tree_node_t from memory.
 */
void tree_node_destroy(tree_node_t *node)
{
    free(node);
}

/* Tree map structure */

/**
 * Create a new tree.
 */
tree_map_t *tree_map_init()
{
    tree_map_t *tree = calloc(1, sizeof(tree_map_t));
    die_on_alloc_fail(tree);

    tree->root = NULL;
    tree->num_nodes = 0;

    return tree;
}

/**
 * Recursive function used by tree_map_destroy() to free each node.
 */
void _tree_map_destroy(tree_node_t *node)
{
    if (node != NULL) {
        _tree_map_destroy(node->left);
        _tree_map_destroy(node->right);
        tree_node_destroy(node);
    }
}

/**
 * Free all nodes in the tree from heap memory.
 */
void tree_map_destroy(tree_map_t *tree)
{
    if (tree != NULL) {
        _tree_map_destroy(tree->root);
        free(tree);
    }
}

/* Tree map functions */

/**
 * Put key into tree with value. If a node with key already exists, the node's
 * value is reset to value.
 */
void tree_map_set(tree_map_t *tree, tree_node_key_t key,
                  tree_node_value_t value)
{
    tree_node_t *root = tree->root;

    /* Empty tree: new node becomes root */
    if (root == NULL) {
        ++tree->num_nodes;
        tree->root = tree_node_init(key, value);
    } else {
        /* Insert as leaf node */
        tree_node_t *curr_node = _insert_pos(root, key);
        tree_node_t *node;

        if (key == curr_node->key) {
            /* Node already in tree: reassign value */
            curr_node->value = value;
        } else {
            /* If node with key not already in tree */
            ++tree->num_nodes;

            /* Create new node */
            if ((node = tree_node_init(key, value)) != NULL) {
                /* Place new node in tree */
                if (key < curr_node->key) {
                    curr_node->left = node;
                } else {
                    curr_node->right = node;
                }

                node->parent = curr_node;
            }
        }
    }
}

/**
 * Get the node with corresponding key. If node with key is not in tree, NULL
 * is returned.
 */
tree_node_t *tree_map_get(tree_map_t *tree, tree_node_key_t key)
{
    tree_node_t *node = _insert_pos(tree->root, key);

    return (node != NULL && node->key == key) ? node : NULL;
}

/**
 * Returns node with next smallest key, or NULL if no such node exits.
 */
/*
tree_node_t *_prev_node(tree_node_t *node)
{
    tree_node_t *prev_node = node->left;

    while (prev_node != NULL && prev_node->right != NULL)
        prev_node = prev_node->right;

    return prev_node;
}
*/

/**
 * Remove node referenced by key from tree.
 */
/*
void tree_map_remove(tree_map_t *tree, tree_node_key_t key)
{
    tree_node_t *node, *prev_node;

    // if key not in tree
    if (tree == NULL || (node = _insert_pos(tree->root, key)) == NULL ||
        node->key != key)
        return;

    //if (node->left != NULL && node->right != NULL)
    if (node->left != NULL) {
        if (node->right != NULL) {
            // Node has two children
            // Get node with next smallest key
            prev_node = _prev_node(node);

            // Copy key and value from prev_node into node
            node->key = prev_node->key;
            node->value = prev_node->value;

            // Transfer prev_node subtree to node
            if (prev_node->parent == node) {
                prev_node->parent->left = prev_node->left;
            } else {
                prev_node->parent->right = prev_node->left;
            }

            if (prev_node->left != NULL) {
                prev_node->left->parent = prev_node->parent;
            }

            // Set prev_node as node to be freed
            node = prev_node;
        } else {
            // Node has only left child
            if (node == tree->root) {
                tree->root = node->left;
            } else {
                node->parent->left = node->left;
            }
            node->left->parent = node->parent;
        }
    } else {
        //else if (node->left == NULL && node->right != NULL)
        if (node->right != NULL) {
            // Node has only right child
            if (node == tree->root) {
                tree->root = node->right;
            } else {
                node->parent->left = node->right;
            }

            node->right->parent = node->parent;
        } else {
            // Node has no children
            if (node == tree->root) {
                tree->root = NULL;
            } else if (node == node->parent->left) {
                node->parent->left = NULL;
            } else {
                node->parent->right = NULL;
            }
        }
    }

    tree_node_destroy(node);
}
*/

/**
 * Keyset helper function.
 */
int _tree_map_keyset(tree_node_t *node, tree_node_key_t *keyset, int i)
{
    int rv = i;

    if (node != NULL) {
        /* Collect keys from left subtree */
        i = _tree_map_keyset(node->left, keyset, i);

        /* Collect this node's key */
        keyset[i++] = node->key;

        /* Collect keys from right subtree and return */
        rv = _tree_map_keyset(node->right, keyset, i);
    }

    return rv;
}

/**
 * Returns new tree_node_key_t array of keys in sorted order.
 */
tree_node_key_t *tree_map_keyset(tree_map_t *tree)
{
    tree_node_key_t *keyset = calloc(tree->num_nodes, sizeof(tree_node_key_t));
    die_on_alloc_fail(keyset);

    _tree_map_keyset(tree->root, keyset, 0);

    return keyset;
}

/**
 * Sets node referenced by *key to **node.
 * Returns whether operation was successful.
 */
bool tree_map_set_node(tree_node_t **node, tree_map_t *tree, tree_node_key_t key)
{
    bool rv = true;
    tree_node_t *tmp_node = tree_map_get(tree, key);

    if (tmp_node != NULL) {
        *node = tmp_node;
    } else {
        log_error("Failed to retrieve node from tree map. This should not happen.");
        rv = false;
    }

    return rv;
}

/**
 * Sets the keyset of *tree to **keyset.
 * Returns whether operation was successful.
 */
bool tree_map_set_keyset(tree_node_key_t **keyset, tree_map_t *tree)
{
    bool rv = true;
    tree_node_key_t *tmp_keyset = tree_map_keyset(tree);

    if (tmp_keyset != NULL) {
        *keyset = tmp_keyset;
    } else {
        log_error("Failed to retrieve keyset for tree map. This should not happen.");
        rv = false;
    }

    return rv;
}
