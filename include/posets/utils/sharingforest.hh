#pragma once

#include <unordered_map>

#include <boost/functional/hash.hpp>
#include <cassert>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <queue>
#include <ranges>
#include <stack>
#include <tuple>
#include <vector>

#include <posets/concepts.hh>

namespace posets::utils {

// We will be (micro)managing (C-style) dynamic memory for the set of nodes we
// keep in each layer. The initial number of nodes allowed is determined by
// the number below (multiplied by the number of layers = the dimension + 1).
// When the buffer is full, we double its capacity and copy everything to the
// new block of reserved memory.
#ifndef SHARINGFOREST_INIT_LAYER_SIZE
# define SHARINGFOREST_INIT_LAYER_SIZE 100UL
#endif
#ifndef SHARINGFOREST_INIT_MAX_CHILDREN
# define SHARINGFOREST_INIT_MAX_CHILDREN 10UL
#endif

  // Forward definition for the operator<<
  template <Vector>
  class sharingforest;

  template <Vector V>
  std::ostream& operator<< (std::ostream& os, const utils::sharingforest<V>& f);

  template <Vector V>
  class sharingforest {
    private:
      template <Vector V2>
      friend std::ostream& operator<< (std::ostream& os, const utils::sharingforest<V2>& f);

      size_t dim;

      struct st_node {
          typename V::value_type label;
          size_t numchild;
          size_t cbuffer_offset;
      };

      // NOTE: We will be keeping a unique table of st_nodes using an unordered
      // map, for this we need a hash function and an equivalence operation.
      // Since nodes keep offsets instead of pointers, we need both the hash and
      // equivalence operations to be able to access the unique table. Hence, we
      // keep a reference to the table in them.
      class st_hash {
          sharingforest* f;

        public:
          st_hash (sharingforest* that) : f {that} {}

          size_t operator() (const st_node& k) const {
            size_t res = std::hash<typename V::value_type> () (k.label);
            size_t* children = f->child_buffer + k.cbuffer_offset;
            for (size_t i = 0; i < k.numchild; i++)
              res ^= std::hash<size_t> () (children[i]) << (i + 1);
            return res;
          }
      };

      class st_equal {
          sharingforest* f;

        public:
          st_equal (sharingforest* that) : f {that} {}

          bool operator() (const st_node& lhs, const st_node& rhs) const {
            if (lhs.label != rhs.label or lhs.numchild != rhs.numchild)
              return false;
            size_t* lhs_children = f->child_buffer + lhs.cbuffer_offset;
            size_t* rhs_children = f->child_buffer + rhs.cbuffer_offset;
            for (size_t i = 0; i < lhs.numchild; i++)
              if (lhs_children[i] != rhs_children[i])
                return false;
            return true;
          }
      };

      std::vector<std::vector<st_node>> layers;
      size_t* child_buffer;
      size_t cbuffer_size;
      size_t cbuffer_nxt;

      // This is a cache/hash-map to get in-layer node identifiers from their
      // signature
      std::vector<std::unordered_map<st_node, size_t, st_hash, st_equal>> inverse;
      // Second cache/hash-map to check for a pair of node(-identifiers) in a layer
      // whether there is a simulation relation in the left-to-right direction
      std::vector<std::unordered_map<std::pair<size_t, size_t>, bool,
                                     boost::hash<std::pair<size_t, size_t>>>>
          simulating;
      // Two more for union and intersection
      std::vector<std::unordered_map<std::pair<size_t, size_t>, size_t,
                                     boost::hash<std::pair<size_t, size_t>>>>
          cached_union, cached_inter;

      void init (size_t dim) {
        this->dim = dim;
        layers.resize (dim + 1);

        for (size_t i = 0; i < dim + 1; i++) {
          inverse.emplace_back (SHARINGFOREST_INIT_LAYER_SIZE, st_hash (this), st_equal (this));
          simulating.emplace_back ();
          cached_union.emplace_back ();
          cached_inter.emplace_back ();
        }

        cbuffer_size = SHARINGFOREST_INIT_LAYER_SIZE * SHARINGFOREST_INIT_MAX_CHILDREN;
        child_buffer = new size_t[cbuffer_size];
        cbuffer_nxt = 0;
      }

      std::optional<size_t> has_son (st_node& node, size_t child_layer, int val) {
        if (node.numchild == 0)
          return std::nullopt;

        size_t left = 0;
        size_t right = node.numchild - 1;
        size_t* children = child_buffer + node.cbuffer_offset;

        while (left <= right) {
          size_t mid = left + ((right - left) / 2);
          assert (mid < node.numchild);
          typename V::value_type mid_val = layers[child_layer][children[mid]].label;

          if (mid_val == val)
            return mid;
          if (mid_val > val)
            left = mid + 1;
          else
            right = mid - 1;
        }

        return std::nullopt;
      }

      bool simulates (size_t n1idx, size_t n2idx, size_t layidx) {
        auto node_pair = std::make_pair (n1idx, n2idx);
        auto cached = simulating[layidx].find (node_pair);
        if (cached != simulating[layidx].end ())
          return cached->second;

        st_node n1 = layers[layidx][n1idx];
        st_node n2 = layers[layidx][n2idx];
        // If the node is in the last layer, we just check the labels
        if (layidx == this->dim)
          return n1.label >= n2.label;
        size_t* n1_children = child_buffer + n1.cbuffer_offset;
        size_t* n2_children = child_buffer + n2.cbuffer_offset;
        if (n1.label < n2.label or layers[layidx + 1][n1_children[0]].label <
                                       layers[layidx + 1][n2_children[n2.numchild - 1]].label) {
          // If the label of n1 is too small or its largest child is already smaller than n2's
          // smallest, we already know it can't simulate
          return false;
        }

        // Stack contains node and child ID of S, node and child ID of T, layer
        std::stack<std::tuple<size_t, size_t, size_t, size_t, size_t>> current_stack;
        current_stack.emplace (n1idx, 0, n2idx, 0, layidx);

        while (not current_stack.empty ()) {
          auto [n1idx, c1, n2idx, c2, layidx] = current_stack.top ();
#ifndef NDEBUG
          std::cout << "Simulation check: " << "Node " << layidx << "." << n1idx << " (c=" << c1
                    << ") vs. " << layidx << "." << n2idx << " (c=" << c2 << ")\n";
#endif
          current_stack.pop ();
          n1 = layers[layidx][n1idx];
          n2 = layers[layidx][n2idx];

          // This is our base case, no more things to check on the n2 side
          // We now claim success for n2 being simulated by n1 and update the top
          // of the stack so that when we go back up in the tree we have one less
          // branch to check on the n2 side
          if (c2 == n2.numchild or layidx == this->dim) {
            assert (c2 == n2.numchild);
            simulating[layidx][std::make_pair (n1idx, n2idx)] = true;
            if (not current_stack.empty ()) {
              auto [m1idx, d1, m2idx, d2, ell] = current_stack.top ();
              current_stack.pop ();
              current_stack.emplace (m1idx, 0, m2idx, d2 + 1, ell);
            }
          }
          // Another base case: we've iterated through all the children on the
          // n1 side and failed to find a simulating one
          else if (c1 == n1.numchild) {
#ifndef NDEBUG
            std::cout << "Another base case (layidx=" << layidx << ")\n";
#endif
            simulating[layidx][std::make_pair (n1idx, n2idx)] = false;
            if (not current_stack.empty ()) {
              auto [m1idx, d1, m2idx, d2, ell] = current_stack.top ();
              current_stack.pop ();
              current_stack.emplace (m1idx, d1 + 1, m2idx, d2, ell);
            }
          }
          // The last case is that we have two valid children indices,
          // then we have to go deeper in the product tree (lest the cache saves
          // us, or we find that the subtree will just certainly not satisfy
          // simulation)
          else {
            n1_children = child_buffer + n1.cbuffer_offset;
            n2_children = child_buffer + n2.cbuffer_offset;
            node_pair = std::make_pair (n1_children[c1], n2_children[c2]);
            cached = simulating[layidx + 1].find (node_pair);
            // Did we get lucky with the cache? then push back an updated node
            // with less obligations or keep searching on the n1 side
            if (cached != simulating[layidx + 1].end ()) {
#ifndef NDEBUG
              std::cout << "Got lucky with simulate cache!\n";
#endif
              if (cached->second)
                current_stack.emplace (n1idx, 0, n2idx, c2 + 1, layidx);
              else
                current_stack.emplace (n1idx, c1 + 1, n2idx, c2, layidx);
            }
            else {
              const st_node c1_node = layers[layidx + 1][n1_children[c1]];
              const st_node c2_node = layers[layidx + 1][n2_children[c2]];
              size_t* c1_children = child_buffer + c1_node.cbuffer_offset;
              size_t* c2_children = child_buffer + c2_node.cbuffer_offset;
              // In some cases, we can eliminate the subtree altogether
              if (c1_node.label < c2_node.label or
                  (layidx + 1 < this->dim and
                   layers[layidx + 2][c1_children[0]].label <
                       layers[layidx + 2][c2_children[c2_node.numchild - 1]].label)) {
                current_stack.emplace (n1idx, c1 + 1, n2idx, c2, layidx);
                // Alright, we have to go deeper now; and push back what we just
                // popped
              }
              else {
                current_stack.emplace (n1idx, c1, n2idx, c2, layidx);
                current_stack.emplace (n1_children[c1], 0, n2_children[c2], 0, layidx + 1);
              }
            }
          }
        }
        node_pair = std::make_pair (n1idx, n2idx);
        cached = simulating[layidx].find (node_pair);
        assert (cached != simulating[layidx].end ());
        return cached->second;
      }

      void add_son (st_node& node, size_t son_layer, size_t son) {
        const int last = node.numchild - 1;
        [[maybe_unused]] // only used in assert below
        const st_node& son_node = layers[son_layer][son];
        size_t* children = child_buffer + node.cbuffer_offset;

        // the new node is smaller than all existing ones
        assert (last == -1 or layers[son_layer][children[last]].label > son_node.label);

        // Insert the new value
        children[last + 1] = son;
        node.numchild++;
      }

      void add_son_unordered (st_node& node, size_t son_layer, size_t son) {
        assert (son_layer <= this->dim);  // Can't add son after the last layer
        // Find the insertion point using binary search
        int left = 0;
        int right = node.numchild - 1;
        const st_node& son_node = layers[son_layer][son];
        size_t* children = child_buffer + node.cbuffer_offset;
        while (left <= right) {
          const int mid = left + ((right - left) / 2);
          assert (mid < static_cast<int> (node.numchild));
          assert (mid >= 0);
          assert (children[mid] < layers[son_layer].size ());
          const typename V::value_type mid_val = layers[son_layer][children[mid]].label;

          if (son_node.label == mid_val) {
            const size_t new_son = node_union (son, children[mid], son_layer);
            children = child_buffer + node.cbuffer_offset;
            children[mid] = new_son;
            return;
          }
          if (mid_val < son_node.label)
            right = mid - 1;
          else
            left = mid + 1;
        }

        // NOLINTBEGIN(bugprone-unchecked-optional-access)
        // Shift elements in the child buffer to make room for the new child
        // Overwrite simulated elements!
        std::optional<size_t> to_insert {son};
        size_t next_insertion {static_cast<size_t> (left)};
        const size_t current_children {node.numchild};
        for (size_t i = left; i < current_children; i++) {
          assert (next_insertion <= i);
          if (not simulates (son, children[i], son_layer)) {
            if (next_insertion < i and not to_insert.has_value ()) {
              // The next insertion point is smaller, so we already looked at it
              // We can directly overwrite the value with our current one
              children[next_insertion] = children[i];
            }
            else if (to_insert.has_value ()) {
              // The space we want to insert is occupied - we save the value for later
              size_t temp = children[i];
              children[next_insertion] = to_insert.value ();
              to_insert = temp;
            }
            next_insertion++;
          }
          else {
            if (to_insert.has_value ()) {
              // There might still be a child we have not inserted so we overwrite
              // No other check necessary because next_insertion is always <= i
              children[next_insertion] = to_insert.value ();
              to_insert.reset ();
              next_insertion++;
            }
            node.numchild--;
          }
        }

        // If we had to shift every element, there will be something left to insert
        // We add it to the new end
        if (to_insert.has_value ())
          children[next_insertion] = to_insert.value ();
        // NOLINTEND(bugprone-unchecked-optional-access)

        // Increase count by the new element
        node.numchild++;
      }

      void add_son_if_not_simulated (size_t node, size_t son_layer, st_node& father) {
        size_t* siblings = child_buffer + father.cbuffer_offset;
        for (size_t s = 0; s < father.numchild; s++)
          if (simulates (siblings[s], node, son_layer))
            return;
        add_son (father, son_layer, node);
      }

      size_t add_node (st_node& node, size_t destination_layer) {
        auto existing_node = inverse[destination_layer].find (node);
        if (existing_node == inverse[destination_layer].end ()) {
          const size_t new_id = layers[destination_layer].size ();
          inverse[destination_layer][node] = new_id;
          layers[destination_layer].push_back (node);
          return new_id;
        }
        return existing_node->second;
      }

      size_t add_children (int num_child) {
        // Double the child buffer if it is full
        if (cbuffer_nxt + num_child >= cbuffer_size) {
          auto* new_buffer = new size_t[cbuffer_size * 2];
          std::memcpy (new_buffer, child_buffer, cbuffer_size * sizeof (size_t));
          cbuffer_size *= 2;
          delete[] child_buffer;
          child_buffer = new_buffer;
        }
        const size_t res = cbuffer_nxt;
        cbuffer_nxt += num_child;
        return res;
      }

      bool is_simulated (size_t nodeidx, st_node& father, size_t destination_layer) {
        size_t* siblings = child_buffer + father.cbuffer_offset;
        for (size_t s = 0; s < father.numchild; s++)
          if (simulates (siblings[s], nodeidx, destination_layer))
            return true;
        return false;
      }

      std::pair<size_t, bool> add_if_not_simulated (st_node& node, size_t destination_layer,
                                                    st_node& father) {
        const size_t res = add_node (node, destination_layer);
        return std::make_pair (res, is_simulated (res, father, destination_layer));
      }

      size_t node_union (size_t ns, size_t nt, size_t destination_layer) {
        // Stack contains node and child ID of S, node and child ID of T, layer
        std::stack<std::tuple<size_t, size_t, size_t, size_t, size_t>> current_stack;

        [[maybe_unused]] const st_node root_nod_e1 = layers[destination_layer][ns];
        [[maybe_unused]] const st_node root_nod_e2 = layers[destination_layer][nt];
        assert (root_nod_e1.numchild > 0 and root_nod_e2.numchild > 0);

        current_stack.emplace (ns, 0, nt, 0, destination_layer);

        while (not current_stack.empty ()) {
          auto [n_s, c_s, n_t, c_t, layer] = current_stack.top ();
          current_stack.pop ();
          assert (n_s < layers[layer].size ());
          const st_node node_s = layers[layer][n_s];
          assert (n_t < layers[layer].size ());
          const st_node node_t = layers[layer][n_t];

          // It's the first time we see this combination of nodes, so let's insert a draft
          // of it, still dirty, will be cleaned later
          if (c_s == 0 and c_t == 0) {
            assert (node_s.label == node_t.label);
            // Before creating a draft node and continuing, let's check the
            // cache
            auto cache_res = cached_union[layer].find (std::make_pair (n_s, n_t));
            if (layer > destination_layer and cache_res != cached_union[layer].end ()) {
#ifndef NDEBUG
              std::cout << "Avoided node in union = cache hit.\n";
#endif
              auto& father = layers[layer - 1].back ();
              if (not is_simulated (cache_res->second, father, layer))
                add_son (father, layer, cache_res->second);
              continue;
            }
            // Not found, so draft a node up
            if (layer < this->dim) {
              layers[layer].emplace_back (node_s.label, 0,
                                          add_children (node_s.numchild + node_t.numchild));
            }
            else {
              layers[layer].emplace_back (node_s.label, 0);
            }
          }

          auto& under_construction = layers[layer].back ();
          // Base case: We are done with this combo and can clean up
          if (c_s == node_s.numchild and c_t == node_t.numchild) {
            if (layer == destination_layer) {
              // We skip this step, we clean the root after the loop, but we are done
              break;
            }
            // The very first thing to do is to check whether our draft of node is
            // not a repetition of something in the table already!
            auto& father = layers[layer - 1].back ();
            layers[layer].pop_back ();
            auto [union_res, domd] = add_if_not_simulated (under_construction, layer, father);
            cached_union[layer][std::make_pair (n_s, n_t)] = union_res;
            if (not domd)
              add_son (father, layer, union_res);
            // Recursive step: Either just add the son to the draft node and
            // continue with the next node, or continue down the tree if
            // necessary
          }
          else if (layer < this->dim) {
            size_t* node_s_children = child_buffer + node_s.cbuffer_offset;
            size_t* node_t_children = child_buffer + node_t.cbuffer_offset;

            // Case 1: One node is "done"
            // We add the existing node (including all its sons!) to the node
            // currently being constructed
            if (c_s == node_s.numchild) {
              add_son_if_not_simulated (node_t_children[c_t], layer + 1, under_construction);
              if (c_t < node_t.numchild)
                current_stack.emplace (n_s, c_s, n_t, c_t + 1, layer);
            }
            else if (c_t == node_t.numchild) {
              add_son_if_not_simulated (node_s_children[c_s], layer + 1, under_construction);
              if (c_s < node_s.numchild)
                current_stack.emplace (n_s, c_s + 1, n_t, c_t, layer);
            }
            else {
              const st_node& son_s = layers[layer + 1][node_s_children[c_s]];
              const st_node& son_t = layers[layer + 1][node_t_children[c_t]];
              // Case 2: One child is "ahead", We add just as when one list is done
              if (son_t.label > son_s.label) {
                add_son_if_not_simulated (node_t_children[c_t], layer + 1, under_construction);
                if (c_t < node_t.numchild)
                  current_stack.emplace (n_s, c_s, n_t, c_t + 1, layer);
              }
              else if (son_s.label > son_t.label) {
                add_son_if_not_simulated (node_s_children[c_s], layer + 1, under_construction);
                if (c_s < node_s.numchild)
                  current_stack.emplace (n_s, c_s + 1, n_t, c_t, layer);
              }
              // Case 3: The values of the sons match, so we need to continue the recursion
              else {
                assert (son_s.label == son_t.label);
                if (c_s < node_s.numchild and c_t < node_t.numchild)
                  current_stack.emplace (n_s, c_s + 1, n_t, c_t + 1, layer);
                current_stack.emplace (node_s_children[c_s], 0, node_t_children[c_t], 0,
                                       layer + 1);
              }
            }
          }
        }

        // Clean up the root, we remove it and re-add it, so it is checked whether an identical
        // node exists
        auto constructed_node = layers[destination_layer].back ();
        layers[destination_layer].pop_back ();
        return add_node (constructed_node, destination_layer);
      }

      std::optional<size_t> node_intersect (size_t n_s, size_t n_t, size_t destination_layer,
                                            std::optional<st_node> father) {
        st_node& node_s = layers[destination_layer][n_s];
        st_node& node_t = layers[destination_layer][n_t];
        st_node new_node {std::min (node_s.label, node_t.label), 0};
        // We haven't reached the bottom of the tree, we need to add children
        if (destination_layer < this->dim) {
          new_node.cbuffer_offset = add_children (node_s.numchild + node_t.numchild);

          size_t* node_s_children = child_buffer + node_s.cbuffer_offset;
          size_t* node_t_children = child_buffer + node_t.cbuffer_offset;
          for (size_t s_s = 0; s_s < node_s.numchild; s_s++) {
            for (size_t s_t = 0; s_t < node_t.numchild; s_t++) {
              auto intersect_res = node_intersect (node_s_children[s_s], node_t_children[s_t],
                                                   destination_layer + 1, new_node);
              if (intersect_res.has_value ()) {
                // Can happen that we insert the same value twice, so we need to check
                auto existing_son =
                    has_son (new_node, destination_layer + 1,
                             layers[destination_layer + 1][intersect_res.value ()].label);
                if (existing_son.has_value ()) {
                  size_t new_son = node_union (existing_son.value (), intersect_res.value (),
                                               destination_layer + 1);
                  size_t* new_node_children = child_buffer + new_node.cbuffer_offset;
                  new_node_children[existing_son.value ()] = new_son;
                }
                else {
                  add_son (new_node, destination_layer + 1, intersect_res.value ());
                }
              }
            }
          }

          if (new_node.numchild == 0)
            return std::nullopt;
        }
        if (father.has_value ()) {
          auto [res, domd] = add_if_not_simulated (new_node, destination_layer, father.value ());
          if (domd)
            return std::nullopt;
          return res;
        }
        return add_node (new_node, destination_layer);
      }

      // NOLINTBEGIN(misc-no-recursion)
      /* Recursive creation of nodes of Trie while using the inverse map to avoid
       * creating duplicate nodes in terms of (residual/right) language. This
       * results on the creation of the minimal DFA for the set of vectors.
       *
       * Complexity: The implementation below has complexity O(n.d.lg(n)) where n is
       * the number of vectors, d is the number of dimensions. The nd factor is
       * nor surprising: it already corresponds to the set of prefixes of the set
       * of vectors. The lg(n) factor comes from the use of an (ordered) map at
       * every recursive step to (re)partition the vectors for the children based
       * on the list itself. It can be traded off by a factor of k (see TODO in
       * code).
       */
      size_t build_node (std::vector<size_t>& vecs, size_t current_layer, const auto& element_vec,
                         bool check_sim = true) {
        assert (not vecs.empty ());
        // If currentLayer is 0, we set the label to the dummy value -1 for the root
        // Else all nodes should have the same value at index currentLayer - 1, so
        // we just use the first
        const typename V::value_type label {current_layer == 0
                                                ? static_cast<typename V::value_type> (-1)
                                                : element_vec[vecs[0]][current_layer - 1]};
        st_node new_node {label, 0};
        // We have not reached the last layer - so add children
        if (current_layer < this->dim) {
          if (vecs.size () == 1) {
            auto& vec = element_vec[vecs[0]];
            st_node last_son_node {vec[vec.size () - 1], 0};
            size_t next_son = add_node (last_son_node, this->dim);
            for (size_t i = vec.size () - 1; i > current_layer; i--) {
              st_node son_node {vec[i - 1], 0};
              son_node.cbuffer_offset = add_children (1);
              add_son (son_node, i + 1, next_son);
              next_son = add_node (son_node, i);
            }
            new_node.cbuffer_offset = add_children (1);
            add_son (new_node, current_layer + 1, next_son);
          }
          else {
            // Partition and order the future children
            // TODO: Try to do constant-access bucketing based on the value of k.
            // Probably won't pay off unless the set of vectors we are adding is
            // dense in most components.
            std::map<typename V::value_type, std::vector<size_t>, std::greater<>> new_partition {};
            for (const auto& vec : vecs)
              new_partition[element_vec[vec][current_layer]].push_back (vec);
            new_node.cbuffer_offset = add_children (new_partition.size ());

            for (auto& [n, children] : new_partition) {
              // Build a new son for each individual value at currentLayer + 1
              const size_t new_son =
                  build_node (children, current_layer + 1, element_vec, check_sim);
              bool found = false;
              if (check_sim) {
                size_t* current_children = child_buffer + new_node.cbuffer_offset;
                for (size_t s = 0; s < new_node.numchild; s++) {
                  if (simulates (current_children[s], new_son, current_layer + 1)) {
                    found = true;
                    break;
                  }
                }
              }
              if (not found)
                add_son (new_node, current_layer + 1, new_son);
            }
          }
        }
        return add_node (new_node, current_layer);
      }
      // NOLINTEND(misc-no-recursion)

    public:
      sharingforest () = delete;

      sharingforest (size_t dim) { this->init (dim); }

      ~sharingforest () {
        if (layers.empty ())
          return;

        delete[] child_buffer;
      }

      [[nodiscard]] std::vector<V> get_all (std::optional<size_t> root = {}) const {
        // Stack with tuples (layer, node id, child id)
        std::stack<std::tuple<size_t, size_t, size_t>> to_visit;

        if (root) {
          const st_node& root_node = layers[0][root.value ()];
          size_t* first_children = child_buffer + root_node.cbuffer_offset;
          for (size_t c = 0; c < root_node.numchild; c++)
            to_visit.emplace (1, first_children[c], 0);
        }
        else {
          // Add all roots at dimension 0
          for (size_t i = 0; i < layers[1].size (); i++)
            to_visit.emplace (1, i, 0);
        }

        std::vector<V> res;
        std::vector<typename V::value_type> temp;
        while (not to_visit.empty ()) {
          const auto [lay, node, child] = to_visit.top ();
          to_visit.pop ();
          const auto parent = layers[lay][node];
          // first time we see parent? get its label
          if (child == 0)
            temp.push_back (parent.label);

          // base case: reached the bottom layer
          if (lay == this->dim) {
            assert (child == 0);
            std::vector<typename V::value_type> cpy {temp};
            res.push_back (V (std::move (cpy)));
            temp.pop_back ();  // done with this node
          }
          else {  // recursive case
            // Either we're done with this node and we just mark it as visited or
            // we need to keep it and we add it's next son
            size_t* children = child_buffer + parent.cbuffer_offset;
            if (child < parent.numchild) {
              to_visit.emplace (lay, node, child + 1);
              to_visit.emplace (lay + 1, children[child], 0);
            }
            else {
              temp.pop_back ();  // done with this parent
            }
          }
        }
        return res;
      }

      void print_children (size_t n, size_t layer) {
#ifndef NDEBUG
        assert (layer <= this->dim);
        st_node& node = layers[layer][n];
        size_t* children = child_buffer + node.cbuffer_offset;
        std::cout << std::string (layer, '\t') << layer << "." << n << " ["
                  << static_cast<int> (node.label) << "] -> (" << (layer == this->dim ? "" : "\n");
        for (size_t i = 0; i < node.numchild; i++)
          print_children (children[i], layer + 1);
        std::cout << (layer == this->dim ? " " : std::string (layer, '\t')) << " )\n";
#endif
      }

      size_t st_union (size_t root1, size_t root2) { return node_union (root1, root2, 0); }

      size_t st_intersect (size_t root1, size_t root2) {
        auto cache_res = cached_inter[0].find (std::make_pair (root1, root2));
        if (cache_res != cached_inter[0].end ())
          return cache_res->second;

        // Stack contains node and child ID of S, node and child ID of T, layer
        std::stack<std::tuple<size_t, size_t, size_t, size_t, size_t>> current_stack;

        // Insert a draft of root node, dirty one without checking if it's there,
        // we'll clean later
        const st_node root_nod_e1 = layers[0][root1];
        const st_node root_nod_e2 = layers[0][root2];
        layers[0].emplace_back (static_cast<typename V::value_type> (-1), 0,
                                add_children (root_nod_e1.numchild + root_nod_e2.numchild));

        // We are ready to start a stack-simulated DFS of the synchronized-product
        // of the trees
        assert (root_nod_e1.numchild > 0 and root_nod_e2.numchild > 0);
        size_t* root1_children = child_buffer + root_nod_e1.cbuffer_offset;
        size_t* root2_children = child_buffer + root_nod_e2.cbuffer_offset;
        for (size_t c_1 = 1; c_1 <= root_nod_e1.numchild; c_1++) {
          for (size_t c_2 = 1; c_2 <= root_nod_e2.numchild; c_2++) {
            current_stack.push ({root1_children[root_nod_e1.numchild - c_1], 0,
                                 root2_children[root_nod_e2.numchild - c_2], 0, 1});
          }
        }

        while (not current_stack.empty ()) {
          auto [n_s, c_s, n_t, c_t, layer] = current_stack.top ();
          current_stack.pop ();
          assert (n_s < layers[layer].size ());
          const st_node node_s = layers[layer][n_s];
          assert (n_t < layers[layer].size ());
          const st_node node_t = layers[layer][n_t];

          // It's the first time we see this product node, so let's insert a draft
          // of it, again dirty
          if (c_s == 0 and c_t == 0) {
            // Before creating a draft node and continuing, let's check the
            // cache
            auto cache_res = cached_inter[layer].find (std::make_pair (n_s, n_t));
            if (cache_res != cached_inter[layer].end ()) {
#ifndef NDEBUG
              std::cout << "Avoided node in intersection = cache hit.\n";
#endif
              auto& father = layers[layer - 1].back ();
              if (not is_simulated (cache_res->second, father, layer))
                add_son_unordered (father, layer, cache_res->second);
              continue;
            }
            // Not found, so draft a node up
            if (layer < this->dim) {
              layers[layer].emplace_back (std::min (node_s.label, node_t.label), 0,
                                          add_children (node_s.numchild + node_t.numchild));
            }
            else {
              layers[layer].emplace_back (std::min (node_s.label, node_t.label), 0);
            }
          }

          // This is our base case, we've iterated through all the children in the
          // product tree, it's finally time to clean up stuff
          if (c_s == node_s.numchild and c_t == node_t.numchild) {
            // The very first thing to do is to check whether our draft of node is
            // not a repetition of something in the table already!
            auto under_construction = layers[layer].back ();
            auto& father = layers[layer - 1].back ();
            layers[layer].pop_back ();
            auto [intersect_res, domd] = add_if_not_simulated (under_construction, layer, father);
            cached_inter[layer][std::make_pair (n_s, n_t)] = intersect_res;
            if (not domd)
              add_son_unordered (father, layer, intersect_res);
            // Below we have the "recursive" step in which we put two elements into
            // the stack to remember what was the next child of the node at this
            // level and to go deeper in the tree if needed
          }
          else if (layer < this->dim) {
            size_t* node_s_children = child_buffer + node_s.cbuffer_offset;
            size_t* node_t_children = child_buffer + node_t.cbuffer_offset;

            if (c_t < node_t.numchild)
              current_stack.emplace (n_s, c_s, n_t, c_t + 1, layer);
            else if (c_s < node_s.numchild)
              current_stack.emplace (n_s, c_s + 1, n_t, 0, layer);
            if (c_t < node_t.numchild and c_s < node_s.numchild)
              current_stack.emplace (node_s_children[c_s], 0, node_t_children[c_t], 0, layer + 1);
          }
        }

        // Clean up the root too! It's added to the layer, so we just need to make
        // sure it's not there already and remove this second copy if it is
        auto under_construction = layers[0].back ();
        layers[0].pop_back ();
        return add_node (under_construction, 0);
      }

      /* Recursive domination check of given vector by vectors in the language of
       * the given tree root.
       *
       * Complexity: The implementation below has complexity O(n.d) where n is
       * the number of vectors, d is the number of dimensions. This is the same as
       * list-based data structures for vectors. However, since it is a DFS of the
       * DFA representation of (bisimulation non-dominated) vectors, it could be
       * exponentially faster than this.
       */
      bool covers_vector (size_t root, const V& covered, bool strict = false) {
        // Stack with tuples (layer, node id, strictness, child id)
        std::stack<std::tuple<size_t, size_t, bool, size_t>> to_visit;
        // Visited cache
        std::vector<std::unordered_map<size_t, bool>> visited;
        visited.reserve (this->dim + 1);
        for (size_t i = 0; i < this->dim + 1; i++)
          visited.emplace_back ();

        // Add all roots at dimension 0 such that their labels cover the first
        // component of the given vector
        const st_node& root_node = layers[0][root];
        size_t* root_children = child_buffer + root_node.cbuffer_offset;
        for (size_t i = 0; i < root_node.numchild; i++) {
          assert (root_children[i] < layers[1].size ());
          if (covered[0] <= layers[1][root_children[i]].label) {
            const bool owe_strict = strict and covered[0] == layers[1][root_children[i]].label;
            to_visit.emplace (1, root_children[i], owe_strict, 0);
          }
        }

        while (not to_visit.empty ()) {
          const auto [lay, node, owe_strict, child] = to_visit.top ();
          assert (lay <= this->dim);
          assert (node < layers[lay].size ());
          to_visit.pop ();
          // if this node has been visited, it means this node is useless
          // NOTE: this works only because we are doing a DFS and not a BFS
          if (child == 0) {
            auto res = visited[lay].find (node);
            // if we visited with strictness then we can safely exit only if
            // we come back with strictness prev_strict -> owe_strict, i.e.
            // not prev_strict or owe_strict
            if (res != visited[lay].end ()) {
              if ((not res->second) or owe_strict) {
#ifndef NDEBUG
                std::cout << "Avoided node in covers check, DFS cache helps.\n";
#endif
                continue;
              }
              // we conjoin with previous result if any to make sure we have
              // more early exits based on implication condition above
              visited[lay][node] = res->second and owe_strict;
            }
            else {
              // no early exit? then mark the node as visited and keep going
              visited[lay][node] = owe_strict;
            }
          }
          const auto parent = layers[lay][node];

          // base case: reached the bottom layer
          if (lay == this->dim) {
            assert (child == 0);
            if (covered[lay - 1] < parent.label or
                ((not owe_strict) and covered[lay - 1] == parent.label))
              return true;
          }
          else {  // recursive case
            // Either we're done with this node and we just mark it as visited or
            // we need to keep it and we add it's next son
            assert (child < parent.numchild);
            size_t* children = child_buffer + parent.cbuffer_offset;
            const size_t c = child;
            auto child_node = layers[lay + 1][children[c]];
            // early exit if the largest child is smaller
            if (covered[lay] > child_node.label)
              continue;

            const bool still_owe_strict = owe_strict and covered[lay] == child_node.label;
            assert (c < parent.numchild);
            if (c + 1 < parent.numchild)
              to_visit.emplace (lay, node, owe_strict, c + 1);
            to_visit.emplace (lay + 1, children[c], still_owe_strict, 0);
          }
        }
        return false;
      }

      template <std::ranges::input_range R>
      size_t add_vectors (R&& elements, bool check_sim = true) {
        assert (not layers.empty ());

        auto element_vec = std::forward<R> (elements);
        // We start a Trie encoded as a map from prefixes to sets of indices of
        // the original vectors, we insert the root too: an empty prefix mapped to
        // the set of all indices
        std::vector<size_t> vector_ids (element_vec.size ());
        // NOLINTBEGIN(boost-use-ranges)
        std::iota (vector_ids.begin (), vector_ids.end (), 0);
        // NOLINTEND(boost-use-ranges)

        const size_t root_id = build_node (vector_ids, 0, element_vec, check_sim);
#ifndef NDEBUG
        size_t maxlayer = 0;
        size_t totlayer = 0;
        for (size_t i = 0; i < this->dim + 1; i++) {
          totlayer += this->layers[i].size ();
          if (maxlayer < this->layers[i].size ())
            maxlayer = this->layers[i].size ();
        }
        std::cout << "[" << this->dim << " Forest stats] Max layer size=" << maxlayer
                  << " total layers' size=" << totlayer << " (bytes per el=" << sizeof (st_node)
                  << ")" << '\n';
        std::cout << "[" << this->dim << " Forest stats] Child buff size=" << this->cbuffer_size
                  << " (bytes per el=" << sizeof (size_t) << ")\n";
#endif
        return root_id;
      }

      /*
        For testing: Check that all children are ordered descending
      */
      bool check_child_order () {
        size_t layer_num = 0;
        for (auto& l : layers) {
          for (st_node& n : l) {
            size_t* children = child_buffer + n.cbuffer_offset;
            for (size_t i = 1; i < n.numchild; i++) {
              // There is a child with a larger label than the previous child
              if (layers[layer_num + 1][children[i]].label >
                  layers[layer_num + 1][children[i - 1]].label) {
                return false;
              }
            }
          }
          layer_num++;
        }
        return true;
      }

      /*
        For testing: Check that one root simulates another
      */
      bool check_simulation (size_t n1, size_t n2) { return simulates (n1, n2, 0); }
  };

  template <Vector V>
  inline std::ostream& operator<< (std::ostream& os, const sharingforest<V>& f) {
    for (auto&& el : f.get_all ())
      os << el << '\n';

    os << "Layers:" << '\n';
    for (auto& l : f.layers) {
      for (auto& n : l)
        os << static_cast<int> (n.label) << " ";
      os << '\n';
    }
    return os;
  }

}  // namespace posets::utils
