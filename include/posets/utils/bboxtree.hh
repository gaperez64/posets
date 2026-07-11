#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>

#include <posets/concepts.hh>

namespace posets::utils {
  // A static bounding-box hierarchy over an antichain.  join_c is an upper
  // corner and meet_c a lower corner of every subtree.  The tree is rebuilt
  // in bulk because Acacia's apply operation changes every stored vector.
  template <Vector V, size_t B = 8>
  class bboxtree {
      static_assert (B >= 2);

      struct node {
          V join_c;
          V meet_c;
          uint32_t begin;
          uint32_t end;
          uint32_t first_child;
          uint32_t nchildren;

          node (V&& join_corner, V&& meet_corner, uint32_t first, uint32_t last, uint32_t child,
                uint32_t children)
            : join_c {std::move (join_corner)},
              meet_c {std::move (meet_corner)},
              begin {first},
              end {last},
              first_child {child},
              nchildren {children} {}
      };

      std::vector<V> vector_set;
      std::vector<node> nodes;
      mutable std::vector<uint32_t> stack_cache;

      void append_node (uint32_t begin, uint32_t end, uint32_t first_child, uint32_t nchildren) {
        V join_corner =
            nchildren == 0 ? vector_set[begin].copy () : nodes[first_child].join_c.copy ();
        V meet_corner =
            nchildren == 0 ? vector_set[begin].copy () : nodes[first_child].meet_c.copy ();
        const uint32_t count = nchildren == 0 ? end - begin : nchildren;
        for (uint32_t i = 1; i < count; ++i) {
          const V& next_join =
              nchildren == 0 ? vector_set[begin + i] : nodes[first_child + i].join_c;
          const V& next_meet =
              nchildren == 0 ? vector_set[begin + i] : nodes[first_child + i].meet_c;
          join_corner = join_corner.join (next_join);
          meet_corner = meet_corner.meet (next_meet);
        }
        nodes.emplace_back (std::move (join_corner), std::move (meet_corner), begin, end,
                            first_child, nchildren);
      }

    public:
      bboxtree () = default;
      bboxtree (const bboxtree&) = delete;
      bboxtree (bboxtree&&) = default;
      bboxtree& operator= (const bboxtree&) = delete;
      bboxtree& operator= (bboxtree&&) = default;

      template <std::ranges::input_range R>
      explicit bboxtree (R&& elements) {
        std::vector<V> values;
        values.reserve (elements.size ());
        for (auto&& value : elements)
          values.push_back (std::move (value));
        assign_and_build (std::move (values));
      }

      void assign_and_build (std::vector<V>&& sorted_antichain) {
        assert (not sorted_antichain.empty ());
        assert (sorted_antichain.size () <= UINT32_MAX);
        vector_set = std::move (sorted_antichain);
        nodes.clear ();
        nodes.reserve ((vector_set.size () * B) / (B - 1) + 1);

        uint32_t level_begin = 0;
        for (uint32_t begin = 0; begin < vector_set.size (); begin += B) {
          const uint32_t end = std::min<uint32_t> (begin + B, vector_set.size ());
          append_node (begin, end, 0, 0);
        }
        uint32_t level_end = nodes.size ();

        while (level_end - level_begin > 1) {
          const uint32_t next_begin = nodes.size ();
          for (uint32_t first = level_begin; first < level_end; first += B) {
            const uint32_t nchildren = std::min<uint32_t> (B, level_end - first);
            append_node (nodes[first].begin, nodes[first + nchildren - 1].end, first, nchildren);
          }
          level_begin = next_begin;
          level_end = nodes.size ();
        }
        stack_cache.reserve (nodes.size ());
      }

      [[nodiscard]] bool dominates (const V& v, bool strict = false) const {
        if (nodes.empty ())
          return false;
        stack_cache.clear ();
        stack_cache.push_back (nodes.size () - 1);

        while (not stack_cache.empty ()) {
          const node& current = nodes[stack_cache.back ()];
          stack_cache.pop_back ();

          auto upper = v.partial_order (current.join_c);
          if (not upper.leq ())
            continue;

          auto lower = v.partial_order (current.meet_c);
          if (lower.leq ()) {
            if (not strict or not lower.geq () or current.end - current.begin >= 2)
              return true;
          }

          if (current.nchildren == 0) {
            for (uint32_t i = current.begin; i < current.end; ++i) {
              auto po = v.partial_order (vector_set[i]);
              if (po.leq () and (not strict or not po.geq ()))
                return true;
            }
          }
          else {
            for (uint32_t i = 0; i < current.nchildren; ++i)
              stack_cache.push_back (current.first_child + i);
          }
        }
        return false;
      }

      [[nodiscard]] size_t size () const { return vector_set.size (); }
      [[nodiscard]] bool empty () const { return vector_set.empty (); }
      [[nodiscard]] auto begin () { return vector_set.begin (); }
      [[nodiscard]] auto begin () const { return vector_set.begin (); }
      [[nodiscard]] auto end () { return vector_set.end (); }
      [[nodiscard]] auto end () const { return vector_set.end (); }
      [[nodiscard]] auto& get_backing_vector () { return vector_set; }
      [[nodiscard]] const auto& get_backing_vector () const { return vector_set; }

      [[nodiscard]] bool is_antichain () const {
        for (auto it = vector_set.begin (); it != vector_set.end (); ++it)
          for (auto jt = it + 1; jt != vector_set.end (); ++jt) {
            auto po = it->partial_order (*jt);
            if (po.leq () or po.geq ())
              return false;
          }
        return true;
      }
  };
}
