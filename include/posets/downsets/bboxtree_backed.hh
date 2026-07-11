#pragma once

#include <algorithm>
#include <cassert>
#include <iterator>
#include <ostream>
#include <vector>

#include <posets/concepts.hh>
#include <posets/utils/bboxtree.hh>
#include <posets/utils/reduce.hh>

namespace posets::downsets {
  template <Vector V>
  class bboxtree_backed {
      utils::bboxtree<V> tree;

      void reset_tree (std::vector<V>&& elements) {
        auto antichain = utils::reduce_to_maxima (std::move (elements));
        assert (not antichain.empty ());
        tree.assign_and_build (std::move (antichain));
        assert (tree.is_antichain ());
      }

    public:
      using value_type = V;

      bboxtree_backed () = delete;
      explicit bboxtree_backed (std::vector<V>&& elements) { reset_tree (std::move (elements)); }
      explicit bboxtree_backed (V&& element) {
        std::vector<V> singleton;
        singleton.push_back (std::move (element));
        tree.assign_and_build (std::move (singleton));
      }
      bboxtree_backed (const bboxtree_backed&) = delete;
      bboxtree_backed (bboxtree_backed&&) = default;
      bboxtree_backed& operator= (const bboxtree_backed&) = delete;
      bboxtree_backed& operator= (bboxtree_backed&&) = default;

      [[nodiscard]] bool contains (const V& value) const { return tree.dominates (value); }
      [[nodiscard]] size_t size () const { return tree.size (); }

      void union_with (bboxtree_backed&& other) {
        std::vector<bool> keep_this (size ());
        std::vector<bool> keep_other (other.size ());
        for (size_t i = 0; i < size (); ++i)
          keep_this[i] = not other.tree.dominates (tree.get_backing_vector ()[i], true);
        for (size_t i = 0; i < other.size (); ++i)
          keep_other[i] = not tree.dominates (other.tree.get_backing_vector ()[i]);

        std::vector<V> lhs;
        std::vector<V> rhs;
        lhs.reserve (size ());
        rhs.reserve (other.size ());
        for (size_t i = 0; i < keep_this.size (); ++i)
          if (keep_this[i])
            lhs.push_back (std::move (tree.get_backing_vector ()[i]));
        for (size_t i = 0; i < keep_other.size (); ++i)
          if (keep_other[i])
            rhs.push_back (std::move (other.tree.get_backing_vector ()[i]));

        std::vector<V> merged;
        merged.reserve (lhs.size () + rhs.size ());
        std::merge (std::make_move_iterator (lhs.begin ()), std::make_move_iterator (lhs.end ()),
                    std::make_move_iterator (rhs.begin ()), std::make_move_iterator (rhs.end ()),
                    std::back_inserter (merged), [] (const V& lhs, const V& rhs) {
                      return utils::sum_key (lhs) > utils::sum_key (rhs);
                    });
        assert (not merged.empty ());
        tree.assign_and_build (std::move (merged));
      }

      void intersect_with (const bboxtree_backed& other) {
        std::vector<V> intersection;
        bool smaller_set = false;
        for (const auto& x : tree) {
          const bool dominated = other.tree.dominates (x);
          if (dominated)
            intersection.push_back (x.copy ());
          else
            for (const auto& y : other)
              intersection.push_back (x.meet (y));
          smaller_set or_eq not dominated;
        }
        if (smaller_set)
          reset_tree (std::move (intersection));
      }

      template <typename F>
      bboxtree_backed apply (const F& lambda) const {
        std::vector<V> result;
        result.reserve (size ());
        for (const auto& value : tree)
          result.push_back (lambda (value));
        return bboxtree_backed (std::move (result));
      }

      [[nodiscard]] auto begin () { return tree.begin (); }
      [[nodiscard]] auto begin () const { return tree.begin (); }
      [[nodiscard]] auto end () { return tree.end (); }
      [[nodiscard]] auto end () const { return tree.end (); }
      [[nodiscard]] auto& get_backing_vector () { return tree.get_backing_vector (); }
      [[nodiscard]] const auto& get_backing_vector () const { return tree.get_backing_vector (); }
  };

  template <Vector V>
  inline std::ostream& operator<< (std::ostream& os, const bboxtree_backed<V>& downset) {
    for (const auto& value : downset)
      os << value << '\n';
    return os;
  }
}
