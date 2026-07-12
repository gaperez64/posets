#pragma once

#include <algorithm>
#include <cassert>
#include <iterator>
#include <ostream>
#include <vector>

#include <posets/concepts.hh>
#include <posets/utils/comptrie.hh>
#include <posets/utils/reduce.hh>

namespace posets::downsets {
  template <Vector V>
  class comptrie_backed {
      utils::comptrie<V> trie;

      void reset_trie (std::vector<V>&& elements) {
        auto antichain = utils::reduce_to_maxima (std::move (elements));
        assert (not antichain.empty ());
        trie.assign_and_build (std::move (antichain));
        assert (trie.is_antichain ());
      }

    public:
      using value_type = V;
      comptrie_backed () = delete;
      explicit comptrie_backed (std::vector<V>&& elements) { reset_trie (std::move (elements)); }
      explicit comptrie_backed (V&& element) {
        std::vector<V> singleton;
        singleton.push_back (std::move (element));
        trie.assign_and_build (std::move (singleton));
      }
      comptrie_backed (const comptrie_backed&) = delete;
      comptrie_backed (comptrie_backed&&) = default;
      comptrie_backed& operator= (const comptrie_backed&) = delete;
      comptrie_backed& operator= (comptrie_backed&&) = default;

      [[nodiscard]] bool contains (const V& value) const { return trie.dominates (value); }
      [[nodiscard]] size_t size () const { return trie.size (); }

      void union_with (comptrie_backed&& other) {
        std::vector<bool> keep_this (size ());
        std::vector<bool> keep_other (other.size ());
        for (size_t i = 0; i < size (); ++i)
          keep_this[i] = not other.trie.dominates (trie.get_backing_vector ()[i], true);
        for (size_t i = 0; i < other.size (); ++i)
          keep_other[i] = not trie.dominates (other.trie.get_backing_vector ()[i]);

        std::vector<V> lhs;
        std::vector<V> rhs;
        lhs.reserve (size ());
        rhs.reserve (other.size ());
        for (size_t i = 0; i < keep_this.size (); ++i)
          if (keep_this[i])
            lhs.push_back (std::move (trie.get_backing_vector ()[i]));
        for (size_t i = 0; i < keep_other.size (); ++i)
          if (keep_other[i])
            rhs.push_back (std::move (other.trie.get_backing_vector ()[i]));
        std::vector<V> merged;
        merged.reserve (lhs.size () + rhs.size ());
        std::merge (std::make_move_iterator (lhs.begin ()), std::make_move_iterator (lhs.end ()),
                    std::make_move_iterator (rhs.begin ()), std::make_move_iterator (rhs.end ()),
                    std::back_inserter (merged), [] (const V& lhs, const V& rhs) {
                      return utils::sum_key (lhs) > utils::sum_key (rhs);
                    });
        trie.assign_and_build (std::move (merged));
      }

      void intersect_with (const comptrie_backed& other) {
        std::vector<V> intersection;
        bool smaller_set = false;
        for (const auto& x : trie) {
          const bool dominated = other.trie.dominates (x);
          if (dominated)
            intersection.push_back (x.copy ());
          else
            for (const auto& y : other)
              intersection.push_back (x.meet (y));
          smaller_set or_eq not dominated;
        }
        if (smaller_set)
          reset_trie (std::move (intersection));
      }

      template <typename F>
      comptrie_backed apply (const F& lambda) const {
        std::vector<V> result;
        result.reserve (size ());
        for (const auto& value : trie)
          result.push_back (lambda (value));
        return comptrie_backed (std::move (result));
      }

      [[nodiscard]] auto begin () { return trie.begin (); }
      [[nodiscard]] auto begin () const { return trie.begin (); }
      [[nodiscard]] auto end () { return trie.end (); }
      [[nodiscard]] auto end () const { return trie.end (); }
      [[nodiscard]] auto& get_backing_vector () { return trie.get_backing_vector (); }
      [[nodiscard]] const auto& get_backing_vector () const { return trie.get_backing_vector (); }
  };

  template <Vector V>
  inline std::ostream& operator<< (std::ostream& os, const comptrie_backed<V>& downset) {
    for (const auto& value : downset)
      os << value << '\n';
    return os;
  }
}
