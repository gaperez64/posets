#pragma once

#include <algorithm>
#include <concepts>
#include <numeric>
#include <utility>
#include <vector>

#include <posets/concepts.hh>
#include <posets/utils/kdtree.hh>

#ifndef POSETS_REDUCE_SCAN_THRESHOLD
# define POSETS_REDUCE_SCAN_THRESHOLD 4096UL
#endif

namespace posets::utils {
  template <Vector V>
  [[nodiscard]] long sum_key (const V& v) {
    if constexpr (requires {
                    { v.cached_sum () } -> std::convertible_to<long>;
                  })
      return static_cast<long> (v.cached_sum ());
    else {
      long sum = 0;
      for (size_t i = 0; i < v.size (); ++i)
        sum += static_cast<long> (v[i]);
      return sum;
    }
  }

  // Sort move-only vectors without repeatedly recomputing uncached sums from
  // all their coordinates in the comparison function.
  template <Vector V>
  void sort_by_sum_desc (std::vector<V>& values, std::vector<long>* sorted_keys = nullptr) {
    std::vector<long> keys;
    keys.reserve (values.size ());
    for (const auto& value : values)
      keys.push_back (sum_key (value));

    std::vector<size_t> order (values.size ());
    std::iota (order.begin (), order.end (), 0);  // NOLINT(boost-use-ranges)
    std::stable_sort (order.begin (), order.end (),
                      [&keys] (size_t lhs, size_t rhs) { return keys[lhs] > keys[rhs]; });

    std::vector<V> sorted;
    sorted.reserve (values.size ());
    std::vector<long> reordered_keys;
    if (sorted_keys)
      reordered_keys.reserve (values.size ());
    for (const size_t i : order) {
      if (sorted_keys)
        reordered_keys.push_back (keys[i]);
      sorted.push_back (std::move (values[i]));
    }
    values = std::move (sorted);
    if (sorted_keys)
      *sorted_keys = std::move (reordered_keys);
  }

  // Componentwise domination is monotone in the coordinate sum.  Once the
  // candidates are in descending-sum order, a new candidate can only be
  // dominated by an element already kept.  Equal sums make non-strict
  // domination equivalent to equality, so this also removes duplicates.
  template <Vector V>
  [[nodiscard]] std::vector<V> reduce_to_maxima (std::vector<V>&& candidates,
                                                 std::vector<long>* maxima_keys = nullptr) {
    if (candidates.size () < 2) {
      if (maxima_keys) {
        maxima_keys->clear ();
        if (not candidates.empty ())
          maxima_keys->push_back (sum_key (candidates.front ()));
      }
      return std::move (candidates);
    }

    std::vector<long> keys;
    sort_by_sum_desc (candidates, &keys);

    if (candidates.size () > POSETS_REDUCE_SCAN_THRESHOLD) {
      // Avoid a quadratic flat scan for the n*m meet clouds produced by
      // intersection.  Equal-sum domination means equality, so first make
      // equal-sum runs deterministic and remove adjacent duplicates.
      size_t begin = 0;
      while (begin < candidates.size ()) {
        size_t end = begin + 1;
        while (end < candidates.size () and keys[end] == keys[begin])
          ++end;
        std::sort (candidates.begin () + begin, candidates.begin () + end);
        begin = end;
      }
      std::vector<V> unique_candidates;
      std::vector<long> unique_keys;
      unique_candidates.reserve (candidates.size ());
      unique_keys.reserve (candidates.size ());
      for (size_t i = 0; i < candidates.size (); ++i)
        if (unique_candidates.empty () or candidates[i] != unique_candidates.back ()) {
          unique_candidates.push_back (std::move (candidates[i]));
          unique_keys.push_back (keys[i]);
        }

      utils::kdtree<V> filter_tree;
      filter_tree.relabel_tree (std::move (unique_candidates));
      std::ranges::reverse (unique_keys);
      std::vector<bool> keep (filter_tree.size ());
      bool removed = false;
      for (size_t i = 0; i < filter_tree.size (); ++i) {
        keep[i] = not filter_tree.dominates (filter_tree.get_backing_vector ()[i], true);
        removed or_eq not keep[i];
      }
      if (not removed) {
        auto maxima = std::move (filter_tree.get_backing_vector ());
        std::ranges::reverse (maxima);
        std::ranges::reverse (unique_keys);
        if (maxima_keys)
          *maxima_keys = std::move (unique_keys);
        return maxima;
      }

      std::vector<V> maxima;
      std::vector<long> kept_keys;
      maxima.reserve (filter_tree.size ());
      kept_keys.reserve (filter_tree.size ());
      for (size_t i = 0; i < keep.size (); ++i)
        if (keep[i]) {
          maxima.push_back (std::move (filter_tree.get_backing_vector ()[i]));
          kept_keys.push_back (unique_keys[i]);
        }
      std::ranges::reverse (maxima);
      std::ranges::reverse (kept_keys);
      if (maxima_keys)
        *maxima_keys = std::move (kept_keys);
      return maxima;
    }

    std::vector<V> maxima;
    std::vector<long> kept_keys;
    maxima.reserve (candidates.size ());
    if (maxima_keys)
      kept_keys.reserve (candidates.size ());
    for (size_t i = 0; i < candidates.size (); ++i) {
      auto& candidate = candidates[i];
      const bool dominated = std::ranges::any_of (
          maxima, [&candidate] (const V& kept) { return candidate.partial_order (kept).leq (); });
      if (not dominated) {
        maxima.push_back (std::move (candidate));
        if (maxima_keys)
          kept_keys.push_back (keys[i]);
      }
    }
    if (maxima_keys)
      *maxima_keys = std::move (kept_keys);
    return maxima;
  }
}
