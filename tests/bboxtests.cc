#include <cassert>
#include <cstdint>
#include <span>
#include <vector>

#include <posets/utils/bboxtree.hh>
#include <posets/utils/reduce.hh>
#include <posets/vectors.hh>

using V = posets::vectors::simd_vector_backed_sum<std::int8_t>;

V make (std::initializer_list<std::int8_t> value) {
  std::vector<std::int8_t> storage (value);
  return V (std::span<const std::int8_t> (storage));
}

int main () {
  std::vector<V> values;
  values.push_back (make ({6, 3, 2}));
  values.push_back (make ({5, 5, 4}));
  values.push_back (make ({2, 6, 2}));
  values = posets::utils::reduce_to_maxima (std::move (values));

  posets::utils::bboxtree<V, 2> tree;
  tree.assign_and_build (std::move (values));
  assert (tree.size () == 3);
  assert (tree.is_antichain ());
  auto low = make ({2, 2, 2});
  auto exact = make ({6, 3, 2});
  auto high = make ({7, 7, 7});
  assert (tree.dominates (low));
  assert (tree.dominates (low, true));
  assert (tree.dominates (exact));
  assert (not tree.dominates (exact, true));
  assert (not tree.dominates (high));
}
