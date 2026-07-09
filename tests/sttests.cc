#include <cstdint>
#include <vector>

#include <posets/utils/sharingforest.hh>
#include <posets/vectors.hh>

namespace utils = posets::utils;

using test_value_type = std::int8_t;
using VType = posets::vectors::vector_backed<test_value_type>;

std::vector<VType> vvtovv (const std::vector<std::vector<test_value_type>>& vv) {
  std::vector<VType> out;
  for (size_t i = 0; i < vv.size (); ++i)
    out.emplace_back (VType (std::move (vv[i])));
  return out;
}

int main (int argc, const char* argv[]) {
  utils::sharingforest<VType> f {3};
  utils::sharingforest<VType> f4 {4};

  // Add some vectors to create a tree
  std::vector<std::vector<test_value_type>> data {{6, 3, 2}, {5, 5, 4}, {2, 6, 2}};
  auto idcs = f.add_vectors (std::move (vvtovv (data)));
  std::cout << f << std::endl;
  f.print_children (idcs, 0);
  std::vector<test_value_type> v = {2, 2, 2};
  assert (f.covers_vector (idcs, VType (std::move (v))));
  v = {6, 3, 1};
  assert (f.covers_vector (idcs, VType (std::move (v))));
  v = {8, 3, 1};
  assert (not f.covers_vector (idcs, VType (std::move (v))));
  v = {9, 4, 8};
  assert (not f.covers_vector (idcs, VType (std::move (v))));

  // We do some longer examples
  utils::sharingforest<VType> f10 {10};
  std::vector<std::vector<test_value_type>> longdata {{2, 4, 1, 8, 7, 4, 1, 10, 2, 8}};
  auto longidcs = f10.add_vectors (std::move (vvtovv (longdata)));
  v = {2, 4, 1, 8, 5, 4, 1, 8, 1, 9};
  assert (not f10.covers_vector (longidcs, VType (std::move (v))));
  v = {2, 4, 1, 6, 7, 4, 1, 10, 0, 4};
  assert (f10.covers_vector (longidcs, VType (std::move (v)), false));
  assert (f10.covers_vector (longidcs, VType (std::move (v)), true));

  // Test adding a second tree to the forest
  data = {{7, 4, 3}, {4, 8, 4}, {2, 5, 6}};
  auto idcs2 = f.add_vectors (std::move (vvtovv (data)));
  std::cout << f << std::endl;
  f.print_children (idcs2, 0);

  // And yet another tree, this time with a shared suffix
  // and with dimension 4
  data = {{3, 2, 2, 1}, {4, 1, 2, 1}, {5, 0, 2, 1}};
  auto fourdim = f4.add_vectors (std::move (vvtovv (data)));
  std::cout << f4 << std::endl;
  f4.print_children (fourdim, 0);
  v = {3, 0, 2, 2};
  assert (not f4.covers_vector (fourdim, VType (std::move (v))));
  v = {1, 0, 1, 1};
  assert (f4.covers_vector (fourdim, VType (std::move (v))));

  // We repeat the first tests above to check we did not break the previous tree
  v = {2, 2, 2};
  assert (f.covers_vector (idcs, VType (std::move (v))));
  v = {6, 3, 1};
  assert (f.covers_vector (idcs, VType (std::move (v))));
  v = {8, 3, 1};
  assert (not f.covers_vector (idcs, VType (std::move (v))));

  // We can now make tests about the new tree
  v = {2, 2, 2};
  assert (f.covers_vector (idcs2, VType (std::move (v))));
  v = {7, 3, 1};
  assert (f.covers_vector (idcs2, VType (std::move (v))));
  v = {8, 3, 1};
  assert (not f.covers_vector (idcs2, VType (std::move (v))));
  v = {8, 5, 3};
  assert (not f.covers_vector (idcs2, VType (std::move (v))));

  // One more test so that covers needs to explore branches with shared
  // prefixes
  data = {{3, 2, 2}, {3, 4, 1}, {3, 2, 3}, {3, 4, 0}};
  auto idcs3 = f.add_vectors (std::move (vvtovv (data)));
  std::cout << f << std::endl;
  v = {3, 2, 3};
  assert (f.covers_vector (idcs3, VType (std::move (v))));
  v = {3, 2, 4};
  assert (not f.covers_vector (idcs3, VType (std::move (v))));
  v = {3, 4, 1};
  assert (f.covers_vector (idcs3, VType (std::move (v))));
  v = {3, 4, 2};
  assert (not f.covers_vector (idcs3, VType (std::move (v))));

  // Test that all children are ordered descending
  assert (f.check_child_order ());

  // Test that the root of a ST is simulated by the
  // root of the original ST
  data = {{3, 4, 1}, {5, 2, 1}, {5, 1, 2}};
  auto idcs4 = f.add_vectors (std::move (vvtovv (data)));
  assert (f.check_simulation (idcs, idcs4));

  // Test that the root of a second ST is not simulated
  // by the root of the original ST
  data = {{5, 2, 2}, {4, 9, 9}};
  auto idcs5 = f.add_vectors (std::move (vvtovv (data)));
  assert (not f.check_simulation (idcs, idcs5));

  // Test that a single vector is dominated by a tree
  data = {{1, 2, 3}, {2, 5, 1}, {4, 1, 1}};
  auto idcs5a = f.add_vectors (std::move (vvtovv (data)));
  data = {{0, 1, 2}};
  auto idcs5b = f.add_vectors (std::move (vvtovv (data)));
  assert (f.check_simulation (idcs5a, idcs5b));

  // Union tests: Check for vectors that are only in one ST
  // and a vector that is in none
  auto uRoot = f.st_union (idcs, idcs2);
  std::cout << "Union result:\n";
  f.print_children (uRoot, 0);

  v = {7, 4, 1};
  assert (f.covers_vector (uRoot, VType (std::move (v))));
  v = {5, 5, 3};
  assert (f.covers_vector (uRoot, VType (std::move (v))));
  v = {8, 3, 1};
  assert (not f.covers_vector (uRoot, VType (std::move (v))));

  // Intersection test: Exact opposite of union
  auto iRoot = f.st_intersect (idcs, idcs2);
  std::cout << "Intersection result:\n";
  f.print_children (iRoot, 0);
  v = {7, 4, 1};
  assert (not f.covers_vector (iRoot, VType (std::move (v))));
  v = {5, 5, 3};
  assert (not f.covers_vector (iRoot, VType (std::move (v))));
  v = {3, 5, 4};
  assert (f.covers_vector (iRoot, VType (std::move (v))));

  return 0;
}
