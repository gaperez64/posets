posets_dep = declare_dependency(include_directories: ['.', boost_inc])

header_files = [
  'posets/downsets/full_set.hh',
  'posets/downsets/kdtree_backed.hh',
  'posets/downsets/set_backed.hh',
  'posets/downsets/sharingtree_backed.hh',
  'posets/downsets/sharingtrie_backed.hh',
  'posets/downsets/simple_sharingtree_backed.hh',
  'posets/downsets/vector_backed.hh',
  'posets/downsets/vector_backed_bin.hh',
  'posets/downsets/vector_backed_one_dim_split_intersection_only.hh',
  'posets/downsets/vector_or_kdtree_backed.hh',
  'posets/downsets.hh',
  'posets/utils/kdtree.hh',
  'posets/utils/sharingforest.hh',
  'posets/utils/sharingtrie.hh',
  'posets/utils/ref_ptr_cmp.hh',
  'posets/utils/simd_traits.hh',
  'posets/utils/vector_mm.hh',
  'posets/vectors/generic.hh',
  'posets/vectors/generic_partial_order.hh',
  'posets/vectors/generic_helpers.hh',
  'posets/vectors/X_and_bitset.hh',
  'posets/vectors.hh'
]

install_headers(header_files, preserve_path: true)
