#pragma once

#include <cassert>
#include <cstring>
#include <experimental/simd>
#include <iostream>

#include <posets/concepts.hh>
#include <posets/utils/simd_traits.hh>
#include <posets/vectors/generic_helpers.hh>
#include <posets/vectors/generic_partial_order.hh>

namespace posets::vectors {
  template <typename Data, bool HasSum, bool EmbedsData>
    requires HasData<Data>
  class generic : private sum_member<HasSum>, private malloc_member<EmbedsData, Data> {
    public:
      using block_type = typename Data::value_type;
      using value_type = block_type::value_type;
      static const auto uses_simd = [] () { return IsDataSimd<Data>; }();

      static const auto items_per_block = sizeof (block_type) / sizeof (value_type);

      static constexpr size_t blocks_for (size_t nelts) {
        return (nelts + items_per_block - 1) / items_per_block;
      }

    private:
      static const bool is_resizable = ([] () { return requires (Data d) { d.resize (42); }; }) ();

      static_assert (not is_resizable or EmbedsData,
                     "Resizable data should be embeded, as they are already managed by pointers.");

      block_type* data () {
        if constexpr (EmbedsData)
          return datap.data ();
        else
          return datap->data ();
      }

      [[nodiscard]] const block_type* data () const {
        if constexpr (EmbedsData)
          return datap.data ();
        else
          return datap->data ();
      }

      [[nodiscard]] size_t data_size () const {
        if constexpr (EmbedsData)
          return datap.size ();
        else
          return datap->size ();
      }

      void clear_back () {
        if (data_size () > blocks_for (k) or k % items_per_block) {
          char* start = reinterpret_cast<char*> (data () + (k / items_per_block));
          char* end = reinterpret_cast<char*> (data () + data_size ());
          memset (start, 0, end - start);
        }
      }

    public:
      generic (size_t k)
        requires is_resizable
        : k {k},
          datap {blocks_for (k)} {
        clear_back ();
      }

      generic (size_t k)
        requires (not is_resizable)
        : k {k} {
        if constexpr (not EmbedsData)
          datap = this->malloc.construct ();
        assert (data_size () >= blocks_for (k));
        clear_back ();
      }

      generic (std::span<const value_type> v) : generic (v.size ()) {
        if constexpr (HasSum) {
          this->sum = 0;
          for (auto&& c : v)
            this->sum += c;
        }
        std::memcpy (reinterpret_cast<value_type*> (data ()), v.data (),
                     v.size () * sizeof (value_type));
      }

      generic () = delete;
      generic (const generic& other) = delete;
      generic (generic&& other) noexcept : k {other.k}, datap {std::move (other.datap)} {
        if constexpr (not EmbedsData)
          other.datap = nullptr;
        if constexpr (HasSum)
          this->sum = other.sum;
      }

      ~generic () {
        if constexpr (not EmbedsData)
          if (datap)
            this->malloc.destroy (datap);
      }

      // explicit copy operator
      [[nodiscard]] generic copy () const {
        if constexpr (EmbedsData) {
          auto res = generic (k);
          res.datap = datap;
          if constexpr (HasSum)
            res.sum = this->sum;
          return res;
        }
        else {
          auto res = generic (std::span ((value_type*) datap, k));
          if constexpr (HasSum)
            res.sum = this->sum;
          return res;
        }
      }

      generic& operator= (generic&& other) noexcept {
        if constexpr (not EmbedsData)
          if (datap)
            this->malloc.destroy (datap);
        k = other.k;
        datap = other.datap;
        if constexpr (not EmbedsData)
          other.datap = nullptr;
        if constexpr (HasSum)
          this->sum = other.sum;

        return *this;
      }

      generic& operator= (const generic& other) = delete;

      void to_vector (std::span<value_type> v) const {
        assert (v.size () >= k);
        std::memcpy (static_cast<void*> (v.data ()), static_cast<const void*> (data ()),
                     k * sizeof (value_type));
      }

      [[nodiscard]] auto partial_order (const generic& rhs) const {
        return generic_partial_order (*this, rhs);
      }

      bool operator== (const generic& rhs) const {
        if constexpr (HasSum)
          if (this->sum != rhs.sum)
            return false;
        // Trust memcmp to DTRT
        return std::memcmp (rhs.data (), data (), k * sizeof (value_type)) == 0;
      }

      bool operator!= (const generic& rhs) const {
        if constexpr (HasSum)
          if (this->sum != rhs.sum)
            return true;
        // Trust memcmp to DTRT
        return std::memcmp (rhs.data (), data (), k * sizeof (value_type)) != 0;
      }

      // Used by Sets, should be a total order.  Do not use.
      bool operator< (const generic& rhs) const {
        for (size_t i = 0; i < data_size (); ++i) {
          if constexpr (uses_simd) {
            auto lhs_lt_rhs = data ()[i] < rhs.data ()[i];
            auto rhs_lt_lhs = rhs.data ()[i] < data ()[i];
            auto p1 = find_first_set (lhs_lt_rhs);
            auto p2 = find_first_set (rhs_lt_lhs);
            if (p1 == p2)
              continue;
            return (p1 < p2);
          }
          else {
            // This is the lexicographical order.
            auto order = data ()[i] <=> rhs.data ()[i];
            if (order == std::weak_ordering::equivalent)
              continue;
            return order == std::weak_ordering::less;
          }
        }
        return false;
      }

      [[nodiscard]] generic meet (const generic& rhs) const {
        auto res = generic (k);
        if constexpr (not EmbedsData)
          res.datap = this->malloc.construct ();

        for (size_t i = 0; i < data_size (); ++i) {
          if constexpr (uses_simd)
            res.data ()[i] = std::experimental::min (data ()[i], rhs.data ()[i]);
          else
            for (size_t j = 0; j < items_per_block; ++j) {
              auto pos = (i * items_per_block) + j;
              res.at (pos) = std::min ((*this)[pos], rhs[pos]);
            }

          // In case of SIMD, this:
          //   res.sum += std::experimental::reduce (res.data ()[i]);
          // should NOT be used since this can lead to overflows over char.
          // instead, we manually loop through:
          if constexpr (HasSum)
            for (size_t j = 0; j < items_per_block; ++j)
              res.sum += res.data ()[i][j];
        }

        return res;
      }

      [[nodiscard]] generic join (const generic& rhs) const {
        auto res = generic (k);
        if constexpr (not EmbedsData)
          res.datap = this->malloc.construct ();

        for (size_t i = 0; i < data_size (); ++i) {
          if constexpr (uses_simd)
            res.data ()[i] = std::experimental::max (data ()[i], rhs.data ()[i]);
          else
            for (size_t j = 0; j < items_per_block; ++j) {
              const auto pos = (i * items_per_block) + j;
              res.at (pos) = std::max ((*this)[pos], rhs[pos]);
            }

          // SIMD reductions over narrow element types can overflow, so keep
          // the scalar accumulation used by meet().
          if constexpr (HasSum)
            for (size_t j = 0; j < items_per_block; ++j)
              res.sum += res.data ()[i][j];
        }

        return res;
      }

      void meet_with (const generic& rhs) {
        assert (k == rhs.k);
        int updated_sum = 0;  // NOLINT(misc-const-correctness): mutated only in HasSum variants.
        for (size_t i = 0; i < data_size (); ++i) {
          if constexpr (uses_simd)
            data ()[i] = std::experimental::min (data ()[i], rhs.data ()[i]);
          else
            for (size_t j = 0; j < items_per_block; ++j) {
              const auto pos = (i * items_per_block) + j;
              at (pos) = std::min ((*this)[pos], rhs[pos]);
            }

          if constexpr (HasSum)
            for (size_t j = 0; j < items_per_block; ++j)
              updated_sum += data ()[i][j];
        }
        if constexpr (HasSum)
          this->sum = updated_sum;
      }

      void join_with (const generic& rhs) {
        assert (k == rhs.k);
        int updated_sum = 0;  // NOLINT(misc-const-correctness): mutated only in HasSum variants.
        for (size_t i = 0; i < data_size (); ++i) {
          if constexpr (uses_simd)
            data ()[i] = std::experimental::max (data ()[i], rhs.data ()[i]);
          else
            for (size_t j = 0; j < items_per_block; ++j) {
              const auto pos = (i * items_per_block) + j;
              at (pos) = std::max ((*this)[pos], rhs[pos]);
            }

          if constexpr (HasSum)
            for (size_t j = 0; j < items_per_block; ++j)
              updated_sum += data ()[i][j];
        }
        if constexpr (HasSum)
          this->sum = updated_sum;
      }

      [[nodiscard]] auto size () const { return k; }

      auto& print (std::ostream& os) const {
        os << "{ ";
        for (size_t i = 0; i < k; ++i)
          os << (int) (*this)[i] << " ";
        os << "}";
        return os;
      }

    private:
      value_type& at (size_t i) { return *(reinterpret_cast<value_type*> (data ()) + i); }

      [[nodiscard]] const value_type& at (size_t i) const {
        return *(reinterpret_cast<const value_type*> (data ()) + i);
      }

    public:
      const value_type& operator[] (size_t i) const { return at (i); }

      using iterator = value_type*;
      using const_iterator = const value_type*;

      iterator begin () { return iterator (data ()); }
      [[nodiscard]] const_iterator begin () const { return const_iterator (data ()); }
      iterator end () { return iterator (data ()) + k; }
      [[nodiscard]] const_iterator end () const { return const_iterator (data ()) + k; }

      [[nodiscard]] auto bin () const {
        if constexpr (HasSum)
          return std::abs (this->sum) / k;
        else
          return std::abs ((*this)[0]) / k;  // NOLINT(clang-diagnostic-absolute-value)
      }

      [[nodiscard]] int cached_sum () const
        requires HasSum
      {
        return this->sum;
      }

    private:
      friend generic_partial_order<generic>;

      // template <template <typename> typename D, typename T>
      // requires IsGenericVector<D<T>> or IsGenericVector<D<T, 42>>
      // friend struct traits<D, T>;

      size_t k;
      std::conditional_t<EmbedsData, Data, Data*> datap;
  };
}
