#ifndef BLOCKCONTAINER_HPP
#define BLOCKCONTAINER_HPP
#include <vector>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/utility/enable_if.hpp>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include "VectorTraits.hpp"

namespace espresso {
  namespace esutil {
    template<class VectorClass, class Traits = VectorTraits<VectorClass> >
    class BlockVector {
    public:
      /* forward declarations */
      class ConstBlock;
      class Block;
      class const_iterator;
      class iterator;
      template<class, class, class, class> friend class BlockBase;

      /** structure containing the start and end of each block. */
      struct BlockBoundaries {
        size_t start, end;

	BlockBoundaries(size_t _start = 0, size_t _end = 0): start(_start), end(_end) {}
      };

      /** reference to a block. Behaves as much as possible like a
          reference to a normal VectorClass object. HOWEVER, after any
          resize of the underlying BlockVector, especially after resizing
	  ANY OTHER block, this reference is invalid. Resizing this block
	  itself is of course safe.

          The block reference is divided into two classes, Block and
          ConstBlock, which both inherit from BlockBase.  BlockBase
          provides the reference information and the const-functions
          of the vector interface. Block then adds the additional
          non-const functions.  ConstBlock adds just the possibilty to
          convert a Block into a ConstBlock.
      */
      template<class BlockVectorType, class BoundariesPtr,
	       class IteratorType, class ThinIteratorType, class ReferenceType>
      class BlockBase {
	template<class, class> friend class BlockVector;

      public:
        typedef typename Traits::value_type value_type;
	typedef typename Traits::size_type size_type;
        typedef typename Traits::difference_type difference_type;
        
        typedef typename Traits::iterator iterator;
        typedef typename Traits::thin_iterator thin_iterator;
        typedef typename Traits::const_iterator const_iterator;

        typedef typename Traits::reference reference;
        typedef typename Traits::thin_reference thin_reference;
        typedef typename Traits::const_reference const_reference;

        typedef typename Traits::pointer pointer;
        typedef typename Traits::thin_pointer thin_pointer;
        typedef typename Traits::const_pointer const_pointer;

      public:
        IteratorType begin() const { return Traits::template make_thick<IteratorType, ThinIteratorType>(&vector.data, vector.data.begin() + boundaries->start); }
        IteratorType end()   const { return Traits::template make_thick<IteratorType, ThinIteratorType>(&vector.data, vector.data.begin() + boundaries->end); }

        ReferenceType operator[](size_type n) const { return *(begin() + n); }
        ReferenceType at(size_type n) const { return vector.data.at(boundaries->start + n); }

        ReferenceType front() const { return *begin(); }
        ReferenceType back()  const { return *(end() - 1); }

        size_type size()     const { return boundaries->end - boundaries->start; }
        bool empty()         const { return size() == size_type(0); }
        size_type max_size() const { return vector.data.max_size(); }
        /// for a Block of a BlockVector, the capacity is the space to the next neighbor
        size_type capacity() const { return (boundaries+1)->start - boundaries->start; }

      protected:
        /// instances can only be obtained through BlockVector.
        BlockBase(BlockVectorType &_vector, size_t _index)
          : vector(_vector), boundaries(vector.boundaries.begin() + _index) {}
        /// instances can only be obtained through BlockVector.
        BlockBase(BlockVectorType &_vector, BoundariesPtr _boundaries)
          : vector(_vector), boundaries(_boundaries) {}

        /// reference to the underlying block vector
        BlockVectorType &vector;
        /// where the boundary information really is inside vector
        BoundariesPtr boundaries;
      };
      
      /** iterator over all blocks. Basically just the index of the
          block.

	  @tparam ReferenceType type we return as reference
	  (BlockReference or ConstBlockReference)
	  @tparam CRTP curiously recurrent template pattern
      */
      template<class BlockVectorType, class ReferenceType, class CRTP>
      class IteratorBase: public boost::iterator_facade<
	CRTP,
	ReferenceType,
	boost::random_access_traversal_tag,
	ReferenceType,
	ptrdiff_t>
      {
      public:
	/// reference type. Fixes up the iterator_facade for a reference only class
	typedef ReferenceType reference;

	/// default constructor, out-of-range iterator
	IteratorBase(): vector(0), index(0) {}

      protected:
	/// regular constructor
	IteratorBase(BlockVectorType *_vector, size_t _index): vector(_vector), index(_index) {}

        /// blockvector we are iterating over
        BlockVectorType *vector;
	/// index of current block
	size_t index;

      private:
	template <class,class,class> friend class IteratorBase;
	friend class boost::iterator_core_access;

	template<class OtherBlockVectorType, class OtherReferenceType, class OtherIteratorType>
	bool equal(const IteratorBase<OtherBlockVectorType, OtherReferenceType, OtherIteratorType> &other) const {
	  return index == other.index;
	}

	void increment()       { index++; }
	void decrement()       { index--; }
	void advance(size_t n) { index += n; }
	ReferenceType dereference() const { return ReferenceType(*vector, index); }
	ptrdiff_t distance_to(const IteratorBase<BlockVectorType, ReferenceType, CRTP> &other) const {
	  return other.index - index;
	}
      };

    public:
      /// constant block. Behaves mostly like a const VectorClass.
      class ConstBlock: public BlockBase< const BlockVector,
					  typename std::vector<BlockBoundaries>::const_iterator,
					  typename Traits::const_iterator,
					  typename Traits::const_iterator,
					  typename Traits::const_reference > {
      public:
	typedef typename BlockVector::const_iterator BlockIterator;
	typedef BlockBase< const BlockVector,
			   typename std::vector<BlockBoundaries>::const_iterator,
			   typename Traits::const_iterator,
			   typename Traits::const_iterator,
			   typename Traits::const_reference > BlockBaseType;
      
	ConstBlock(const BlockVector &_vector, size_t _index):
	  BlockBaseType(_vector, _index) {};

        // const->non-const conversion
        ConstBlock(const Block &_block):
          BlockBaseType(_block.vector, _block.boundaries) {};
      };

      /// a block. Behaves mostly like a VectorClass.
      class Block: public BlockBase< BlockVector,
				     typename std::vector<BlockBoundaries>::iterator,
                                     typename Traits::iterator,
                                     typename Traits::thin_iterator,
                                     typename Traits::reference > {
        // for const->nonconst conversion
        friend class ConstBlock;

      public:
        typedef typename BlockVector::iterator BlockIterator;
	typedef typename Traits::iterator Iterator;
	typedef typename Traits::thin_iterator ThinIterator;
	typedef typename Traits::const_reference ConstReference;
	typedef BlockBase< BlockVector,
			   typename std::vector<BlockBoundaries>::iterator,
			   typename Traits::iterator,
			   typename Traits::thin_iterator,
			   typename Traits::reference > BlockBaseType;
        typedef typename Traits::size_type size_type;

        Block(BlockVector &_vector, size_t _index): BlockBaseType(_vector, _index) {};

        /// copy one block to another
        Block &operator=(const Block &_block) {
          resize(_block.size());
          assign(_block.begin(), _block.end());
          return *this;
        }

        /// copy one block to another
        Block &operator=(const ConstBlock &_block) {
          resize(_block.size());
          assign(_block.begin(), _block.end());
          return *this;
        }

        /// set to identical values
        void assign(size_type n, ConstReference v) { clear(); insert(begin(), n, v); }

        /** assign a sequence. Note that the iterators need to be
	    random iterators, in contrast to standard STL containers.
            Moreover, this uses std::copy - some vectors like TupleVector
            impose further restrictions on std::copy, for example that
            It needs to be one of its own iterators.

            Remark regarding the second argument: both start and end
            should be iterators of type It. The strange form of the
            second argument comes from an ambiguity in case
            VectorClass::const_reference is a reference to an integer
            (e.g. VectorClass = std::vector<int>). Now, if end would
            just be of type It and you would want to use the other
            assign(3,42) to set this vector to 3 elements of value 42
            (the other assign-variant), this would result in an
            error. Why? Because the compiler will also try this
            template-function with type It=int or similar, which of
            course fails. Therefore, we need to disable this variant of
            assign in case It is an integer, which is precisely what
            disable_if_c does for us.
        */
	template<class It>
	void assign(It start,
                    typename boost::disable_if_c<std::numeric_limits<It>::is_integer, It>::type end) {
          resize(end - start);
          std::copy(start, end, begin());
        }

        void push_back(ConstReference v) {
          resize(size() + 1);
          back() = v;
        }
        void pop_back() { resize(size() - 1); }

        ThinIterator insert(ThinIterator pos, ConstReference v) {
          size_type index = pos - this->begin(); // resize safety
          insert(pos, 1, v);
          return this->begin() + index;
        }

        void insert(ThinIterator pos, size_type n, ConstReference v) {
          pos = makeSpace(pos, n);
          typename Traits::iterator fpos = Traits::template make_thick<Iterator, ThinIterator>(&(this->vector.data), pos);
          std::fill(fpos, fpos + n, v);
        }
        /** insert a sequence. Note that the iterators need to be
	    random iterators, in contrast to standard STL containers.
	    Moreover, this uses std::copy - some vectors like
	    TupleVector impose further restrictions on std::copy, for
	    example that It needs to be one of its own iterators.

            For explanations of the second iterators type (which in
            fact is just It), see the assign template.
        */
        template<class It>
        void insert(ThinIterator pos, It begin,
                    typename boost::disable_if_c<std::numeric_limits<It>::is_integer, It>::type end) {
          pos = makeSpace(pos, (end - begin));
          std::copy(begin, end, Traits::template make_thick<Iterator, ThinIterator>(&(this->vector.data), pos));
        }

        ThinIterator erase(ThinIterator pos) { return erase(pos, pos+1); }
        ThinIterator erase(ThinIterator begin, ThinIterator end) {
          size_type index = begin - this->begin(); // resize safety
          std::copy(end, ThinIterator(this->end()), this->begin() + index);
          resize(size() - (end - begin));
          return this->begin() + index;
        }

        void clear() {resize(0); }

        /** To allow for containers that do no have a value_type, e.g. TupleVector,
            this resize does not take a value_type argument. Could be added, if necessary. */
        void resize(size_type newsize) {
          if (newsize > BlockBaseType::capacity()) {
            reserve(newsize);
          }
          // assign block its new size, now there is space
	  BlockBaseType::boundaries->end =
	    BlockBaseType::boundaries->start + newsize;
	}

        void reserve(size_type newsize) {
          /* check to reduce the number of function calls to the not
             inlined vector.reserveForBlock. Most of the time, we do
             not need to reserve here.*/
          if (newsize > BlockBaseType::capacity()) {
            BlockBaseType::vector.reserveForBlock(BlockBaseType::boundaries, newsize);
          }
        }

      private:

        /** make space for n elements starting at pos. Will shift all the following data by
            n positions. Returns the new position of pos, which my change due to resizing. */
        ThinIterator makeSpace(ThinIterator pos, size_type n);
      };
 
      /// a random iterator over blocks.
      class iterator: public IteratorBase<BlockVector, Block, iterator> {
	// class that can generate an iterator
	template<class, class> friend class BlockVector;
	friend class BlockVector::Block;
	// for const->nonconst conversion
	friend class const_iterator;

      public:
	typedef IteratorBase<BlockVector, Block, iterator> IteratorBaseType;

	/// default constructor, out-of-range iterator
	iterator() {}

      private:
	/// constructable only through BlockVector
	iterator(BlockVector *_vector, size_t _index)
	  : IteratorBaseType(_vector, _index) {}
	/// for &Block
	iterator(BlockVector &_vector, size_t _index)
	  : IteratorBaseType(&_vector, _index) {}
      };

      /// a random iterator over constant blocks.
      class const_iterator: public IteratorBase<const BlockVector, ConstBlock, const_iterator> {
	// class that can generate an iterator
	template<class, class> friend class BlockVector;
	friend class ConstBlock;
      public:
	typedef IteratorBase<const BlockVector, ConstBlock, const_iterator> IteratorBaseType;

	/// default constructor
	const_iterator() {}

	/// non-const->const conversion
	const_iterator(const iterator &_it)
	  : IteratorBaseType(_it.vector, _it.index) {}

      private:
	/// constructable only through BlockVector
	const_iterator(const BlockVector *_vector, size_t _index)
	  : IteratorBaseType(_vector, _index) {}
	/// for &ConstBlock
	const_iterator(const BlockVector &_vector, size_t _index)
	  : IteratorBaseType(&_vector, _index) {}
      };
      
      typedef iterator pointer;
      typedef const_iterator const_pointer;
      
      /// type of block indices, can in theory differ from vectorClass::size_type
      typedef size_t size_type;
      /// type of block index differences, can in theory differ from vectorClass::difference_type
      typedef ptrdiff_t difference_type;

      /******************************************************************
       * Methods of BlockVector
       ******************************************************************/

      BlockVector(size_t nBlocks = 0): boundaries(1), gapSize(2)
      { resize(nBlocks); }

      void resize(size_t newsize);

      iterator       begin()       { return       iterator(this, 0); }
      const_iterator begin() const { return const_iterator(this, 0); }

      iterator       end()         { return       iterator(this, size()); }
      const_iterator end()   const { return const_iterator(this, size()); }

      Block      operator[](size_t n)       { return *(begin() + n); }
      ConstBlock operator[](size_t n) const { return *(begin() + n); }

      Block      at(size_t n)       {
        if (n >= size()) {
	  throw std::out_of_range("BlockVector::at");
	}
        return *(begin() + n);
      }
      ConstBlock at(size_t n) const {
        if (n >= size()) {
	  throw std::out_of_range("BlockVector::at");
	}
        return *(begin() + n);
      }

      Block      front()       { return *begin(); }
      ConstBlock front() const { return *begin(); }

      Block      back()       { return *(end() - 1); }
      ConstBlock back() const { return *(end() - 1); }

      /// insert an empty block
      iterator insert(iterator pos) {
        insert(pos, 1);
        return pos; // our block iterator is resize-safe
      }
      /// insert a copy of another block
      void insert(iterator pos, ConstBlock block) {
        pos = insert(pos);
        *pos = block;
      }
      /// insert n empty blocks
      void insert(iterator pos, size_t n) {
        resize(size() + n);
        const_iterator stop = end() - n;
        std::copy_backward(const_iterator(pos), stop, end());
        for (iterator it = pos; it != stop; ++it) {
          (*pos).clear();
        }
      }
      
      iterator erase(iterator pos) { return erase(pos, pos + 1); }
      iterator erase(iterator begin, iterator end) {
        size_t n = end - begin;
        std::copy(end, this->end(), begin);
        resize(size() - n);
        return begin;
      }

      void clear() { resize(0); }

      size_type size() const { return boundaries.size() - 1; }
      bool empty() const { return size() == size_type(0); }
      size_type max_size() const { return boundaries.max_size(); }

      /// get total number of elements in all blocks
      size_t data_size() const;
      size_t target_gap_size() const { return gapSize; }

      /** access to the raw data storage. Use this e.g. to add properties to
	  a TupleVector.  */
      VectorClass &raw() { return data; }
      const VectorClass &raw() const { return data; }

    private:
      /** block boundaries. The last entry is not really a block; it simply serves to speed up
	  the size calculation for blocks, and contains as start and end the size of data
      */
      std::vector<BlockBoundaries> boundaries;
      VectorClass data;
      int gapSize;

      void reserveForBlock(typename std::vector<BlockBoundaries>::iterator block, typename Traits::size_type newsize);

      void resizeDataBuffer(size_type newsize) {
        data.resize(newsize);
	BlockBoundaries &endbuffer = boundaries.back();
	endbuffer.start = endbuffer.end = data.size();
      }
    };

    /***********************************************************************
     * TEMPLATE IMPLEMENTATION
     ***********************************************************************/

    template<class VectorClass, class Traits>
    size_t BlockVector<VectorClass, Traits>::data_size() const {
      size_t s = 0;
      const_iterator stop = end();
      for (const_iterator it = begin(); it != stop; ++it) {
	s += it->size();
      }
      return s;
    }

    template<class VectorClass, class Traits>
    void
    BlockVector<VectorClass, Traits>::resize
    (size_t newsize) {
      size_t oldSize = size();

      boundaries.resize(newsize + 1);

      if (oldSize < newsize) {
	/* Where to start with the new blocks. They have size 0, but
	   the gap size is the default.  If there were blocks before,
	   leave the default gap to the last one, otherwise we start
	   at 0 */
	size_t blockLocation = oldSize > 0 ? boundaries[oldSize - 1].end + gapSize : 0;
        typename std::vector<BlockBoundaries>::const_iterator stop = boundaries.end();
	for (typename std::vector<BlockBoundaries>::iterator it = boundaries.begin() + oldSize;
	     it != stop; ++it) {
	  it->start = it->end = blockLocation;
	  blockLocation += gapSize;
	}
	resizeDataBuffer(blockLocation - gapSize);
      }
      else {
	// maybe free some memory
	size_t newDataSize = (boundaries.begin() + newsize)->end + gapSize;
	if (data.size() > newDataSize) {
	  resizeDataBuffer(newDataSize);
	}
      }
    }

    template<class VectorClass, class Traits>
    void
    BlockVector<VectorClass, Traits>::reserveForBlock
    (typename std::vector<BlockVector::BlockBoundaries>::iterator block,
     typename Traits::size_type newsize) {
      /* Double the gap size, in the style of a std::vector.
	 There, however, not the gapsize is doubled, but the capacity.
      */
      gapSize *= 2;

      size_t index = block - boundaries.begin();
      if (newsize > (*this)[index].capacity()) {
	/* now it gets messy, we need to shift away all the rest blocks.
	   At this occasion, we also create space in all the other blocks */

	// determine how much space we need to for the remaining blocks
	size_t spaceNeeded = newsize;
	for (iterator it = begin() + index + 1; it != end(); ++it) {
	  spaceNeeded += it->size();
	}
	// add space for gaps - all gaps are reset to the default gap size
	spaceNeeded += (size() - index)*gapSize;

	resizeDataBuffer(boundaries[index].start + spaceNeeded);

	typename Traits::size_type space = data.size();

	// and now, copy all other blocks to their new location, and update
	// location information
	for (size_t b = size() - 1; b > index; --b) {
	  Block block  = (*this)[b];
	  size_t size = block.size();
	  std::copy_backward(block.begin(), block.end(), data.begin() + space - gapSize);

	  space -= gapSize;
	  boundaries[b].end = space;
	  space -= size;
	  boundaries[b].start = space;
	}
      }
    }

    template<class VectorClass, class Traits>
    typename Traits::thin_iterator
    BlockVector<VectorClass, Traits>::Block::makeSpace
    (typename Traits::thin_iterator pos, size_type n) {
      // Iterator->index since STL iterators are not resize-safe
      size_type startId = pos - begin();
      resize(size() + n);
      // pos after resizing
      pos = begin() + startId;
      std::copy_backward(pos, typename Traits::thin_iterator(end() - n), end());
      return pos;
    }
  }
}

#endif
