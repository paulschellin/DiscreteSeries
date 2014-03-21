#ifndef DISCRETESERIES_HPP
#define DISCRETESERIES_HPP 1
#pragma once

/*
	For the redesigned DiscreteSeries class, a few changes will be made:
		
		[ ]	There won't be any iterators defined specifically for the class (no benefit)
		
		[ ]	No preprocessor macros will be used for optional compile-time things (too messy)
		[ ]	The template parameters will be simply typenames, not nested templates.

		[ ] Sampling rate most likely does not belong as a separate object, if multiple domain dimensions are to be allowed.

		[ ]	Delegate constructors will be used to limit annoying boilerplate code.

		[ ] pair iterator

		[ ] ostream_iterator


 */

#include <iostream>
#include <vector>
#include <complex>
#include <functional>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iterator>

#include <exception>
#include <stdexcept>

#include <boost/range.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/iterator/zip_iterator.hpp>


#include <boost/operators.hpp>

namespace PS {

	using ::boost::tuples::operator<<;


	//template <typename T1 = double, typename T2 = T1>
	
	template <typename SinglePassRange1, typename SinglePassRange2 = SinglePassRange1>
	class PlaceholderInterpolator {
		public:
			typedef typename std::iterator_traits< typename SinglePassRange1::iterator >::value_type	T1;
			typedef typename std::iterator_traits< typename SinglePassRange2::iterator >::value_type	T2;
		
			PlaceholderInterpolator (void) {}

			PlaceholderInterpolator (const SinglePassRange1& /*rng1*/, const SinglePassRange2& /*rng2*/)
			{ }

			~PlaceholderInterpolator (void) {}
			
			T2
			eval (const T1 input)
			{
				return T2(input);
			}

			T2
			operator() (const T1 input)
			{
				return eval(input);
			}
	};



	template <typename SinglePassRange1, typename SinglePassRange2 = SinglePassRange1>
	class NaiveLinearInterpolator {
		public:
			typedef typename std::iterator_traits< typename SinglePassRange1::iterator >::value_type	T1;
			typedef typename std::iterator_traits< typename SinglePassRange2::iterator >::value_type	T2;
		

			const SinglePassRange1 range_domain;
			const SinglePassRange2 range_codomain;


			NaiveLinearInterpolator (void) {}

			//template <typename SinglePassRange1, typename SinglePassRange2>
			NaiveLinearInterpolator (const SinglePassRange1& rng1, const SinglePassRange2& rng2)
				: range_domain(rng1)
				, range_codomain(rng2)
			{ }


			~NaiveLinearInterpolator (void) {}


			T2
			eval (const T1 input)
			{
				//	Check if input value is already equal to a bin
				if (boost::find(range_domain, input) != range_domain.end()) {
					return T2( input );
				}

				

				//auto lower = boost::lower_bound (range_domain, [input](T1 x){std::less_equal<T1>(input, x);});
				//auto upper = lower;
				//boost::advance(upper, 1);

				//auto upper = boost::upper_bound ();
				

				auto it1 = boost::adjacent_find(range_domain, [input](T1 x, T1 y){ return x <= input && input <= y; });


				T1 domain_lower (*it1);
				T1 domain_upper (*(it1+1));
				T2 codomain_lower (*(boost::begin(range_codomain) + std::distance(boost::begin(range_domain), it1)));
				T2 codomain_upper (*(boost::begin(range_codomain) + std::distance(boost::begin(range_domain), it1 + 1)));

				return T2(codomain_lower + (input - domain_lower) * (codomain_upper - codomain_lower) / (domain_upper - domain_lower));
			}


			T2
			operator() (const T1 input)
			{
				return eval(input);
			}

	};


	template <typename ContainerT>
	ContainerT
	make_iota (std::size_t n, typename ContainerT::value_type x0, typename ContainerT::value_type delta_x = typename ContainerT::value_type(1))
	{
		ContainerT result (n, x0);
		for (std::size_t i(1); i < n; ++i) {
			result[i] = result[i-1] + delta_x;
		}
		return result;
	}

/*
	template <typename ContainerT>
	ContainerT
	make_iota (ContainerT::value_type x_min, ContainerT::value_type x_max, ContainerT::value_type delta_x = ContainerT::value_type(1))
	{
		ContainerT result (n, x0);
		for (std::size_t i(1); i < n; ++i) {
			result[i] = result[i-1] + delta_x;
		}
		return result;
	}
*/



	template< typename DomainT //= double
			, typename CodomainT //= DomainT
			, typename ContainerT1 //= std::vector
			, typename ContainerT2 = ContainerT1
			//, typename template<typename ...ArgsX1> class ContainerT1 //= std::vector
			//, typename template<typename ...ArgsX2> class ContainerT2 //= ContainerT1
			//, typename ...Args1
			//, typename ...Args2
			, typename InterpolatorT = PlaceholderInterpolator<ContainerT1, ContainerT2>
			>

	/*
	template< typename ContainerDomainT
			, typename ContainerCodomainT = ContainerDomainT
			, typename InterpolatorT = PlaceholderInterpolator<ContainerDomainT, ContainerCodomainT>
			>
	*/
	class DiscreteSeries
	/*
		: boost::arithmetic1< DiscreteSeries<ContainerDomainT, ContainerCodomainT, InterpolatorT>
		, boost::arithmetic2< DiscreteSeries<ContainerDomainT, ContainerCodomainT, InterpolatorT>, ContainerCodomainT>
		>
	 */
	 	: boost::arithmetic1< DiscreteSeries< DomainT, CodomainT, ContainerT1, ContainerT2, InterpolatorT>
		, boost::arithmetic2< DiscreteSeries< DomainT, CodomainT, ContainerT1, ContainerT2, InterpolatorT>
							, ContainerT2 >
		>
	{
	public:
		
		typedef ContainerT1 ContainerDomainT;
		typedef ContainerT2 ContainerCodomainT;

		/*
		typedef ContainerT1<DomainT, Args1...>		ContainerDomainT;
		typedef ContainerT2<CodomainT, Args2...>	ContainerCodomainT;
		*/

		typedef DomainT		domain_type;
		typedef CodomainT	codomain_type;

		typedef ContainerT2	container_type;



		//typedef typename ContainerDomainT::value_type		domain_type;
		//typedef typename ContainerCodomainT::value_type		codomain_type;

		typedef std::pair <domain_type, codomain_type>		value_type;

		typedef typename container_type::reference			reference;
		typedef typename container_type::const_reference		const_reference;

		typedef typename container_type::pointer				pointer;
		typedef typename container_type::const_pointer		const_pointer;

		typedef typename container_type::iterator			iterator;
		typedef typename container_type::const_iterator		const_iterator;

		typedef typename container_type::reverse_iterator		reverse_iterator;
		typedef typename container_type::const_reverse_iterator	const_reverse_iterator;

		typedef typename container_type::difference_type		difference_type;
		typedef typename container_type::size_type			size_type;


		
		//
		//	Typedefs for the zip iterator (iteration by pairs)
		//

		typedef typename ContainerT1::iterator domain_iterator;
		typedef typename ContainerT1::const_iterator domain_const_iterator;

		typedef boost::tuple<domain_iterator, iterator> tuple_iterator;
		typedef boost::tuple<domain_const_iterator, const_iterator> tuple_const_iterator;

		typedef boost::zip_iterator<tuple_iterator> zipped_iterator;
		typedef boost::zip_iterator<tuple_const_iterator> zipped_const_iterator;




		//domain_type			samplingRate_;
		ContainerDomainT	domainValues_;
		ContainerCodomainT	codomainValues_;

		InterpolatorT	interpFn_;


		DiscreteSeries() = delete;


		struct _NoInit{};

		DiscreteSeries(std::size_t n)
			: DiscreteSeries(ContainerDomainT(n), ContainerCodomainT(n))
		{}

		DiscreteSeries(std::size_t n, _NoInit)
			: DiscreteSeries(ContainerDomainT(n), ContainerCodomainT(n))
		{}

		DiscreteSeries(std::size_t n, const domain_type samplingrate)		 
			: DiscreteSeries(ContainerDomainT(make_iota<ContainerDomainT>(n, 0., 1./samplingrate)), ContainerCodomainT(n))
		{}
		
		DiscreteSeries(std::size_t n, const domain_type samplingrate, _NoInit)
			: DiscreteSeries(ContainerDomainT(make_iota<ContainerDomainT>(n, 0., 1./samplingrate)), ContainerCodomainT(n))
		{}


		DiscreteSeries(const ContainerCodomainT& inR, const domain_type samplingrate)
			: DiscreteSeries(ContainerDomainT(make_iota<ContainerDomainT>(inR.size(), 0., 1./samplingrate)), ContainerCodomainT(inR))
		{}

		DiscreteSeries(const ContainerDomainT& inD, const ContainerCodomainT& inR)
			: domainValues_(inD)
			, codomainValues_(inR)
			, interpFn_(inD, inR)
		{
			if (!boost::is_sorted(domainValues_)) {
				throw std::invalid_argument("DiscreteSeries: The domain array was not strictly increasing (interpolators require this)!");
			}
		}


		~DiscreteSeries(void) { }




		//...



		//
		//	Iterators
		//
		/*
			Consider writing using the free functions instead?
			Also, should it expect cbegin and cend?
		 */
		iterator begin() { return codomainValues_.begin(); }
		const_iterator begin() const { return codomainValues_.begin(); }

		iterator end() { return codomainValues_.end(); }
		const_iterator end() const { return codomainValues_.end(); }

		reverse_iterator rbegin() { return reverse_iterator(end()); }
		const_reverse_iterator rbegin() const { return const_reverse_iterator(end()); }

		reverse_iterator rend() { return reverse_iterator(begin()); }
		const_reverse_iterator rend() const { return const_reverse_iterator(begin()); }


		const_iterator cbegin() const { return const_iterator(codomainValues_.begin()); }

		const_iterator cend() const { return const_iterator(codomainValues_.end()); }

		const_reverse_iterator crbegin() const { return const_reverse_iterator(end()); }

		const_reverse_iterator crend() const { return const_reverse_iterator(begin()); }


		zipped_iterator
		begin_zip () {
			return boost::make_zip_iterator (boost::make_tuple ( domainValues_.begin()
														, codomainValues_.begin() ) );
		}

		zipped_iterator
		end_zip () {
			return boost::make_zip_iterator (boost::make_tuple ( domainValues_.end()
														, codomainValues_.end() ) );
		}

		zipped_const_iterator
		begin_zip () const {
			return boost::make_zip_iterator (boost::make_tuple ( domainValues_.begin()
														, codomainValues_.begin() ) );
		}

		zipped_const_iterator
		end_zip () const {
			return boost::make_zip_iterator (boost::make_tuple ( domainValues_.end()
														, codomainValues_.end() ) );
		}

		//
		//	Capacity
		//

		size_type size() const { return codomainValues_.size(); }

		bool empty() const { return !size(); }


		
		//
		//	Element access
		//

		reference front() { return codomainValues_.front(); }
		const_reference front() const { return codomainValues_.front(); }

		reference back() { return codomainValues_.back(); }
		const_reference back() const { return codomainValues_.back(); }

		reference operator[](size_type n) { return codomainValues_[n]; }
		const_reference operator[](size_type n) const { return codomainValues_[n]; }

		reference at(size_type n) { return codomainValues_.at(n); }
		const_reference at(size_type n) const { return codomainValues_.at(n); }


		pointer data() { return codomainValues_.data(); }
		const_pointer data() const { return codomainValues_.data(); }	
		//value_type* data() { return codomainValues_.data(); }
		//const value_type* data() const { return codomainValues_.data(); }

		
		ContainerDomainT&
		get_domain (void) { return domainValues_; }

		ContainerCodomainT&
		get_codomain (void) { return codomainValues_; }

	
		codomain_type operator() (const domain_type x)
		{
			if (x < domainValues_.front() || x > domainValues_.back()) {
				throw std::domain_error("DiscreteSeries attempted to interpolate outside of the domain array!");
			}
			
			return interpFn_(x);
		}


		//
		//	Modifiers
		//

		//void swap()



		//
		//	Valarray compatability
		//

		value_type sum() const
		{
			return std::accumulate(codomainValues_.begin()+1, codomainValues_.end(), codomainValues_.front());
		}

		value_type min() const
		{
			return *boost::min_element(codomainValues_);
		}

		value_type max() const
		{
			return *boost::max_element(codomainValues_);
		}

		
		//	Shift elements left (returns copy)
		/*
		T
		shift(int n) const
		{
			T result(this->size());

			if (n >= 0) {
				if (n < this->size())
					std::copy(codomainValues_.begin() + n, codomainValues_.end(), result.begin());
			} else {
				if (-n < this->size())
					std::copy(codomainValues_.begin(), codomainValues_.end() + n, result.begin - n);
			}

			return result;
		}
		*/

		//	Circularly shift elements left (returns copy)
		/*
		T
		cshift(int m) const
		{
			T result(this->size());
			
			// Reduce m to an equivalent number in the range [0, size()).  We
			// have to be careful with negative numbers, since the sign of a % b
			// is unspecified when a < 0.
			
			long n(m);
			
			if (this->size() < numeric_limits<long>::max())
				n %= long(this->size());
			if (n < 0)
				n += this->size();


			std::copy(codomainValues_.begin(), codomainValues_.end(), result.begin() + (this->size() - n));
			
			std::copy(codomainValues_.begin() + n, codomainValues_.end(), result.begin());

			return result;
		}
		*/
		
		DiscreteSeries apply(value_type fn(value_type)) const
		{
			DiscreteSeries result(this->size());
			boost::transform(codomainValues_, result.begin(), fn);
			return result;
		}


		DiscreteSeries apply(value_type fn(const value_type&)) const
		{
			DiscreteSeries result(this->size());
			boost::transform(codomainValues_, result.begin(), fn);
			return result;
		}


		//
		//	Unary operators
		//

		DiscreteSeries operator+() const { return *this; }

		DiscreteSeries operator-() const
		{
			DiscreteSeries result(*this);
			boost::transform(codomainValues_, result.begin(), std::negate<value_type>());
			return result;
		}

		DiscreteSeries operator~() const
		{
			DiscreteSeries result(*this);
			boost::transform(codomainValues_, result.begin(), 
			[](const value_type x) { return ~x; }
			//std::bit_not<value_type>()
			);
			return result;
		}

		/*
		DiscreteSeries<ContainerDomainT, >
		operator!() const
		{
			DiscreteSeries<> result(*this);
			boost::transform(codomainValues_, result.begin(), std::not_equal_to<value_type>());
			return result;
		}
		*/

		
		//
		//	Scalar computed assignment
		//

		#define DISCRETESERIES_MAKE_SCALAR_OPERATOR( _op_ )		\
		DiscreteSeries& operator _op_ (const value_type& rhs) {	\
			for (auto val : codomainValues_)					\
				val _op_ rhs;									\
			return *this;}

		DISCRETESERIES_MAKE_SCALAR_OPERATOR( *= );
		DISCRETESERIES_MAKE_SCALAR_OPERATOR(/=);
		DISCRETESERIES_MAKE_SCALAR_OPERATOR(%=);
		DISCRETESERIES_MAKE_SCALAR_OPERATOR(+=);
		DISCRETESERIES_MAKE_SCALAR_OPERATOR(-=);
		DISCRETESERIES_MAKE_SCALAR_OPERATOR(^=);
		DISCRETESERIES_MAKE_SCALAR_OPERATOR(&=);
		DISCRETESERIES_MAKE_SCALAR_OPERATOR(|=);
		DISCRETESERIES_MAKE_SCALAR_OPERATOR(<<=);
		DISCRETESERIES_MAKE_SCALAR_OPERATOR(>>=);

		//
		//	Array computed assignment
		//

		#define DISCRETESERIES_MAKE_ARRAY_OPERATOR( _op_ )		\
		template <typename T>									\
		DiscreteSeries& operator _op_ (const T& rhs) {			\
			/*boost::transform(codomainValues_, rhs,			\
				codomainValues_.begin(), _op_ );	*/			\
			for (std::size_t i(0), N(codomainValues_.size()); i < N; ++i)	\
				codomainValues_[i] _op_ rhs[i];					\
			return *this;}

		DISCRETESERIES_MAKE_ARRAY_OPERATOR(*=);
		DISCRETESERIES_MAKE_ARRAY_OPERATOR(/=);
		DISCRETESERIES_MAKE_ARRAY_OPERATOR(%=);
		DISCRETESERIES_MAKE_ARRAY_OPERATOR(+=);
		DISCRETESERIES_MAKE_ARRAY_OPERATOR(-=);
		DISCRETESERIES_MAKE_ARRAY_OPERATOR(^=);
		DISCRETESERIES_MAKE_ARRAY_OPERATOR(&=);
		DISCRETESERIES_MAKE_ARRAY_OPERATOR(|=);
		DISCRETESERIES_MAKE_ARRAY_OPERATOR(<<=);
		DISCRETESERIES_MAKE_ARRAY_OPERATOR(>>=);

		


	};


	//
	//	Relational operators
	//
/*
	template <...>
	bool operator== (const & lhs, const & rhs)
	{
		//	Because boost::equal will use the .begin() and .end() on
		//	both lhs and rhs, it doesn't matter whether they are
		//	DiscreteSeries or any other STL container.
		return lhs.size() == rhs.size() && boost::equal(lhs, rhs);
	}

	template <...>
	bool operator!= (const & lhs, const & rhs)
	{
		return !(lhs == rhs);
	}
	
	template <...>
	bool operator< (const & lhs, const & rhs)
	{
		return boost::lexicographical_compare(lhs, rhs);
	}

	template <...>
	bool operator<= (const & lhs, const & rhs)
	{
		return !(rhs < lhs);
	}

	template <...>
	bool operator> (const & lhs, const & rhs)
	{
		return rhs < lhs;
	}

	template <...>
	bool operator>= (const & lhs, const & rhs)
	{
		return !(lhs < rhs);
	}



	template <...>
	void swap (& x, & y)
	{

	}
*/


	#define DISCRETESERIES_MAKE_MATH_FUNCTION(_math_fn_)	\
	template<typename ...Ts>								\
	inline DiscreteSeries<Ts...>							\
	_math_fn_ (const DiscreteSeries<Ts...>& xs) {			\
		using std:: _math_fn_ ;								\
		DiscreteSeries<Ts...> result (xs.size());			\
		boost::transform(xs, result.begin(), _math_fn_());	\
		return result;}


	DISCRETESERIES_MAKE_MATH_FUNCTION(abs);
	DISCRETESERIES_MAKE_MATH_FUNCTION(acos);
	DISCRETESERIES_MAKE_MATH_FUNCTION(asin);
	DISCRETESERIES_MAKE_MATH_FUNCTION(atan);
	DISCRETESERIES_MAKE_MATH_FUNCTION(atan2);
	DISCRETESERIES_MAKE_MATH_FUNCTION(cos);
	DISCRETESERIES_MAKE_MATH_FUNCTION(cosh);
	DISCRETESERIES_MAKE_MATH_FUNCTION(exp);
	DISCRETESERIES_MAKE_MATH_FUNCTION(log);
	DISCRETESERIES_MAKE_MATH_FUNCTION(log10);
	DISCRETESERIES_MAKE_MATH_FUNCTION(pow);
	DISCRETESERIES_MAKE_MATH_FUNCTION(sin);
	DISCRETESERIES_MAKE_MATH_FUNCTION(sinh);
	DISCRETESERIES_MAKE_MATH_FUNCTION(sqrt);
	DISCRETESERIES_MAKE_MATH_FUNCTION(tan);
	DISCRETESERIES_MAKE_MATH_FUNCTION(tanh);

	//	Note that the scalar atan2 functions are absent, as well
	//	as the scalar pow() functions.

/*
	template<typename Type, unsigned N, unsigned Last>
	struct tuple_printer {

    	static void print(std::ostream& out, const Type& value) {
    	    out << boost::get<N>(value) << ", ";
    	    tuple_printer<Type, N + 1, Last>::print(out, value);
		}
	};

	template<typename Type, unsigned N>
	struct tuple_printer<Type, N, N> {

		static void print(std::ostream& out, const Type& value) {
			out << boost::get<N>(value);
		}
	};
	
	template<typename ...Ts>
	::std::ostream&
	operator<<  ( ::std::ostream& os
				, const boost::tuple<Ts...>& rhs
				//, const typename DiscreteSeries<Ts...>::zipped_const_iterator::value_type& rhs
				//, const typename DiscreteSeries<Ts...>::zipped_const_iterator::value_type& rhs
				)
	{
		//os << boost::get<0>(rhs) << "\t" << boost::get<1>(rhs) << std::endl;
		//os << rhs. template get<0>() << "\t" << rhs. template get<1>() << std::endl;
		
		tuple_printer<boost::tuple<Ts...>, 0, sizeof...(Ts) - 1>::print(os, rhs);

		os << std::endl;
		
		return os;
	}
*/


	template<typename ...Ts>
	std::ostream& operator<< (std::ostream& os, const DiscreteSeries<Ts...>& rhs)
	{
		
		return os;
	}
	



}	// namespace PS

#endif
