#ifndef DISCRETESERIES_HPP
#define DISCRETESERIES_HPP 1

#include <iostream>
#include <vector>
#include <complex>
#include <functional>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iterator>

#include <boost/range.hpp>

#include <boost/operators.hpp>

//#include <inplace_transform/inplace_transform.hpp>

#define DISCRETESERIES_INTERPOLATION 1

#ifdef DISCRETESERIES_INTERPOLATION
	//#include <GSLInterpolationExtension/CubicInterpFunc.hpp>
#endif

//	Options (define preprocessor macros to enable the following):
//	DISCRETESERIES_INTERPOLATION	Enables interpolation. (currently broken)
//	DISCRETESERIES_BOOST_UNITS		Use boost::units. Not yet implemented.


/******************************************************************************
 *	Due to the unpredictable availability of C++11 support, during			  *
 *	compilation the preprocessor macro "__cplusplus" is checked to see if the *
 *	compiler supports the C++11 standard library. If __cplusplus > 199711L,   *
 *	the compiler supports C++11 and so the standard library is used. 		  *
 *	Otherwise, the equivalent Boost libraries are required.					  *
 *																			  *
 *	To force the use of the Boost libraries regardless of C++11 compiler 	  *
 *	support, define the "DISCRETESERIES_USE_BOOST_LIBS" preprocessor macro at *
 *	compile time by passing -DDISCRETESERIES_USE_BOOST_LIBS as a flag to	  *
 *	the compiler's invocation.												  *
 *****************************************************************************/


//#define DISCRETESERIES_USE_STD_FUNCTION 1
#define DISCRETESERIES_USE_BOOST_FUNCTION 1

#ifdef DISCRETESERIES_USE_BOOST_FUNCTION
#include <boost/function.hpp>
#endif


#define DISCRETESERIES_USE_BOOST_SERIALIZE 1

#ifdef DISCRETESERIES_USE_BOOST_SERIALIZE

	#include <boost/archive/text_oarchive.hpp>
	#include <boost/archive/text_iarchive.hpp>
	#include <boost/serialization/version.hpp>
	#include <boost/serialization/vector.hpp>
	#include <boost/serialization/string.hpp>
	#include <boost/archive/xml_oarchive.hpp>
	#include <boost/archive/xml_iarchive.hpp>
	#include <boost/serialization/nvp.hpp>
	#include <boost/serialization/utility.hpp>


#endif

namespace PS {
	using std::cout;
	using std::endl;
	using std::vector;
	using std::complex;
	
	//using std::bind;
	using std::transform;
	using std::accumulate;
	using std::inner_product;
	using std::generate_n;
	using std::for_each;
	//using namespace std::placeholders;
	
	#ifdef DISCRETESERIES_USE_STD_FUNCTION
		using std::function;
	#elif defined (DISCRETESERIES_USE_BOOST_FUNCTION)
		using ::boost::function;
	#endif
	
	//	In order to allow complex functions to return both complex and real
	//	values, we must define the complex operators for real values.
	//	But we probably don't need them anyway.
	/*
	double real (const double input) { return input; }
	double imag (const double input) { return 0.; }
	double abs (const double input) { return input; }
	double arg (const double input) { return 0.; }
	double norm (const double input) { return input * input; }
	double conj (const double input) { return input; }
	//	??? polar()	???
	*/
	
#ifdef DISCRETESERIES_CUSTOM_ITERATORS
	
	template<class Item, int N>
	class DiscreteSeriesIterator
		: public std::iterator<std::output_iterator_tag,	//	Category
								Item,						//	value_type
								ptrdiff_t,					//	difference_type
								Item*,						//	pointer
								Item&						//	reference
		> {
			
			pointer itemPtr_;
		public:
			
			DiscreteSeriesIterator(pointer toCopy)
				: itemPtr_(toCopy)
			{ }
			
			DiscreteSeriesIterator(const DiscreteSeriesIterator& toCopy)
				: itemPtr_(toCopy)
			{ }
			
			
			DiscreteSeriesIterator&
			operator++(void)
			{
				++itemPtr_;
				return *this;
			}
			
			
			DiscreteSeriesIterator
			operator++(int)
			{
				DiscreteSeriesIterator tmp(*this);
				operator++();
				return tmp;
			}
			
			bool
			operator==(const DiscreteSeriesIterator& rhs) const
			{ return (itemPtr_ == rhs.itemPtr_); }
			
			bool
			operator!=(const DiscreteSeriesIterator& rhs) const
			{ return (itemPtr_ != rhs.itemPtr_); }
			
			reference
			operator*(void)
			{ return *itemPtr_; }
			
			const reference
			operator*(void) const
			{ return *static_cast<const pointer>(itemPtr_); }
			
			reference
			operator->(void)
			{ return *itemPtr_; }
			
			const reference
			operator->(void) const
			{ return *static_cast<const pointer>(itemPtr_); }
	};
	
	
#endif




	template<typename T1 = double, typename T2 = T1>
	class FakeInterpolator {
	private:
		
		

	
	public:
	
		FakeInterpolator (void) {}

		template<typename SinglePassRange1, typename SinglePassRange2>
		FakeInterpolator (const SinglePassRange1& rng1, const SinglePassRange2& rng2)
		{
			
		}

		~FakeInterpolator (void) {}

		T2
		eval (const T1 input)
		{
			
		}


		T2
		operator() (const T1 input)
		{
			return eval(input);
		}
	};





/*
	Using partial template specialization, the template parameter of
	std::complex<> can be known, since operators are overloaded 
 
 */

template< template<typename... /*T1*/> class ContainerDomainT
		, template<typename... /*T2*/> class ContainerCodomainT = ContainerDomainT
		, typename T1 = double
		, typename T2 = T1
		, typename InterpT = FakeInterpolator<double> //CubicInterpFunc<T1, T2>
		>

	//template<typename Domain, typename Range, typename AllocD = std::allocator<Domain>, typename AllocR = std::allocator<Range> >
	class DiscreteSeries
		:	boost::arithmetic1 < DiscreteSeries<ContainerDomainT, ContainerCodomainT, T1, T2, InterpT>
		,	boost::arithmetic2 < DiscreteSeries<ContainerDomainT, ContainerCodomainT, T1, T2, InterpT>, ContainerCodomainT<T2>
		> >
	{
	public:
		
		typedef ContainerDomainT<T1>	domain_container_type;
		typedef ContainerCodomainT<T2> codomain_container_type;

		typedef typename domain_container_type::value_type	domain_type;
		typedef typename codomain_container_type::value_type codomain_type;

		typedef DiscreteSeries DSeriesT;
		
		typedef typename codomain_container_type::value_type			value_type;
		
		typedef typename codomain_container_type::reference 			reference;
		typedef typename codomain_container_type::const_reference	const_reference;

		typedef typename codomain_container_type::pointer			pointer;
		typedef typename codomain_container_type::const_pointer		const_pointer;
		
		typedef typename codomain_container_type::iterator			iterator;
		typedef typename codomain_container_type::const_iterator		const_iterator;
		
		typedef typename codomain_container_type::reverse_iterator	reverse_iterator;
		typedef typename codomain_container_type::const_reverse_iterator const_reverse_iterator;
		
		typedef typename codomain_container_type::difference_type	difference_type;
		typedef typename codomain_container_type::size_type			size_type;
		



		domain_type		samplingRate_;
		
		domain_container_type	domainValues_;
		codomain_container_type	codomainValues_;
		

		InterpT interpFn_;
		//CubicInterpFunc<domain_type, codomain_type> interpFn_;


	private:
		
		#ifdef DISCRETESERIES_USE_BOOST_SERIALIZE

		friend class boost::serialization::access;
		template <class Archive>
		void
		serialize(Archive& ar, const unsigned int /*version*/)
		{
			ar & BOOST_SERIALIZATION_NVP(samplingRate_)
				& BOOST_SERIALIZATION_NVP(domainValues_)
				& BOOST_SERIALIZATION_NVP(codomainValues_);
		}

		#endif


	public:



		DiscreteSeries (void)
			: samplingRate_(0.)
			, domainValues_(0)
			, codomainValues_(0)
			, interpFn_()
		{
			// Don't init the interpolation function since the length may change.
			//initInterpFn();
		}
		
		DiscreteSeries (unsigned length)
			: samplingRate_(0.)
			, domainValues_(length)
			, codomainValues_(length)
			, interpFn_()
		{
			initInterpFn();
		}
		
		DiscreteSeries (const unsigned length, const domain_type samplingrate)
			: samplingRate_(samplingrate)
			, domainValues_(length)
			, codomainValues_(length)
			, interpFn_()
		{
			initInterpFn();
		}
		
		DiscreteSeries (const domain_container_type& inD, const codomain_container_type& inR)
			: samplingRate_(0.)
			, domainValues_(inD)
			, codomainValues_(inR)
			, interpFn_()
		{
			initInterpFn();
		}
		
		DiscreteSeries (const codomain_container_type& inR, const domain_type& samplingrate)
			: samplingRate_(samplingrate)
			, domainValues_(inR.size())
			, domainValues_(make_transformed_sequence_iterator(domain_type(0), std::bind2nd(std::divides<domain_type>(), samplingrate)), make_transformed_sequence_iterator(domain_type(inR.size()), std::bind2nd(std::divides<domain_type>(), samplingrate)))
			, codomainValues_(inR)
			, interpFn_()
		{
			initInterpFn();
		}
		
		DiscreteSeries (const DSeriesT& toCopy)
			: samplingRate_(toCopy.samplingRate_)
			, domainValues_(toCopy.domainValues_)
			, codomainValues_(toCopy.codomainValues_)
			, interpFn_()
		{
			initInterpFn();
		}
		
		~DiscreteSeries (void)
		{ }
		
		
		
		void
		initInterpFn (void)
		{
			#ifdef DISCRETESERIES_INTERPOLATION
			//interpFn_.setAndInit(domainValues_, codomainValues_);
			#endif
		}
		
		int
		SetEvenCombSpacing (void)
		{
			for (unsigned i = 0; i < domainValues_.size(); ++i) {
				domainValues_[i] = double(i) / samplingRate_;
			}
			return 0;
		}
		
		size_type
		size(void) const
		{
			return size_type (domainValues_.size());
		}
		
		
		iterator
		begin(void)
		{
			return iterator(codomainValues_.begin());
		}
		
		iterator
		end(void)
		{
			return iterator(codomainValues_.end());
		}
		
		const_iterator
		begin(void) const
		{
			return const_iterator(codomainValues_.begin());
		}
		
		const_iterator
		end(void) const
		{
			return const_iterator(codomainValues_.end());
		}
		
		const_iterator
		cbegin(void) const
		{
			#ifndef DISCRETESERIES_USE_CBEGIN_CEND
				return const_iterator(codomainValues_.begin());
			#else
				return const_iterator(codomainValues_.cbegin());
			#endif
		}
		
		const_iterator
		cend(void) const
		{
			#ifndef DISCRETESERIES_USE_CBEGIN_CEND
				return const_iterator(codomainValues_.end());
			#else
				return const_iterator(codomainValues_.cend());
			#endif
		}
		
		////////////////////////////////////////////////////////////////////////
		//					Arithmetic Operator Overloads
		////////////////////////////////////////////////////////////////////////
		//	Because many of the operators upon the series will be with another
		//	series over an identical domain, for the sake of performance there
		//	are two versions of each operation: one for identical domains,
		//	and one for non-identical domains.
		//
		//	Identical Domains:
		//		Perform an element-wise operation on the DiscreteSeries object.
		//
		//	Different Domains:
		//		The procedure is as follows:
		//		1.	Find the minimum domain value of both *this/lhs and rhs.
		//		2.	Find the maximum domain value of both *this/lhs and rhs.
		//		3.	Determine the step size of total domain (Use minimum step
		//				size or maximum??).
		//		4.	Perform the operator over the new domain for the two ranges.
		//		
		//
		////////////////////////////////////////////////////////////////////////
		
		DSeriesT&
		operator+= (const DSeriesT& rhs)
		{
			if (DomainsEqual(rhs)) {
				boost::transform(codomainValues_, rhs, codomainValues_.begin(), std::plus<codomain_type>());
				//PS::inplace_transform(codomainValues_, rhs, std::plus<Range>());
			} else {
				CombineDomainsWith(rhs);

				codomain_container_type copyOfOldRange (codomainValues_);
				
				codomainValues_.resize(domainValues_.size());
				
				#ifdef DISCRETESERIES_INTERPOLATION
				for (unsigned i = 0; i < domainValues_.size(); ++i) {
					codomainValues_[i] = codomainValues_.InterpValueAt(domainValues_[i]);
					codomainValues_[i] += rhs.InterpValueAt(domainValues_[i]);
				}
				#endif
			}
			return *this;
		}
		
		DSeriesT&
		operator-= (const DSeriesT& rhs)
		{
			if (DomainsEqual(rhs)) {
				boost::transform(codomainValues_, rhs, codomainValues_.begin(), std::minus<codomain_type>());
			} else {
				#ifdef DISCRETESERIES_INTERPOLATION
				for (unsigned i = 0; i < domainValues_.size(); ++i) {
					codomainValues_[i] -= rhs.InterpValueAt(domainValues_[i]);
				}
				#endif
			}
			return *this;
		}
		
		DSeriesT&
		operator*= (const DSeriesT& rhs)
		{
			if (DomainsEqual(rhs)) {
				boost::transform(codomainValues_, rhs, codomainValues_.begin(), std::multiplies<codomain_type>());
			} else {
				#ifdef DISCRETESERIES_INTERPOLATION
				for (unsigned i = 0; i < domainValues_.size(); ++i) {
					codomainValues_[i] *= rhs.InterpValueAt(domainValues_[i]);
				}
				#endif
			}
			return *this;
		}
		
		DSeriesT&
		operator/= (const DSeriesT& rhs)
		{
			if (DomainsEqual(rhs)) {
				boost::transform(codomainValues_, rhs, codomainValues_.begin(), std::divides<codomain_type>());
			} else {
				#ifdef DISCRETESERIES_INTERPOLATION
				for (unsigned i = 0; i < domainValues_.size(); ++i) {
					codomainValues_[i] /= rhs.InterpValueAt(domainValues_[i]);
				}
				#endif
			}
			return *this;
		}
		
		DSeriesT&
		operator+= (const codomain_type& rhs)
		{
			//PS::range::inplace_transform(codomainValues_, std::bind2nd(std::plus<Range>(), rhs));
			for (unsigned i = 0; i < codomainValues_.size(); ++i) {
				codomainValues_.at(i) += rhs;
			}
			return *this;
		}
		
		DSeriesT&
		operator-= (const codomain_type& rhs)
		{
			//PS::range::inplace_transform(codomainValues_, std::bind2nd(std::minus<Range>(), rhs));
			for (unsigned i = 0; i < codomainValues_.size(); ++i) {
				codomainValues_.at(i) -= rhs;
			}
			return *this;
		}
		
		DSeriesT&
		operator*= (const codomain_type& rhs)
		{
			//PS::range::inplace_transform(codomainValues_, std::bind2nd(std::multiplies<Range>(), rhs));
			for (unsigned i = 0; i < codomainValues_.size(); ++i) {
				codomainValues_.at(i) *= rhs;
			}
			return *this;
		}
		
		DSeriesT&
		operator/= (const codomain_type& rhs)
		{
			//PS::range::inplace_transform(codomainValues_, std::bind2nd(std::divides<Range>(), rhs));
			for (unsigned i = 0; i < codomainValues_.size(); ++i) {
				codomainValues_.at(i) /= rhs;
			}
			return *this;
		}
		
		reference
		operator[] (const size_type n)
		{
			return codomainValues_[n];
		}
		
		const_reference
		operator[] (const size_type n) const
		{
			return codomainValues_[n];
		}
		
		reference
		at(const size_type n)
		{
			return codomainValues_.at(n);
		}
		
		const_reference
		at(const size_type n) const
		{
			return codomainValues_.at(n);
		}
		
		domain_type&
		atDomain (const std::size_t n)
		{
			return (domainValues_.at(n));
		}
		
		const domain_type&
		atDomain (const std::size_t n) const
		{
			return (domainValues_.at(n));
		}
		
#ifdef DISCRETESERIES_INTERPOLATION
		codomain_type
		InterpValueAt (const domain_type& x) const
		{
			return interpFn_.eval(x);
		}
		
		codomain_container_type
		InterpValues (const domain_container_type& xs) const
		{
			codomain_container_type result (xs.size());
			boost::transform(xs, result.begin(), interpFn_);
			return result;
		}
#endif

		codomain_type&
		operator() (const domain_type& x)
		{
			#ifdef DISCRETESERIES_INTERPOLATION
				return InterpValueAt(x);
			#else
				return codomainValues_[unsigned(x)];
			#endif
		}
		
		const codomain_type&
		operator() (const domain_type& x) const
		{
			#ifdef DISCRETESERIES_INTERPOLATION
				return InterpValueAt(x);
			#else
				return codomainValues_[unsigned(x)];
			#endif
		}
		
		codomain_type&
		normalize (const codomain_type normfactor)
		{
			codomainValues_ *= normfactor;
			//for (unsigned i = 0; i < codomainValues_.size(); ++i) {
			//	codomainValues_[i] *= normfactor;
			//}
			
			//PS::inplace_transform(codomainValues_, std::bind2nd(std::divides<Range>(), normfactor));
			return codomainValues_;
		}
		
		bool
		DomainsEqual (const DSeriesT& rhs) const
		{
			if ((domainValues_.size() != rhs.domainValues_.size()) && (samplingRate_ != rhs.samplingRate_)){
				return false;
			}
			//return std::equal(domainValues_.begin(), domainValues_.end(), rhs.domainValues_.begin());
			return boost::equal(domainValues_, rhs.domainValues_);
		}
		
		
		domain_container_type
		GetCombinedDomain (const DSeriesT& lhs, const DSeriesT& rhs) const
		{
			//	Get the two domain sizes 
			std::size_t lhsSize = lhs.size() - 1;
			std::size_t rhsSize = rhs.size() - 1;
			
			domain_type domainMin = std::min(lhs.atDomain(0), rhs.atDomain(0));
			domain_type domainMax = std::max(lhs.atDomain(lhsSize), rhs.atDomain(rhsSize));
			
			//	Get the sampling rate
			//	???	- Should it be the minimum rate, the maximum rate, the least
			//		common factor? Should an offset be included to prevent
			//		issues where the minimum and maximum values are not
			//		multiples of the sampling rate?
			domain_type samplingRate = std::max(lhs.GetSamplingRate(), rhs.GetSamplingRate());
			//Domain samplingRate = std::min(lhs.GetSamplingRate(), rhs.GetSamplingRate());
			
			std::size_t numSamples ((domainMax - domainMin) / samplingRate + 1);
			
			domain_container_type result (numSamples);
			
			for (std::size_t i = 0; i != numSamples; ++i) {
				result[i] = domain_type(i) / samplingRate_;
			}
			
			return result;
		}
		
		void
		CombineDomainsWith (const DSeriesT& rhs) const
		{
			domain_type domainMin = std::min(domainValues_.first(), rhs.domainValues_.first());
			domain_type domainMax = std::max(domainValues_.last(), rhs.domainValues_.last());
			
			//	Get the sampling rate
			samplingRate_ = std::max(GetSamplingRate(), rhs.GetSamplingRate());
			
			std::size_t numSamples ((domainMax - domainMin) / samplingRate_ + 1);
			
			domainValues_.resize(numSamples);
			
			for (std::size_t i = 0; i != numSamples; ++i) {
				domainValues_[i] = domain_type(i) / samplingRate_;
			}
		}
		
		/*
		 vR& GetHilbertTransform(void) const {
		 vR* result = new vR(codomainValues_.size());
		 
		 //	There are two ways to attempt this:
		 //		1)	Doing a straightforward convolution, or
		 //		2)	Using the convolution theorem to transform, multiply,
		 //			then inverse-transform.
		 
		 //
		 for (unsigned n = 0; n < codomainValues_.size(); ++n) {
		 Range sum();
		 
		 for (unsigned m = 1; m < codomainValues_.size(); m+=2) {
		 results += codomainValues_[n] * 2. / M_PI / m;
		 }
		 
		 results
		 }
		 return *result;
		 }
		 */
		
		domain_type
		GetSamplingRate (void) const
		{
			return samplingRate_;
		}
		
	
	//private:
	};
	
	template< template<typename /*T1*/> class ContainerDomainT
		, template<typename /*T2*/> class ContainerCodomainT = ContainerDomainT
		, typename T1 = double
		, typename T2 = T1
		, typename InterpT = FakeInterpolator<double> //CubicInterpFunc<T1, T2>
		>
	//template<typename Domain, typename Range, typename AllocD, typename AllocR>
	
	bool operator== (const DiscreteSeries<ContainerDomainT, ContainerCodomainT, T1, T2, InterpT>& lhs, const DiscreteSeries<ContainerDomainT, ContainerCodomainT, T1, T2, InterpT>& rhs);
	
	//friend
	template< template<typename /*T1*/> class ContainerDomainT
		, template<typename /*T2*/> class ContainerCodomainT = ContainerDomainT
		, typename T1 = double
		, typename T2 = T1
		, typename InterpT = FakeInterpolator<double> //CubicInterpFunc<T1, T2>
		>

//	template<typename Domain, typename Range, typename AllocD, typename AllocR>
	std::ostream& operator<< (std::ostream& os, const DiscreteSeries<ContainerDomainT, ContainerCodomainT, T1, T2, InterpT>& rhs) {
		
		os << "# Bin,\tTime" << ",\tValue" << std::endl;
		
		for (unsigned index (0); index < rhs.size(); ++index) {
			os << index << ",\t" << rhs.DiracComb.at(index) << ",\t" << rhs.codomainValues_.at(index) << std::endl;
		}
		
		return os;
	}
	
}


#ifdef DISCRETESERIES_USE_BOOST_SERIALIZE
	/*
	namespace boost {
		namespace serialization {
			template <typename , >
			struct version < PS::DiscreteSeries<T,U> >
			{

				typedef mpl::int_<1> type;
				typedef mpl::integral_c_tag tag;
				BOOST_STATIC_CONSTANT(unsigned int, value = version::type::value);
			};
		}
	}
	*/
	//BOOST_CLASS_VERSION(PS::DiscreteSeries,0)
#endif


#endif
