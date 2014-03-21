

#include <cmath>

#include <iostream>
#include <exception>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

/*

Interpolation Type: gsl_interp_linear
Linear interpolation. This interpolation method does not require any additional memory.

Interpolation Type: gsl_interp_polynomial
Polynomial interpolation. This method should only be used for interpolating small numbers of points because polynomial interpolation introduces large oscillations, even for well-behaved datasets. The number of terms in the interpolating polynomial is equal to the number of points.

Interpolation Type: gsl_interp_cspline
Cubic spline with natural boundary conditions. The resulting curve is piecewise cubic on each interval, with matching first and second derivatives at the supplied data-points. The second derivative is chosen to be zero at the first point and last point.

Interpolation Type: gsl_interp_cspline_periodic
Cubic spline with periodic boundary conditions. The resulting curve is piecewise cubic on each interval, with matching first and second derivatives at the supplied data-points. The derivatives at the first and last points are also matched. Note that the last point in the data must have the same y-value as the first point, otherwise the resulting periodic interpolation will have a discontinuity at the boundary.

Interpolation Type: gsl_interp_akima
Non-rounded Akima spline with natural boundary conditions. This method uses the non-rounded corner algorithm of Wodicka.

Interpolation Type: gsl_interp_akima_periodic
Non-rounded Akima spline with periodic boundary conditions. This method uses the non-rounded corner algorithm of Wodicka.











 */

//extern gsl_interp_linear


struct interp_linear_struct {};
struct interp_polynomial_struct {};
struct interp_cspline_struct {};
struct interp_cspline_periodic_struct {};
struct interp_akima_struct {};
struct interp_akima_periodic_struct {};

enum interp_type_enum {
	type_linear, type_polynomial, type_cspline, type_cspline_periodic, type_akima, type_akima_periodic};

template <typename T>
struct GSLInterpType {
	static const interp_type_enum interp_type;
};

template <>
struct GSLInterpType<interp_linear_struct>
{ static const interp_type_enum interp_type = type_linear; };

template <>
struct GSLInterpType<interp_polynomial_struct>
{ static const interp_type_enum interp_type = type_polynomial; };

template <>
struct GSLInterpType<interp_cspline_struct>
{ static const interp_type_enum interp_type = type_cspline; };

template <>
struct GSLInterpType<interp_cspline_periodic_struct>
{ static const interp_type_enum interp_type = type_cspline_periodic; };

template <>
struct GSLInterpType<interp_akima_struct>
{ static const interp_type_enum interp_type = type_akima; };

template <>
struct GSLInterpType<interp_akima_periodic_struct>
{ static const interp_type_enum interp_type = type_akima_periodic; };




template< typename RandomAccessRange1
		, typename RandomAccessRange2 = RandomAccessRange1
		, typename InterpTypeT = interp_cspline_struct
		>
class GSLInterpolator {
	private:
	typedef typename std::iterator_traits< typename RandomAccessRange1::iterator >::value_type	T1;
	typedef typename std::iterator_traits< typename RandomAccessRange2::iterator >::value_type	T2;
		

	gsl_interp* 	  interp_;
	gsl_interp_accel* accel_;

	const T1*	xFirst_;
	const T2* yFirst_;

	public:

	//template<typename RandomAccessIterator1, typename RandomAccessIterator2>
	/*
	GSLInterpolator (RandomAccessIterator1 xFirst, RandomAccessIterator1 xLast, RandomAccessIterator2 yFirst, const gsl_interp_type* interpType)
		: xFirst_(&(*xFirst))
		, yFirst_(&(*yFirst))
	{
		accel_  = gsl_interp_accel_alloc();
		interp_ = gsl_interp_alloc(interpType, std::distance(xFirst, xLast));

		gsl_interp_init(interp_, xFirst_, yFirst_, std::distance(xFirst, xLast));


	}
	*/
	
	GSLInterpolator (const RandomAccessRange1& rng1, const RandomAccessRange2& rng2)
		: xFirst_ (&(*boost::begin(rng1)))
		, yFirst_ (&(*boost::begin(rng2)))
		, accel_ (gsl_interp_accel_alloc())
		//, interp_ (gsl_interp_alloc(InterpTypeT::interp_type, boost::distance(rng1)))
	{
		switch (GSLInterpType<InterpTypeT>::interp_type) {
			case type_linear:
				if (boost::distance(rng1) < gsl_interp_type_min_size(gsl_interp_linear))
					throw std::length_error("GSLInterpolator: Array has fewer elements than is required by this interpolation method!");
				interp_ = gsl_interp_alloc(gsl_interp_linear, boost::distance(rng1));
				break;

			case type_polynomial:
				if (boost::distance(rng1) < gsl_interp_type_min_size(gsl_interp_polynomial))
					throw std::length_error("GSLInterpolator: Array has fewer elements than is required by this interpolation method!");
				interp_ = gsl_interp_alloc(gsl_interp_polynomial, boost::distance(rng1));
				break;

			case type_cspline:
				if (boost::distance(rng1) < gsl_interp_type_min_size(gsl_interp_cspline))
					throw std::length_error("GSLInterpolator: Array has fewer elements than is required by this interpolation method!");
				interp_ = gsl_interp_alloc(gsl_interp_cspline, boost::distance(rng1));
				break;

			case type_cspline_periodic:
				if (boost::distance(rng1) < gsl_interp_type_min_size(gsl_interp_cspline_periodic))
					throw std::length_error("GSLInterpolator: Array has fewer elements than is required by this interpolation method!");
				interp_ = gsl_interp_alloc(gsl_interp_cspline_periodic, boost::distance(rng1));
				break;

			case type_akima:
				if (boost::distance(rng1) < gsl_interp_type_min_size(gsl_interp_akima))
					throw std::length_error("GSLInterpolator: Array has fewer elements than is required by this interpolation method!");
				interp_ = gsl_interp_alloc(gsl_interp_akima, boost::distance(rng1));
				break;

			case type_akima_periodic:
				if (boost::distance(rng1) < gsl_interp_type_min_size(gsl_interp_akima_periodic))
					throw std::length_error("GSLInterpolator: Array has fewer elements than is required by this interpolation method!");
				interp_ = gsl_interp_alloc(gsl_interp_akima_periodic, boost::distance(rng1));
				break;

			default:
				if (boost::distance(rng1) < gsl_interp_type_min_size(gsl_interp_cspline))
					throw std::length_error("GSLInterpolator: Array has fewer elements than is required by this interpolation method!");
				interp_ = gsl_interp_alloc(gsl_interp_cspline, boost::distance(rng1));
				break;
		}
		

		
		gsl_interp_init(interp_, xFirst_, yFirst_, boost::distance(rng1));
	}
	
/*
	GSLInterpolator (const RandomAccessRange1& rng1, const RandomAccessRange2& rng2, const gsl_interp_type* interpType = gsl_interp_cspline)
		: xFirst_ (&(*boost::begin(rng1)))
		, yFirst_ (&(*boost::begin(rng2)))
		, accel_ (gsl_interp_accel_alloc())
		, interp_ (gsl_interp_alloc(interpType, boost::distance(rng1)))
	{
		gsl_interp_init(interp_, xFirst_, yFirst_, boost::distance(rng1));
	}
*/

	~GSLInterpolator (void)
	{
		gsl_interp_free(interp_);
		gsl_interp_accel_free(accel_);
	}


	T2 eval (const T1 input)
	{
		/*
		T2* result;
		int erVal =gsl_interp_eval_e(interp_, xFirst_, yFirst_, input, accel_, result);
		
		if (erVal == GSL_EDOM)
			std::cout << "EDOM: " << input << "\t("<< *xFirst_<<", "<< *xFirst_<<")"<< std::endl;


		//if (erVal == GSL_NAN)
		//	std::cout << "NAN" << std::endl;

		return *result;
		*/


		return gsl_interp_eval(interp_, xFirst_, yFirst_, input, accel_);
	}

	T2 operator() (const T1 input)
	{ return eval(input); }

};


//	Now alias the templates using the different GSL interpolation types:


template <typename RandomAccessRange1, typename RandomAccessRange2 = RandomAccessRange1>
using GSLInterpolator_linear = GSLInterpolator<RandomAccessRange1, RandomAccessRange2, interp_linear_struct>;

template <typename RandomAccessRange1, typename RandomAccessRange2 = RandomAccessRange1>
using GSLInterpolator_polynomial = GSLInterpolator<RandomAccessRange1, RandomAccessRange2, interp_polynomial_struct>;

template <typename RandomAccessRange1, typename RandomAccessRange2 = RandomAccessRange1>
using GSLInterpolator_cspline = GSLInterpolator<RandomAccessRange1, RandomAccessRange2, interp_cspline_struct>;

template <typename RandomAccessRange1, typename RandomAccessRange2 = RandomAccessRange1>
using GSLInterpolator_cspline_periodic = GSLInterpolator<RandomAccessRange1, RandomAccessRange2, interp_cspline_periodic_struct>;

template <typename RandomAccessRange1, typename RandomAccessRange2 = RandomAccessRange1>
using GSLInterpolator_akima = GSLInterpolator<RandomAccessRange1, RandomAccessRange2, interp_akima_struct>;

template <typename RandomAccessRange1, typename RandomAccessRange2 = RandomAccessRange1>
using GSLInterpolator_akima_periodic = GSLInterpolator<RandomAccessRange1, RandomAccessRange2, interp_akima_periodic_struct>;


