//
//
//		Set:
//
//	$ export GTEST_HOME=gtest-1.7.0/
//	$ export LD_LIBRARY_PATH=$GTEST_HOME/lib:$LD_LIBRARYPATH
//
//		Compile with:
//
//	g++ -std=c++11 -I. -I $GTEST_HOME/include -L $GTEST_HOME/lib -lgtest -lgtest_main -lpthread test_DiscreteSeries.cpp -o test_DiscreteSeries
//
//
//

#include <cmath>
#include <vector>
#include <algorithm>


#include <fstream>

#include <DiscreteSeries.hpp>
#include <gtest/gtest.h>


#include <GSLInterpolator.hpp>

using namespace PS;


//	Modified EXPECT_ITERABLE_BASE so that it doesn't need type specifiers
#define EXPECT_ITERABLE_BASE_NOTYPE( PREDICATE, ref, target) \
    { \
	typedef decltype(ref) REFTYPE; \
	typedef decltype(target) TARTYPE; \
    const REFTYPE& ref_(ref); \
    const TARTYPE& target_(target); \
    auto /*typename REFTYPE::const_iterator*/ refIter = ref_.begin(); \
    auto /*typename TARTYPE::const_iterator*/ tarIter = target_.begin(); \
    unsigned int i = 0; \
    while(refIter != ref_.end()) { \
        if ( tarIter == target_.end() ) { \
            ADD_FAILURE() << #target " has a smaller length than " #ref ; \
            break; \
        } \
        PREDICATE(* refIter, * tarIter) \
            << "Containers " #ref  " (refIter) and " #target " (tarIter)" \
               " differ at index " << i; \
        ++refIter; ++tarIter; ++i; \
    } \
    EXPECT_TRUE( tarIter == target_.end() ) \
        << #ref " has a smaller length than " #target ; \
    }

#define EXPECT_ITERABLE_EQ_NOTYPE( ref, target) \
	EXPECT_ITERABLE_BASE_NOTYPE(EXPECT_EQ, ref, target )

#define EXPECT_ITERABLE_FLOAT_EQ_NOTYPE( ref, target) \
	EXPECT_ITERABLE_BASE_NOTYPE(EXPECT_FLOAT_EQ, ref, target )
#define EXPECT_ITERABLE_DOUBLE_EQ_NOTYPE( ref, target) \
	EXPECT_ITERABLE_BASE_NOTYPE(EXPECT_DOUBLE_EQ, ref, target )


namespace {



class DiscreteSeriesTest : public ::testing::Test {
	protected:

	DiscreteSeriesTest()
		: xs (length_)
		, ys (length_)
		, quarterSteps (length_)
	{

	}

	virtual
	~DiscreteSeriesTest()
	{

	}

	virtual
	void
	SetUp()
	{

		for (unsigned idx (0); idx < xs.size(); ++idx)
			xs.at(idx);

		
		for (unsigned idx (0); idx < ys.size(); ++idx)
			ys.at(idx) = sin (double(idx)/4.);

		
		for (unsigned idx (0); idx < quarterSteps.size(); ++idx)
			quarterSteps.at(idx) = double(idx) / 4.;
		
	}

	virtual
	void
	TearDown()
	{

	}

	public:
	
	const static unsigned length_ = 512;
	std::vector<double> xs;//(length_);
	std::vector<double> ys;//(length_);

	std::vector<double> quarterSteps;//(length_);

};

typedef DiscreteSeries< double, double, std::vector<double> >	DS_T;


TEST_F(DiscreteSeriesTest, FillConstructor)
{
	DS_T myDS (length_);

	
	//DiscreteSeries<> myDS2 (length_);
}


TEST_F(DiscreteSeriesTest, FillWithRateConstructor)
{


	DS_T myDS (length_, 1./4.);

}


TEST_F(DiscreteSeriesTest, NaiveLinear)
{
	DiscreteSeries<double,double,std::vector<double>,std::vector<double>,NaiveLinearInterpolator<std::vector<double> > >	myDS (ys, 4.);

	//EXPECT_DOUBLE_EQ(ys.at(100), myDS.at(100));

	//EXPECT_DOUBLE_EQ(ys.at(100), myDS(quarterSteps.at(100)));
	//EXPECT_DOUBLE_EQ(myDS(quarterSteps[100]), 0.);

	EXPECT_ITERABLE_DOUBLE_EQ_NOTYPE(ys, myDS);

	EXPECT_ITERABLE_DOUBLE_EQ_NOTYPE(quarterSteps, myDS.get_domain());
}

TEST_F(DiscreteSeriesTest, GSL_linear)
{
	DiscreteSeries<double,double,std::vector<double>,std::vector<double>,::GSLInterpolator_linear<std::vector<double> > >	myDS (ys, 4.);

	EXPECT_ITERABLE_DOUBLE_EQ_NOTYPE(ys, myDS);
	EXPECT_ITERABLE_DOUBLE_EQ_NOTYPE(quarterSteps, myDS.get_domain());
	//EXPECT_DOUBLE_EQ(ys.at(100), myDS.at(100));

	//EXPECT_DOUBLE_EQ(ys.at(100), myDS(quarterSteps.at(100)));
	//EXPECT_DOUBLE_EQ(myDS(quarterSteps[100]), 0.);
}

TEST_F(DiscreteSeriesTest, GSL_cspline)
{
	DiscreteSeries<double,double,std::vector<double>,std::vector<double>,::GSLInterpolator_cspline<std::vector<double> > >	myDS (ys, 4.);

	EXPECT_ITERABLE_DOUBLE_EQ_NOTYPE(ys, myDS);
	EXPECT_ITERABLE_DOUBLE_EQ_NOTYPE(quarterSteps, myDS.get_domain());

/*
	EXPECT_DOUBLE_EQ(ys.at(100), myDS.at(100));

	EXPECT_DOUBLE_EQ(ys.at(100), myDS(quarterSteps.at(100)));
	//EXPECT_DOUBLE_EQ(myDS(quarterSteps[100]), 0.);

	std::cout << myDS.get_domain().back() << std::endl;

	std::ofstream ofs ("cspline_test.dat");
	for (std::size_t i(0); i < myDS.size(); ++i) {
		ofs << myDS.get_domain().at(i) <<"\t"<< myDS.get_codomain().at(i) << "\t" 
		
			<< myDS(myDS.get_domain().at(i)) 
			//<< myDS(myDS.get_domain().at(i)+(2.2 * double(i >= myDS.size() - 100))) 
		
			<< std::endl;
	}

	

	ofs.close();
	
	std::ofstream ofs2 ("cspline_reference.dat");
	for (std::size_t i(0); i < myDS.size(); ++i) {
		ofs2 << quarterSteps.at(i) <<"\t"<< ys.at(i) << std::endl;
	}
	ofs2.close();
*/
}

TEST_F(DiscreteSeriesTest, GSL_polynomial)
{
	DiscreteSeries<double,double,std::vector<double>,std::vector<double>,::GSLInterpolator_polynomial<std::vector<double> > >	myDS (ys, 4.);

	EXPECT_ITERABLE_DOUBLE_EQ_NOTYPE(ys, myDS);
	EXPECT_ITERABLE_DOUBLE_EQ_NOTYPE(quarterSteps, myDS.get_domain());

	//EXPECT_DOUBLE_EQ(ys.at(100), myDS.at(100));

	//EXPECT_DOUBLE_EQ(ys.at(100), myDS(quarterSteps.at(100)));
	//EXPECT_DOUBLE_EQ(myDS(quarterSteps[100]), 0.);
}

TEST_F(DiscreteSeriesTest, ContainersCopyConstructor)
{
	ADD_FAILURE() << "Not Yet Implemented!";	
	
	//DS_T myDS (xs, ys);
	//DiscreteSeries<vecD> myDS (xs, ys);
}


TEST_F(DiscreteSeriesTest, ZipIterator)
{
	ADD_FAILURE() << "Not Yet Implemented!";

	DS_T myDS (xs, ys);

	//std::copy(myDS.begin_zip(), myDS.begin_zip() + 32, std::ostream_iterator<const typename decltype(myDS)::zipped_iterator::value_type>(std::cout,"\n"));

}

TEST_F(DiscreteSeriesTest, RangeCopyConstructor)
{

	DS_T myDS (ys, 4.);
	//DiscreteSeries<vecD> myDS (ys, 1./4.);

	EXPECT_EQ(myDS.size(), ys.size());

	for (unsigned idx (0); idx < myDS.size(); ++idx)
	{
		EXPECT_FLOAT_EQ(quarterSteps.at(idx), myDS.get_domain().at(idx));
		EXPECT_FLOAT_EQ(ys.at(idx), myDS.at(idx));
	}


}


TEST_F(DiscreteSeriesTest, CopyConstructor)
{
	DS_T myDS (xs, ys);
	//DiscreteSeries<vecD> myDS (xs, ys);

	
	DS_T myDScopy (myDS);
	//DiscreteSeries<vecD> myDScopy (myDS);

	
	EXPECT_EQ(myDS.size(), myDScopy.size());

	for (unsigned idx (0); idx < myDS.size(); ++idx)
	{
		EXPECT_FLOAT_EQ(xs.at(idx), myDS.get_domain().at(idx));
		EXPECT_FLOAT_EQ(myDS.at(idx), myDScopy.at(idx));
		EXPECT_FLOAT_EQ(xs.at(idx), myDScopy.get_domain().at(idx));
	}
}



TEST_F(DiscreteSeriesTest, AdditionAssignOperator)
{ADD_FAILURE() << "Not Yet Implemented!";}


TEST_F(DiscreteSeriesTest, SubtractionAssignOperator)
{ADD_FAILURE() << "Not Yet Implemented!";}


TEST_F(DiscreteSeriesTest, MultiplicationAssignOperator)
{ADD_FAILURE() << "Not Yet Implemented!";}


TEST_F(DiscreteSeriesTest, DivisionAssignOperator)
{ADD_FAILURE() << "Not Yet Implemented!";}


TEST_F(DiscreteSeriesTest, HeteroAdditionAssignOperator)
{ADD_FAILURE() << "Not Yet Implemented!";}


TEST_F(DiscreteSeriesTest, HeteroSubtractionAssignOperator)
{ADD_FAILURE() << "Not Yet Implemented!";}


TEST_F(DiscreteSeriesTest, HeteroMultiplicationAssignOperator)
{ADD_FAILURE() << "Not Yet Implemented!";}


TEST_F(DiscreteSeriesTest, HeteroDivisionAssignOperator)
{ADD_FAILURE() << "Not Yet Implemented!";}

/*
TEST_F(DiscreteSeriesTest, FillConstructor){}


TEST_F(DiscreteSeriesTest, FillConstructor){}


TEST_F(DiscreteSeriesTest, FillConstructor){}


TEST_F(DiscreteSeriesTest, FillConstructor){}
*/

}	//	namespace

int
main (int argc, char** argv)
{
	::testing::InitGoogleTest(&argc,argv);
	return RUN_ALL_TESTS();
}

