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


#include <DiscreteSeries.hpp>
#include <gtest/gtest.h>


using namespace PS;


namespace {



class DiscreteSeriesTest : public ::testing::Test {
	protected:

	DiscreteSeriesTest()
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

	}

	virtual
	void
	TearDown()
	{

	}


};


TEST_F(DiscreteSeriesTest, FillConstructor){}


TEST_F(DiscreteSeriesTest, FillWithRateConstructor){}


TEST_F(DiscreteSeriesTest, ContainersCopyConstructor){}


TEST_F(DiscreteSeriesTest, RangeCopyConstructor){}


TEST_F(DiscreteSeriesTest, CopyConstructor){}



TEST_F(DiscreteSeriesTest, AdditionAssignOperator){}


TEST_F(DiscreteSeriesTest, SubtractionAssignOperator){}


TEST_F(DiscreteSeriesTest, MultiplicationAssignOperator){}


TEST_F(DiscreteSeriesTest, DivisionAssignOperator){}


TEST_F(DiscreteSeriesTest, HeteroAdditionAssignOperator){}


TEST_F(DiscreteSeriesTest, HeteroSubtractionAssignOperator){}


TEST_F(DiscreteSeriesTest, HeteroMultiplicationAssignOperator){}


TEST_F(DiscreteSeriesTest, HeteroDivisionAssignOperator){}

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

