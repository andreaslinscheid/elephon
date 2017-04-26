#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Input_test
#include <boost/test/unit_test.hpp>
#include <string>
#include "IOMethods/Input.h"

BOOST_AUTO_TEST_CASE( Default_Args )
{
	std::string filename("test.dat");
	int argc = 1;
	char const * argv[argc+1];
	argv[0] = "unused binary name";
	argv[1] = "../IOMethods/test_input_file.dat";
	elephon::IOMethods::Input input(argc,argv);
	auto sc = input.get_super_cell_size();
	BOOST_CHECK( sc.size() == 3 );
	BOOST_CHECK( sc[0] == 2 );
	BOOST_CHECK( sc[1] == 2 );
	BOOST_CHECK( sc[2] == 1 );

	BOOST_CHECK( input.get_numFS() == 2000 );
}
