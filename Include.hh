#ifndef INCLUDE 
#define INCLUDE
typedef Eigen::ArrayXXf Grid;
typedef std::vector<Grid> StateType;
/* define error messages */
#define ALARM(message)  std::cout << message << std::endl 
#define ASSERT_FILE_NOT_FOUND(AFILE) \
        std::cerr << "***** ERROR *****" << std::endl; \
        std::cerr << "can not open file \"" << AFILE << "\" ..."<< std::endl; \
        std::cerr << "*****************" << std::endl; \
        assert(0) 
#define ASSERT_VARIABLE_OUT_OF_RANGE(AVAR) \
        std::cerr << "***** ERROR *****" << std::endl; \
        std::cerr << "variable \"" << AVAR << "\" out of range ..."<< std::endl; \
        std::cerr << "*****************" << std::endl; \
        assert(0) 
#define ASSERT_DIMENSION_MISMATCH(AVAR1, AVAR2) \
        std::cerr << "***** ERROR *****" << std::endl; \
        std::cerr << "dimension mismatch for variable \"" << AVAR1 << "\" and \"" << AVAR2 << "\" ..."<< std::endl; \
        std::cerr << "*****************" << std::endl; \
        assert(0) 

/* define constants */
#define ZERO1(n) Eigen::ArrayXf::Zero(n)
#define ZERO2(n,m) Eigen::ArrayXXf::Zero(n,m)
#endif
