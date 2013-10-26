#ifndef INCLUDE 
#define INCLUDE
typedef Eigen::ArrayXXf Grid;
typedef std::vector<Grid> StateType;
#define _alarm(message)  std::cout << message << std::endl 
#define ERROR_1  assert(0); _alarm("axis out of range")
#endif
