#ifndef BOUNDARY
#define BOUNDARY
#include "Include.hh"

enum Type{Neumann, Dirichlet, Periodic, Undefined};
struct Edge{
	Edge(){ type = Undefined;}

	Edge& operator<< (const Type &_type){ type = _type; return *this;}
	Edge& operator| (const Grid &_value){ value = _value; return *this;}

	friend std::ostream& operator<< (std::ostream &os, const Edge &other){	
		os << "\tType : ";
		switch (other.type){
			case Neumann:
				os << "Neumann  \t";
				break;
			case Dirichlet:
				os << "Dirichlet\t";
				break;
			case Periodic:
				os << "Periodic \t";
				break;
			case Undefined:
				os << "Undefined\t";
				break;
		}
		os << " | ";
		if (other.value.size() > 4){
			os << "[" << other.value(0) << ", "
				<< other.value(1) << ", "
				<< "... "
				<< other.value(other.value.size()-2) << ", "
				<< other.value(other.value.size()-1) << "] ";
		} else{
            if (other.value.rows() == 1) os << other.value;
            else os << other.value.transpose();
		}
		return os;
	}

	Type type;
	Grid value;
};

struct Boundary{
	friend std::ostream& operator<< (std::ostream &os, const Boundary &other){	
		os << "Left" << other.left << std::endl;
		os << "Right"<< other.right << std::endl;
		os << "Bottom"<< other.bottom << std::endl;
		os << "Top" << other.top;
		return os;
	}

	Edge left;
	Edge right;
	Edge bottom;
	Edge top;
};
#endif
