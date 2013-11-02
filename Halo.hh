#ifndef HALO
#define HALO
#include "Include.hh"

enum BoundaryType{Neumann, Dirichlet, Periodic, Undefined};

struct Boundary{
	BoundaryType type;
	Grid value;

	Boundary(){ type = Undefined;}

	Boundary& operator<< (const BoundaryType &_type){ type = _type; return *this;}
	Boundary& operator| (const Grid &_value){ value = _value; return *this;}

	friend std::ostream& operator<< (std::ostream &os, const Boundary &other){	
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
};

class Halo{
public:
    Grid left, right, bottom, top;

	friend std::ostream& operator<< (std::ostream &os, const Halo &other){	
		os << "Left" << other.l << std::endl;
		os << "Right"<< other.r << std::endl;
		os << "Bottom"<< other.b << std::endl;
		os << "Top" << other.t;
		return os;
	}
    void set_left_right(BoundaryType type, Grid value){r = l << type | value;}
    void set_left_right(BoundaryType type){r = l << type;}
    void set_bottom_top(BoundaryType type, Grid value){b = t << type | value;}
    void set_bottom_top(BoundaryType type){b = t << type;}
    void set_left(BoundaryType type, Grid value){ l << type | value;}
    void set_left(BoundaryType type){ l << type;}
    void set_right(BoundaryType type, Grid value){ r << type | value;}
    void set_right(BoundaryType type){ r << type;}
    void set_bottom(BoundaryType type, Grid value){ b << type | value;}
    void set_bottom(BoundaryType type){ b << type;}
    void set_top(BoundaryType type, Grid value){ t << type | value;}
    void set_top(BoundaryType type){ t << type;}
    void set_periodic(){r = l = b = t << Periodic;}

    void update(const Grid &var){
        switch (l.type){
            case Neumann:
                left = var.row(0) - l.value;
                break;
            case Dirichlet:
                left = l.value;
                break;
            case Periodic:
                left = var.row(var.rows() - 1);
                break;
            case Undefined:
                left = ZERO2(1, var.cols());
                break;
        }
        switch (r.type){
            case Neumann:
                right = var.row(var.rows() - 1) - r.value;
                break;
            case Dirichlet:
                right = r.value;
                break;
            case Periodic:
                right = var.row(0);
                break;
            case Undefined:
                right = ZERO2(1, var.cols());
                break;
        }
        switch (b.type){
            case Neumann:
                bottom = var.col(0) - b.value;
                break;
            case Dirichlet:
                bottom = b.value;
                break;
            case Periodic:
                bottom = var.col(var.cols() - 1);
                break;
            case Undefined:
                bottom = ZERO2(var.rows(), 1);
                break;
        }
        switch (t.type){
            case Neumann:
                top = var.col(var.cols() - 1) - b.value;
                break;
            case Dirichlet:
                top = t.value;
                break;
            case Periodic:
                top = var.col(0);
                break;
            case Undefined:
                top = ZERO2(var.rows(), 1);
                break;
        }
    }
protected:
    Boundary l, r, b, t;
};
#endif
