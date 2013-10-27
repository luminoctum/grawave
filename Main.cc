#include <iostream>
#include <iomanip>
#include <boost/numeric/odeint.hpp>
#include "Cartesian2d.hh"

using namespace std;
using namespace boost::numeric::odeint;

int main(){
	/* system instance
     * */
	Cartesian2d sys("control.in");
	//ShallowWater sys("control.in"); 
    sys.set_boundary_conditions();
	sys.init_variables();
	cout << sys
		<< setw(8) << left << "steps:" 
		<< setw(15) << left << "Model Time(s):"
		<< setw(15) << left << "Real Time(s):" << endl;
	/* stepper 
     * */
	runge_kutta4<StateType> stepper;
    //adams_bashforth<4, StateType> stepper;
	float t;
	for (t = sys.start(); t < sys.end(); t += sys.step()){
		sys.observe(t);
        // somehow, do_step do not update other variables in sys
        // and you should use sys.var explicitly
		stepper.do_step(sys, sys.var, t, sys.step());
		sys.update(t);
		//sys.debug();
	}
	sys.observe(t);
}

