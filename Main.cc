#include <iostream>
#include <iomanip>
#include <boost/numeric/odeint.hpp>
#include "ShallowWater.hh"

using namespace std;
using namespace boost::numeric::odeint;

int main(){
	/* system instance
     * */
	//ColumnHeating sys("control.in");
	ShallowWater sys("control.in"); 
    sys.set_boundary_conditions();
	sys.init_variables();
    //sys.set_NH3_mixr(4.E-4);
	cout << sys
		<< setw(8) << left << "steps:" 
		<< setw(15) << left << "Model Time(s):"
		<< setw(15) << left << "Real Time(s):" << endl;
	/* stepper 
     * */
	runge_kutta4<State> stepper;
    //adams_bashforth<4, State> stepper;
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

