/*********************************************
 * OPL 12.7.0.0 Model
 * Author: Administrator
 * Creation Date: 2018��11��16�� at ����12:54:01
 *********************************************/
 
 int n=2;
float objectiveForXEqualsStart=300;
float breakpoint[1..n]=[100,200];
float slope[1..n+1]=[1,2,-3];
dvar int x;

maximize piecewise(i in 1..n) 
{slope[i] -> breakpoint[i]; -3}(0,objectiveForXEqualsStart) x;


subject to{
	x == 250;
}


